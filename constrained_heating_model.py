"""
Heating model that constrains total flux according to Withbroe and Noyes 
estimates of coronal flux.
"""
import numpy as np


class CustomHeatingModel(object):
    """
    Computes magnitudes and waiting times between successive events
    given a desired ratio between the cooling time and waiting period. The
    total energy output is also constrained by the estimates of Withbroe and
    Noyes (1978)
    """
    def __init__(self,heating_options):
        self.heating_options = heating_options
    
    def calculate_event_properties(self,loop):
        self.number_events = self._calculate_number_events(loop)
        rates = (self.power_law_distributions[loop.name]
                 /(self.heating_options['duration'] 
                   - 0.5*(self.heating_options['duration_rise'] + self.heating_options['duration_decay'])))
        delays = (self.base_config['total_time'] - self.number_events*self.heating_options['duration'])*rates/rates.sum()
        delays *= np.random.uniform(low=1. - 1./(self.number_events**2),high=1.)
        running_total = 0.0
        rise_start = np.empty(rates.shape)
        for i in range(self.number_events):
            running_total += delays[i]
            rise_start[i] = i*self.heating_options['duration'] + running_total
        rise_end = rise_start + self.heating_options['duration_rise']
        decay_start = rise_end
        decay_end = rise_start + self.heating_options['duration']
        return {'magnitude':rates,
                'rise_start':rise_start,
                'rise_end':rise_end,
                'decay_start':decay_start,
                'decay_end':decay_end
               }
    
    def _calculate_number_events(self, loop):
        mean_waiting_period = self.heating_options['frequency_parameter']*self.cooling_time(loop)
        num_events = int(np.floor((self.base_config['total_time'] + mean_waiting_period)/
                                  (self.heating_options['duration'] + mean_waiting_period))) 
        return num_events
    
    def power_law(self, a0, a1, alpha, x):
        """
        Truncated power-law distribution
        """
        return ((a1**(alpha + 1.) - a0**(alpha + 1.))*x + a0**(alpha + 1.))**(1./(alpha + 1.))
    
    def max_strand_energy(self,loop):
        return ((self.heating_options['stress_level']*loop.field_strength.value.max())**2)/8./np.pi
    
    def constrain_distribution(self, field, iterator_options):
        # Initialize quantities
        power_law_distributions = {}
        num_iterations = 0
        error = 1e300
        cross_sections = self.calculate_cross_sections(field)
        upper_bounds = np.array([self.max_strand_energy(loop) for loop in field.loops])
        lower_bounds = upper_bounds/100.
        
        # Iteratively adjust lower power-law bound
        while error > iterator_options['tolerance'] and num_iterations < iterator_options['max_iterations']:
            weights = np.empty(len(field.loops))
            # Calculate power-law distributions for all loops
            for i,loop in enumerate(field.loops):
                pl = self.power_law(lower_bounds[i],upper_bounds[i],self.heating_options['power_law_slope'],
                                    np.random.rand(self._calculate_number_events(loop)))
                power_law_distributions[loop.name] = pl
                weights[i] = pl.sum()*loop.full_length.value*cross_sections[i]
            
            # Update error terms
            phi = weights.sum()/cross_sections.sum()/iterator_options['total_ar_flux']/self.base_config['total_time']
            error = np.fabs(1. - 1./phi)
            # Update lower-bounds
            lower_bounds = np.minimum(np.maximum(lower_bounds + lower_bounds*(1. - phi),
                                                 iterator_options['delta0']*upper_bounds),
                                      iterator_options['delta1']*upper_bounds)
            num_iterations += 1
            
        if num_iterations == iterator_options['max_iterations']:
            warnings.warn('Max number of iterations reached with error {}'.format(error))
        
        self.power_law_distributions = power_law_distributions
    
    def cooling_time(self, loop):
        """
        Estimate loop cooling time.
        """
        half_length = loop.full_length.value/2.
        average_heating_rate_max = self.max_strand_energy(loop)/(self.heating_options['duration']/2.)#*u.erg/(u.cm**3)/u.s
        # set some constants
        alpha = -0.5
        chi = 6e-20#*(u.erg*(u.cm**3)/u.s*u.K**(0.5))
        kappa_0 = 1e-6#*(u.erg/u.cm/u.s*(u.K**(-7/2)))
        c1,c2,c3 = 2.0,0.9,0.6
        gamma = 5./3.
        # estimate max n0T0
        T0 = c2*(7.*half_length**2*average_heating_rate_max/2./kappa_0)**(2./7.)
        top_term = average_heating_rate_max - 2.*kappa_0*(T0**(3.5))/(7.*(c2**2.5)*c3*(half_length**2)*gamma)
        bottom_term = c1*chi*(T0**alpha)*(1. - c2/c3/gamma)
        n0 = np.sqrt(top_term/bottom_term)
        n0T0 = n0*T0
        # Cargill cooling expression
        term1 = (2. - alpha)/(1. - alpha)
        term2 = (kappa_0**(4. - 2.*alpha))*(chi**7)
        term3 = ((half_length)**(8. - 4.*alpha))/(n0T0**(3+2.*alpha))
        return term1*3.*const.k_B.cgs.value*(1/term2*term3)**(1/(11. - 2.*alpha))
    
    def calculate_cross_sections(self, field):
        """
        Estimate loop cross-sectional area
        """
        footpoints = np.array([loop.coordinates[0,:2].value for loop in field.loops])
        footpoint_dist,bin_edges = np.histogramdd(
            footpoints,
            bins = (field.clipped_hmi_map.dimensions.x.value,field.clipped_hmi_map.dimensions.y.value),
            range = (field._convert_angle_to_length(field.clipped_hmi_map.xrange).value,
                     field._convert_angle_to_length(field.clipped_hmi_map.yrange).value))
        area_per_pixel = (field._convert_angle_to_length(field.clipped_hmi_map.scale.x*1*u.pixel)
                          *field._convert_angle_to_length(field.clipped_hmi_map.scale.y*1*u.pixel)).value
        ix = np.digitize(footpoints[:,0],bin_edges[0]) - 1
        iy = np.digitize(footpoints[:,1],bin_edges[1]) - 1
        
        return (area_per_pixel/footpoint_dist)[ix,iy]