"""
Heat bundles of fieldlines together
"""
import warnings

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

from synthesizAR.util import heeq_to_hcc_coord


class BundleHeatingModel(object):
    """
    Heating model for heating many strands simultaneously
    """
    
    def __init__(self, heating_options, field, storm_duration, num_bins=50, footpoints=None):
        self.field = field
        self.storm_duration = storm_duration
        self.heating_options = heating_options
        self.num_bins = num_bins
        if footpoints is None:
            self.footpoints = self.get_loop_footpoints()
        else:
            self.footpoints = footpoints
        
    def calculate_event_properties(self, loop):
        """
        Get the event start times and rates for each event on the loop
        """
        try:
            bundle = [b for b in self.bundles if loop in b][0]
        except IndexError:
            raise IndexError(f'{loop.name} not found in any bundles')
        start_times = bundle.storm_start_times + bundle.get_offset(loop)
        rates = (self.stress**2) * (bundle.storm_fluxes / bundle.number_loops) / loop.full_length.to(u.cm) / (self.heating_options['duration'] / 2.) 
        
        return {'magnitude': rates.value,
                'rise_start': start_times.to(u.s).value,
                'rise_end': (start_times + self.heating_options['duration']/2.).to(u.s).value,
                'decay_start': (start_times + self.heating_options['duration']/2).to(u.s).value,
                'decay_end': (start_times + self.heating_options['duration']).to(u.s).value,}
    
    def constrain_number_storms(self, total_time, number_storms=100, flux_constraint=1e7*u.erg/u.cm**2/u.s,
                                error_tolerance=1e-2, safety_decrease=0.5, safety_increase=2., max_tries=100):
        """
        For a given constraint on the total AR flux, find the total number of storms
        over all bundles which satisfies this constraint
        
        This needs to be run before configuring any of the simulations
        """
        # Choose starting value for Nstorms
        # Select start times from uniform distribution
        # Assign each time to a bundle
        # Compute total flux from all bundles and all storms
        # Compare to constraint and correct
        # Continue if not satisfied
        # Else stop and set the starting times of storms per bundle
        # Warning if not met in finite number of tries
        error = 1e300
        tries = 0
        total_flux_constraint = flux_constraint * total_time
        while error > error_tolerance and tries < max_tries:
            # Get storm times
            storm_times = np.random.uniform(low=0., high=total_time.to(u.s).value, size=number_storms)
            # Create bundles
            bundles = self.create_bundles(self.footpoints, self.num_bins)
            p_b = self.get_bundle_probabilities(bundles)
            # Assign each storm to a bundle based on the relative field strength
            storm_distribution = np.random.choice(len(bundles), size=number_storms, replace=True, p=p_b)
            total_flux = 0. * u.erg / (u.cm**2)
            for i,b in enumerate(bundles):
                storm_indices, = np.where(storm_distribution == i)
                try:
                    b.storm_start_times = np.sort(storm_times[storm_indices])
                    total_flux += b.total_flux
                except IndexError:
                    warnings.warn(f'No storms found for bundle {i}')
            total_flux *= self.heating_options['stress']**2
            phi = total_flux / total_flux_constraint
            error = np.fabs(1. - phi)
            number_storms = int(min(max(number_storms + (1. - phi)*number_storms, safety_decrease*number_storms), safety_increase*number_storms))
            tries += 1
            print(f'Error = {error}, phi = {phi}, n_storms = {number_storms}')
            
        if tries >= max_tries:
            warnings.warn('Exceeded maximum number of tries')
        self.bundles = bundles
    
    def get_loop_footpoints(self):
        """
        Make a list of the HPC coordinates of all loop footpoints
        """
        footpoints = np.zeros((len(self.field.loops), 2)) * u.arcsec
        for i,l in enumerate(self.field.loops):
            l_hpc = (heeq_to_hcc_coord(*l.coordinates.T, self.field.magnetogram.observer_coordinate)
                     .transform_to(self.field.magnetogram.coordinate_frame))
            footpoints[i,0] = l_hpc.Tx[0]
            footpoints[i,1] = l_hpc.Ty[0]
            
        return footpoints
    
    def create_bundles(self, footpoints, num_bins,):
        """
        Create bundles by dividing up the magnetogram into equally sized boxes and binning
        the loop footpoints
        """
        # Bin footpoints into regular grid
        hist, bin_edges_x, bin_edges_y = np.histogram2d(
            footpoints[:,0],footpoints[:,1], bins=(num_bins, num_bins),
            range=((self.field.magnetogram.bottom_left_coord.Tx.to(u.arcsec).value,
                    self.field.magnetogram.top_right_coord.Tx.to(u.arcsec).value),
                   (self.field.magnetogram.bottom_left_coord.Ty.to(u.arcsec).value,
                    self.field.magnetogram.top_right_coord.Ty.to(u.arcsec).value)))
        # Find index of cell that the footpoints fall into
        bundle_indices = np.stack([np.digitize(footpoints[:,0], bin_edges_x),
                                   np.digitize(footpoints[:,1], bin_edges_y)], axis=1)
        # Create bundles from non-empty cells
        bundles = []
        for i_x, i_y in np.stack(np.where(hist != 0), axis=1):
            loop_indices, = np.where(np.logical_and(
                bundle_indices[:,0]-1 == i_x, bundle_indices[:,1]-1 == i_y))
            loops = np.array(self.field.loops)[loop_indices].tolist()
            # Define bundle bounding box
            bottom_left_corner = SkyCoord(bin_edges_x[i_x]*u.arcsec, bin_edges_y[i_y]*u.arcsec,
                                          frame=self.field.magnetogram.coordinate_frame)
            top_right_corner = SkyCoord(bin_edges_x[i_x+1]*u.arcsec, bin_edges_y[i_y+1]*u.arcsec,
                                        frame=self.field.magnetogram.coordinate_frame)
            bundles.append(Bundle(loops, bottom_left_corner, top_right_corner, self.storm_duration))
            
        return bundles
    
    def get_bundle_probabilities(self, bundles):
        """
        Compute probability of a storm occuring on a given bundle
        """
        field_strengths = np.array([b.footpoint_field_strength.value for b in bundles])
        return field_strengths / field_strengths.sum()
    
    
class Bundle(object):
    
    def __init__(self, loops, bottom_left_corner, top_right_corner, storm_duration):
        self.loops = loops
        self.corners = (bottom_left_corner, top_right_corner)
        self.storm_duration = storm_duration
        self.storm_start_times = u.Quantity([], 's')
    
    def __contains__(self, loop):
        return loop in self.loops
    
    def __getitem__(self, key):
        return self.loops[key]
    
    @property
    def number_loops(self):
        return len(self.loops)
    
    @property
    def number_storms(self):
        """
        Total number of storms occuring on the bundle
        """
        return self.storm_start_times.shape[0]
        
    @property
    def total_flux(self):
        """
        Total flux into bundle over all storms. Does not include stress factor
        """
        total_flux = 0.
        for l in self:
            # Set apex where the radial distance is max
            s_apex = l.field_aligned_coordinate[np.argmax(np.sqrt((l.coordinates**2).sum()))].to(u.cm)
            # Use the "coronal" portion of the field (upper 50%) to compute the total storm flux
            delta_s = 0.25 * l.full_length.to(u.cm)
            i_corona, = np.where(np.logical_and(
                l.field_aligned_coordinate.to(u.cm) >= s_apex - delta_s,
                l.field_aligned_coordinate.to(u.cm) <= s_apex + delta_s))
            B_corona = np.average(l.field_strength[i_corona], weights=np.gradient(l.field_aligned_coordinate[i_corona]))
            total_flux += ((B_corona.to(u.Gauss))**2 * l.full_length.to(u.cm) / (8.*np.pi)).value
            
        return total_flux * self.number_storms * u.erg / u.cm**2
            
    def get_offset(self, loop):
        """
        Find the offset of the loop event time from the storm starting time
        """
        i = self.loops.index(loop)
        return self.storm_duration / len(self.loops) * i
    
    @property
    def storm_wait_times(self):
        """
        Delays between storms. The first storm on each bundle is relative to t=0
        """
        if self.number_storms == 0:
            return u.Quantity([], 's')
        else:
            return np.diff(np.insert(self.storm_start_times, 0, 0*u.s,))
        
    @property
    def storm_fluxes(self):
        """
        Flux for each storm. Proportional to delay prior to event
        """
        if self.number_storms == 0:
            return u.Quantity([], 'erg / (cm2 s)')
        else:
            return self.storm_wait_times / self.storm_wait_times.sum() * self.total_flux
        
    @property
    def footpoint_field_strength(self):
        """
        Average field strength of all loops at one footpoint
        """
        return u.Quantity([l.field_strength.to(u.Gauss).value[0] for l in self], 'Gauss').sum() / self.number_loops
        