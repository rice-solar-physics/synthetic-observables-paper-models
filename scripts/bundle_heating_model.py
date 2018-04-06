"""
Heat bundles of fieldlines together
"""
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

from synthesizAR.util import heeq_to_hcc_coord


class BundleHeatingModel(object):
    """
    Heating model for heating many strands simultaneously
    """
    
    def __init__(self, heating_options, field, storm_duration, num_bins=50,):
        self.field = field
        self.storm_duration = storm_duration
        self.heating_options = heating_options
        self.num_bins = num_bins
        # footpoints = self.get_loop_footpoints()
        # self.bundles = self.create_bundles(footpoints,num_bins)
        
    def calculate_event_properties(self, loop):
        """
        Get the event start times and rates for each event on the loop
        """
        # Find which bundle the loop is in
        
        # For each storm occuring on bundle b, the start time is the storm start time + the loop offset
        
        # The rate is the flux for the storm s times the stress parameter divided by the number of loops on the bundle, 
        # divided by the loop length, divided by half the duration of the event (for triangular pulses)
        
        return {'magnitude':rates,
                'rise_start':rise_start,
                'rise_end':rise_end,
                'decay_start':decay_start,
                'decay_end':decay_end,
               }
    
    def constrain_number_storms(self, total_time, number_storms=  flux_constraint=1e7*u.erg/u.cm**2/u.s):
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
        # Else stop and assign attributes appropriately
        # * total number of storms
        # * Starting times of storms per bundle
        # Warning if not met in finite number of tries
    
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
        for i_x, i_y for np.stack(np.where(hist != 0), axis=1):
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
        
    
class Bundle(object):
    
    def __init__(self, loops, bottom_left_corner, top_right_corner, storm_duration):
        self.loops = loops
        self.corners = (bottom_left_corner, top_right_corner)
        self.storm_duration = storm_duration
        self.storm_start_times = u.Quantity([], 's')
        
    def get_event_offset(self, loop):
        """
        Find the offset of the loop event time from the storm starting time
        """
        i = self.loops.index(loop)
        return self.storm_duration / len(self.loops) * i
    
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
            B_corona = np.mean(l.field_strength[i_corona], weights=np.gradient(l.field_aligned_coordinate[i_corona]))
            total_flux += ((B_corona.to(u.Gauss))**2 * l.full_length.to(u.cm) / (8.*np.pi)).value
            
        return total_flux * self.number_storms * u.erg / u.cm**2
    
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
        return u.Quantity([l.field_strength.to(u.Gauss).value[0] for l in self], 'Gauss') / self.number_loops
        
    