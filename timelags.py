"""
Class for computing timelags from AIA maps. Uses Dask to compute timelags in each pixel in parallel
"""
import os
import warnings

import h5py
import numpy as np
import dask.array
from sunpy.map import Map,GenericMap
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.utils.console import ProgressBar


class AIATimeLags(object):
    """
    Compute AIA timelag maps from many synthesized AIA observations
    
    Parameters
    ----------
    instr : synthesizAR.instruments.InstrumentSDOAIA
    fits_root_path : str
    """
    
    def __init__(self, instr, fits_root_path):
        self.instr = instr
        self.ff_format = os.path.join(fits_root_path, f'{instr.name}','{channel}','map_t{i_time:06d}.fits')
        delta_t = np.diff(instr.observing_time.value).cumsum()
        self.timelags = np.hstack([-delta_t[::-1],np.array([0]),delta_t])*self.instr.observing_time.unit
        self.cubes,self.meta_templates = self.load_data()
        
    def load_data(self, **kwargs):
        """
        Create data cubes over desired time range for all channels
        """
        cubes = {}
        meta_templates = {}
        with h5py.File(self.instr.counts_file,'r') as hf:
            instr_time = u.Quantity(hf['time'],hf['time'].attrs['units'])
        with ProgressBar(len(self.instr.channels)*len(self.instr.observing_time),
                         ipython_widget=kwargs.get('notebook',True)) as progress:
            for channel in self.instr.channels:
                cubes[channel['name']] = None
                for i,t in enumerate(self.instr.observing_time):
                    i_time = np.where(t == instr_time)[0][0]
                    tmp = Map(self.ff_format.format(channel=channel['name'], i_time=i_time))
                    if cubes[channel['name']] is None:
                        cubes[channel['name']] = np.empty(tmp.data.shape + self.instr.observing_time.shape)
                    cubes[channel['name']][:, :, i] = tmp.data
                    progress.update()
                meta_templates[channel['name']] = tmp.meta.copy()
            
        return cubes,meta_templates
    
    def make_timeseries(self, cubes, channel, left_corner, right_corner):
        tmp = Map(cubes[channel][:,:,0],self.meta_templates[channel])
        blc = tmp.world_to_pixel(SkyCoord(*left_corner, frame=tmp.coordinate_frame))
        trc = tmp.world_to_pixel(SkyCoord(*right_corner, frame=tmp.coordinate_frame))
        ts = cubes[channel][round(blc.y.value):round(trc.y.value),
                            round(blc.x.value):round(trc.x.value),:].mean(axis=(0,1))
        return ts
    
    def correlation_1d(self, channel_a, channel_b, left_corner, right_corner):
        ts_a = self.make_timeseries(self.cubes, channel_a, left_corner, right_corner)
        ts_b = self.make_timeseries(self.cubes, channel_b, left_corner, right_corner)
        ts_a /= ts_a.max()
        ts_b /= ts_b.max()
        cc = np.fft.irfft(np.fft.rfft(ts_a[::-1], n=self.timelags.shape[0])
                          * np.fft.rfft(ts_b, n=self.timelags.shape[0]), n=self.timelags.shape[0])
        return cc
    
    def correlation_2d(self, channel_a, channel_b, **kwargs):
        """
        Compute cross-correlation between two channels for each pixel in an AIA map 
        """
        # Normalize
        max_a = self.cubes[channel_a].max(axis=2)
        v_a = self.cubes[channel_a]/np.where(max_a==0,1,max_a)[:,:,np.newaxis]
        max_b = self.cubes[channel_b].max(axis=2)
        v_b = self.cubes[channel_b]/np.where(max_b==0,1,max_b)[:,:,np.newaxis]
        return self._dask_correlation(v_a, v_b, chunks=kwargs.get('chunks', None))
    
    def _dask_correlation(self, array_a, array_b, chunks=None):
        """
        Create Dask task graph to compute cross-correlation using FFT
        """
        # Create Dask arrays
        if chunks is None:
            chunks = (int(array_a.shape[0]/20),int(array_a.shape[1]/20))+self.instr.observing_time.shape
        darray_a = dask.array.from_array(array_a[:,:,::-1], chunks=chunks)
        darray_b = dask.array.from_array(array_b, chunks=chunks)
        # Build task graph
        fft_a = dask.array.fft.rfft(darray_a, axis=2, n=self.timelags.shape[0])
        fft_b = dask.array.fft.rfft(darray_b, axis=2, n=self.timelags.shape[0])
        cc = dask.array.fft.irfft(fft_a*fft_b, axis=2, n=self.timelags.shape[0])
        return cc
    
    def make_timelag_map(self, channel_a, channel_b, correlation_threshold=1.):
        """
        Compute map of timelag associated with maximum cross-correlation between
        two channels in each pixel of an AIA map.
        """
        cc = self.correlation_2d(channel_a,channel_b).compute()
        max_cc = np.max(cc,axis=2)
        max_timelag = self.timelags[np.argmax(cc,axis=2)]
        max_timelag = np.where(max_cc<correlation_threshold,0,max_timelag)
        # Metadata
        meta = self.meta_templates[channel_a].copy()
        del meta['instrume']
        del meta['t_obs']
        del meta['wavelnth']
        meta_cc = meta.copy()
        meta_cc['bunit'] = ''
        meta_cc['comment'] = f'{channel_a}-{channel_b} cross-correlation'
        meta_timelag = meta.copy()
        meta_timelag['unit'] = 's'
        meta_timelag['comment'] = f'{channel_a}-{channel_b} timelag'

        plot_settings = {'cmap':'plasma'}
        correlation_map = GenericMap(max_cc,meta_cc,plot_settings=plot_settings)
        plot_settings = {'cmap':'RdBu','vmin':self.timelags.value.min(),'vmax':self.timelags.value.max()}
        timelag_map = GenericMap(max_timelag,meta_timelag,plot_settings=plot_settings)
        
        return correlation_map,timelag_map
    
    @staticmethod
    def timelag_map(filename):
        m = Map(filename)
        if 'timelag' in m.meta['comment']:
            m.plot_settings.update({'cmap':'RdBu','vmin':self.timelags.value.min(),'vmax':self.timelags.value.max()})
        elif 'correlation' in m.meta['comment']:
            m.plot_Settings.update({'cmap':'plasma'})
        else:
            warnings.warn('Map does not seem to be either a correlation or timelag map')
                
        return m
        
        