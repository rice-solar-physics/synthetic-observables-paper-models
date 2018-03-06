"""
Class for computing timelags from AIA maps. Uses Dask to compute timelags in each pixel in parallel
"""
import os
import warnings

import h5py
import numpy as np
import dask.array
from sunpy.map import Map,GenericMap
from sunpy.util.metadata import MetaDict
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
    
    def __init__(self, instr, hdf5_filename, fits_root_path=None, **kwargs):
        self.instr = instr
        delta_t = np.diff(instr.observing_time.value).cumsum()
        self.timelags = np.hstack([-delta_t[::-1],np.array([0]),delta_t])*self.instr.observing_time.unit
        self.hdf5_filename = hdf5_filename
        if fits_root_path:
            ff_format = os.path.join(fits_root_path, f'{instr.name}','{channel}','map_t{i_time:06d}.fits')
            self.load_data(ff_format, **kwargs)
    
    def get_metadata(self, channel_name):
        with h5py.File(self.hdf5_filename, 'r') as hf:
            meta = dict(hf[channel_name]['meta'].attrs)
        return MetaDict(meta)
    
    def get_data(self, channel_name):
        hf = h5py.File(self.hdf5_filename, 'r')
        return dask.array.from_array(hf[channel_name]['data'], hf[channel_name]['data'].chunks)
            
        
    def load_data(self, ff_format, **kwargs):
        """
        Create data cubes over desired time range for all channels
        """
        with h5py.File(self.instr.counts_file,'r') as hf:
            instr_time = u.Quantity(hf['time'],hf['time'].attrs['units'])
        with ProgressBar(len(self.instr.channels)*len(self.instr.observing_time),
                         ipython_widget=kwargs.get('notebook',True)) as progress:
            with h5py.File(self.hdf5_filename, 'w') as hf:
                for channel in self.instr.channels:
                    grp = hf.create_group(channel['name'])
                    # Read in tmp file to get shape
                    tmp = Map(ff_format.format(channel=channel['name'], i_time=0))
                    shape = tmp.data.shape + self.instr.observing_time.shape
                    chunks = kwargs.get('chunks', (shape[0]//20, shape[1]//20, shape[2]))
                    dset = grp.create_dataset('data', shape, chunks=chunks)
                    tmp_array = np.zeros(shape)
                    for i,t in enumerate(self.instr.observing_time):
                        i_time = np.where(t == instr_time)[0][0]
                        tmp = Map(ff_format.format(channel=channel['name'], i_time=i_time))
                        tmp_array[:, :, i] = tmp.data
                        progress.update()
                    dset[:, :, :] = tmp_array
                    # Add map metadata as attributes
                    meta = grp.create_group('meta')
                    for k in tmp.meta:
                        try:
                            meta.attrs[k] = tmp.meta[k]
                        except TypeError:
                            continue
                
    def make_timeseries(self, channel, left_corner, right_corner):
        darray = self.get_data(channel)
        tmp = Map(np.array(darray[:,:,0]), self.get_metadata(channel))
        blc = tmp.world_to_pixel(SkyCoord(*left_corner, frame=tmp.coordinate_frame))
        trc = tmp.world_to_pixel(SkyCoord(*right_corner, frame=tmp.coordinate_frame))
        ts = darray[round(blc.y.value):round(trc.y.value),
                    round(blc.x.value):round(trc.x.value), :].mean(axis=(0,1)).compute()
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
        Create Dask task graph to compute cross-correlation using FFT for each pixel in an AIA map 
        """
        darray_a = self.get_data(channel_a)
        darray_b = self.get_data(channel_b)
        # Normalize
        max_a = darray_a.max(axis=2)
        v_a = darray_a / dask.array.where(max_a==0, 1, max_a)[:,:,np.newaxis]
        max_b = darray_b.max(axis=2)
        v_b = darray_b / dask.array.where(max_b==0, 1, max_b)[:,:,np.newaxis]
        # Fast Fourier Transform of both channels
        fft_a = dask.array.fft.rfft(v_a, axis=2, n=self.timelags.shape[0])
        fft_b = dask.array.fft.rfft(v_b, axis=2, n=self.timelags.shape[0])
        # Inverse of product of FFTS to get cross-correlation
        return dask.array.fft.irfft(fft_a*fft_b, axis=2, n=self.timelags.shape[0])
    
    def make_timelag_map(self, channel_a, channel_b, correlation_threshold=1., **kwargs):
        """
        Compute map of timelag associated with maximum cross-correlation between
        two channels in each pixel of an AIA map.
        """
        cc = self.correlation_2d(channel_a, channel_b, **kwargs)
        max_cc = cc.max(axis=2).compute()
        i_max_cc = cc.argmax(axis=2).compute()
        max_timelag = self.timelags[i_max_cc]
        max_timelag = np.where(max_cc < correlation_threshold, np.nan, max_timelag)
        # Metadata
        meta = self.get_metadata(channel_a).copy()
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
        correlation_map = GenericMap(max_cc, meta_cc, plot_settings=plot_settings)
        plot_settings = {'cmap':'RdBu', 'vmin':self.timelags.value.min(), 'vmax':self.timelags.value.max()}
        timelag_map = GenericMap(max_timelag, meta_timelag, plot_settings=plot_settings)
        
        return correlation_map, timelag_map
    
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
        
        