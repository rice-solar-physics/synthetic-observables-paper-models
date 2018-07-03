"""
Class for computing timelags from AIA maps. Uses Dask to compute timelags in each pixel in parallel
"""
import os
import glob
import warnings

import h5py
import numpy as np
from scipy.interpolate import interp1d
import dask.array
from sunpy.map import Map, GenericMap
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
        self.timelags = (np.hstack([-delta_t[::-1], np.array([0]), delta_t])
                         * self.instr.observing_time.unit)
        self.hdf5_filename = hdf5_filename
        if fits_root_path:
            ff_format = os.path.join(fits_root_path,
                                     f'{instr.name}', '{channel}', 'map_t{i_time:06d}.fits')
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
        with h5py.File(self.instr.counts_file, 'r') as hf:
            instr_time = u.Quantity(hf['time'], hf['time'].attrs['units'])
        with ProgressBar(len(self.instr.channels) * len(self.instr.observing_time),
                         ipython_widget=kwargs.get('notebook', True)) as progress:
            with h5py.File(self.hdf5_filename, 'w') as hf:
                for channel in self.instr.channels:
                    grp = hf.create_group(channel['name'])
                    # Read in tmp file to get shape
                    i_time = np.where(self.instr.observing_time[0] == instr_time)[0][0]
                    tmp = Map(ff_format.format(channel=channel['name'], i_time=i_time))
                    shape = tmp.data.shape + self.instr.observing_time.shape
                    chunks = kwargs.get('chunks', (shape[0]//20, shape[1]//20, shape[2]))
                    dset = grp.create_dataset('data', shape, chunks=chunks)
                    tmp_array = np.zeros(shape)
                    for i, t in enumerate(self.instr.observing_time):
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
        tmp = Map(np.array(darray[:, :, 0]), self.get_metadata(channel))
        blc = tmp.world_to_pixel(SkyCoord(*left_corner, frame=tmp.coordinate_frame))
        trc = tmp.world_to_pixel(SkyCoord(*right_corner, frame=tmp.coordinate_frame))
        ts = darray[round(blc.y.value):round(trc.y.value),
                    round(blc.x.value):round(trc.x.value), :].mean(axis=(0, 1)).compute()
        return ts

    def correlation_1d(self, channel_a, channel_b, left_corner, right_corner):
        ts_a = self.make_timeseries(channel_a, left_corner, right_corner)
        ts_b = self.make_timeseries(channel_b, left_corner, right_corner)
        ts_a = (ts_a - ts_a.mean()) / ts_a.std()
        ts_b = (ts_b - ts_b.mean()) / ts_b.std()
        cc = np.fft.irfft(np.fft.rfft(ts_a[::-1], n=self.timelags.shape[0])
                          * np.fft.rfft(ts_b, n=self.timelags.shape[0]), n=self.timelags.shape[0])
        return cc

    def correlation_2d(self, channel_a, channel_b):
        """
        Create Dask task graph to compute cross-correlation using FFT for each pixel in an AIA map
        """
        darray_a = self.get_data(channel_a)[:, :, ::-1]
        darray_b = self.get_data(channel_b)
        # Normalize
        std_a = darray_a.std(axis=2)
        std_a = dask.array.where(std_a == 0, 1, std_a)
        v_a = (darray_a - darray_a.mean(axis=2)[:, :, np.newaxis]) / std_a[:, :, np.newaxis]
        std_b = darray_b.std(axis=2)
        std_b = dask.array.where(std_b == 0, 1, std_b)
        v_b = (darray_b - darray_b.mean(axis=2)[:, :, np.newaxis]) / std_b[:, :, np.newaxis]
        # Fast Fourier Transform of both channels
        fft_a = dask.array.fft.rfft(v_a, axis=2, n=self.timelags.shape[0])
        fft_b = dask.array.fft.rfft(v_b, axis=2, n=self.timelags.shape[0])
        # Inverse of product of FFTS to get cross-correlation (by convolution theorem)
        return dask.array.fft.irfft(fft_a * fft_b, axis=2, n=self.timelags.shape[0])

    def make_correlation_map(self, channel_a, channel_b, **kwargs):
        """
        Build map of max correlation value between two AIA channels
        """
        cc = self.correlation_2d(channel_a, channel_b)
        bounds = kwargs.get('timelag_bounds', None)
        if bounds is not None:
            indices, = np.where(np.logical_and(self.timelags >= bounds[0],
                                               self.timelags <= bounds[1]))
            start = indices[0]
            stop = indices[-1] + 1
        else:
            start = 0
            stop = self.timelags.shape[0] + 1
        max_cc = cc[:,:,start:stop].max(axis=2).compute()
        # Metadata
        meta = self.get_metadata(channel_a).copy()
        del meta['instrume']
        del meta['t_obs']
        del meta['wavelnth']
        meta['bunit'] = ''
        meta['comment'] = f'{channel_a}-{channel_b} cross-correlation'

        plot_settings = {'cmap': 'plasma'}
        plot_settings.update(kwargs.get('plot_settings', {}))
        correlation_map = GenericMap(max_cc, meta, plot_settings=plot_settings)

        return correlation_map

    def make_timelag_map(self, channel_a, channel_b, **kwargs):
        """
        Compute map of timelag associated with maximum cross-correlation between
        two channels in each pixel of an AIA map.
        """
        cc = self.correlation_2d(channel_a, channel_b)
        bounds = kwargs.get('timelag_bounds', None)
        if bounds is not None:
            indices, = np.where(np.logical_and(self.timelags >= bounds[0],
                                               self.timelags <= bounds[1]))
            start = indices[0]
            stop = indices[-1] + 1
        else:
            start = 0
            stop = self.timelags.shape[0] + 1
        if kwargs.get('return_correlation_map'):
            _cc = cc.compute()
            max_cc = _cc[:,:,start:stop].max(axis=2)
            i_max_cc = _cc[:,:,start:stop].argmax(axis=2)
            del _cc
        else:
            i_max_cc = cc[:,:,start:stop].argmax(axis=2).compute()
        max_timelag = self.timelags[start:stop][i_max_cc]
        # Metadata
        meta = self.get_metadata(channel_a).copy()
        del meta['instrume']
        del meta['t_obs']
        del meta['wavelnth']
        meta['bunit'] = 's'
        meta['comment'] = f'{channel_a}-{channel_b} timelag'
        plot_settings = {'cmap': 'RdBu_r', 'vmin': self.timelags[start:stop].value.min(),
                         'vmax': self.timelags[start:stop].value.max()}
        plot_settings.update(kwargs.get('plot_settings', {}))
        timelag_map = GenericMap(max_timelag, meta.copy(), plot_settings=plot_settings.copy())
        if kwargs.get('return_correlation_map'):
            meta['bunit'] = ''
            meta['comment'] = f'{channel_a}-{channel_b} cross-correlation'
            plot_settings['cmap'] = 'plasma'
            del plot_settings['vmin']
            del plot_settings['vmax']
            correlation_map = GenericMap(max_cc, meta, plot_settings=plot_settings)
            return timelag_map,correlation_map
        else:
            return timelag_map
    
    
class AIATimeLagsObserved(AIATimeLags):
    
    def __init__(self, hdf5_filename, fits_root_path=None, **kwargs):
        self.hdf5_filename = hdf5_filename
        if fits_root_path is not None:
            self.load_data(fits_root_path, **kwargs)
    
    def load_data(self, fits_root_path, **kwargs):
        """
        Save all AIA maps to HDF5 file. Note that these maps must already be prepped,
        derotated, and cropped.
        """
        channels = [94, 131, 171, 193, 211, 335]
        cadence = 12.0
        # Find time interval and reference time for each channel
        ref_time = {}
        time_min, time_max = None, None
        for chan in channels:
            tmp_list = sorted(glob.glob(os.path.join(fits_root_path, f'*_{chan}_*.fits')))
            tmp = Map(tmp_list[0], tmp_list[-1])
            if time_min is None or tmp[0].date < time_min:
                time_min = tmp[0].date
            if time_max is None or tmp[1].date > time_max:
                time_max = tmp[1].date
            ref_time[f'{chan}'] = tmp[0].date
            shape = tmp[0].data.shape
        # Save common time
        common_time = np.arange(0., (time_max - time_min).seconds, cadence) * u.s
        with h5py.File(self.hdf5_filename, 'w') as hf:
            dset = hf.create_dataset('observing_time', data=common_time.value)
            dset.attrs['units'] = common_time.unit.to_string()
        chunks = kwargs.get('chunks', (shape[0]//20, shape[1]//20, common_time.shape[0]))
        # Load all maps from all channels
        with ProgressBar(len(glob.glob(os.path.join(fits_root_path, '*'))), ipython_widget=kwargs.get('notebook', True)) as progress:
            for chan in channels:
                files = sorted(glob.glob(os.path.join(fits_root_path, f'*_{chan}_*.fits')))
                # Load all data into memory
                tmp_time = np.zeros((len(files),))
                tmp_array = np.empty(shape + (len(files),))
                for i,f in enumerate(files):
                    # load file
                    tmp_map = Map(f)
                    tmp_array[:, :, i] = tmp_map.data
                    # get time
                    tmp_time[i] = (tmp_map.date - ref_time[f'{chan}']).seconds
                    progress.update()
                # interpolate to common time
                tmp_array = interp1d(tmp_time, tmp_array, axis=2, fill_value='extrapolate', 
                                     kind='linear', bounds_error=False, assume_sorted=True)(common_time.value)
                # Load all data into an HDF5 file
                with h5py.File(self.hdf5_filename, 'a') as hf:
                    grp = hf.create_group(f'{chan}')
                    dset = grp.create_dataset('data', shape+common_time.shape, chunks=chunks)
                    dset[:,:,:] = tmp_array
                    meta = grp.create_group('meta')
                    for k in tmp_map.meta:
                        try:
                            meta.attrs[k] = tmp_map.meta[k]
                        except TypeError:
                            continue
                    
    @property
    def timelags(self):
        delta_t = np.diff(self.observing_time.value).cumsum()
        return (np.hstack([-delta_t[::-1], np.array([0]), delta_t]) 
                * self.observing_time.unit)
    
    @property
    def observing_time(self):
        with h5py.File(self.hdf5_filename,'r') as hf:
            obs_time = u.Quantity(hf['observing_time'], hf['observing_time'].attrs['units'])
        return obs_time
