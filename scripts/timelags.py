"""
Class for computing timelags from AIA maps. Uses Dask to compute timelags in each pixel in parallel
"""
import os
import glob
import warnings

import h5py
import numpy as np
from scipy.interpolate import interp1d
import dask.array as da
from sunpy.map import Map, GenericMap
from sunpy.util.metadata import MetaDict
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.utils.console import ProgressBar

from synthesizAR.util import get_keys

from aiacube import DistributedAIACollection


class AIATimeLags(DistributedAIACollection):
    """
    Compute AIA timelag maps
    """
    @property
    def timelags(self):
        c = self.channels[0]
        delta_t = np.diff(self[c].time.value).cumsum()
        return np.hstack([-delta_t[::-1], np.array([0]), delta_t]) * self[c].time.unit

    def make_timeseries(self, channel, left_corner, right_corner, **kwargs):
        tmp = self[channel].maps[0].compute()
        x_l, y_l = tmp.world_to_pixel(SkyCoord(*left_corner, frame=tmp.coordinate_frame))
        x_u, y_u = tmp.world_to_pixel(SkyCoord(*right_corner, frame=tmp.coordinate_frame))
        x_l, y_l, x_u, y_u = np.round([x_l.value, y_l.value, x_u.value, y_u.value]).astype(np.int)
        chunks = kwargs.get('chunks', (tmp.data.shape[0]//10, tmp.data.shape[1]//10))
        return (self[channel].rechunk(self[channel].time.shape+chunks)[:, y_l:y_u, x_l:x_u]
                .mean(axis=(1, 2)))

    def correlation_1d(self, channel_a, channel_b, left_corner, right_corner, **kwargs):
        ts_a = self.make_timeseries(channel_a, left_corner, right_corner, **kwargs)
        ts_b = self.make_timeseries(channel_b, left_corner, right_corner, **kwargs)
        ts_a = (ts_a - ts_a.mean()) / ts_a.std()
        ts_b = (ts_b - ts_b.mean()) / ts_b.std()
        cc = da.fft.irfft(da.fft.rfft(ts_a[::-1], n=self.timelags.shape[0])
                          * da.fft.rfft(ts_b, n=self.timelags.shape[0]), n=self.timelags.shape[0])
        return cc

    def correlation_2d(self, channel_a, channel_b, **kwargs):
        """
        Create Dask task graph to compute cross-correlation using FFT for each pixel in an AIA map
        """
        shape = self[channel_a].maps[0].data.shape
        chunks = kwargs.get('chunks', (shape[0]//10, shape[1]//10))
        cube_a = self[channel_a].rechunk(self[channel_a].time.shape+chunks)[::-1, :, :]
        cube_b = self[channel_b].rechunk(self[channel_b].time.shape+chunks)
        # Normalize
        std_a = cube_a.std(axis=0)
        std_a = da.where(std_a == 0, 1, std_a)
        v_a = (cube_a - cube_a.mean(axis=0)[np.newaxis, :, :]) / std_a[np.newaxis, :, :]
        std_b = cube_b.std(axis=0)
        std_b = da.where(std_b == 0, 1, std_b)
        v_b = (cube_b - cube_b.mean(axis=0)[np.newaxis, :, :]) / std_b[np.newaxis, :, :]
        # Fast Fourier Transform of both channels
        fft_a = da.fft.rfft(v_a, axis=0, n=self.timelags.shape[0])
        fft_b = da.fft.rfft(v_b, axis=0, n=self.timelags.shape[0])
        # Inverse of product of FFTS to get cross-correlation (by convolution theorem)
        return da.fft.irfft(fft_a * fft_b, axis=0, n=self.timelags.shape[0])

    def make_correlation_map(self, channel_a, channel_b, **kwargs):
        """
        Build map of max correlation value between two AIA channels
        """
        cc = self.correlation_2d(channel_a, channel_b, **kwargs)
        bounds = kwargs.get('timelag_bounds', None)
        if bounds is not None:
            indices, = np.where(np.logical_and(self.timelags >= bounds[0],
                                               self.timelags <= bounds[1]))
            start = indices[0]
            stop = indices[-1] + 1
        else:
            start = 0
            stop = self.timelags.shape[0] + 1
        max_cc = cc[start:stop, :, :].max(axis=0).compute()
        # Metadata
        meta = self[channel_a].headers[0].copy()
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

        #NOTE: Extra block for computing correlation map is to save time
        """
        cc = self.correlation_2d(channel_a, channel_b, **kwargs)
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
            max_cc = _cc[start:stop, :, :].max(axis=0)
            i_max_cc = _cc[start:stop, :, :].argmax(axis=0)
            del _cc
        else:
            i_max_cc = cc[start:stop, :, :].argmax(axis=0).compute()
        max_timelag = self.timelags[start:stop][i_max_cc]
        # Metadata
        meta = self[channel_a].headers[0].copy()
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
            return timelag_map, correlation_map
        else:
            return timelag_map
