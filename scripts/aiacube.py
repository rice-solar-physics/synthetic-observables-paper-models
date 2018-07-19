"""
Object for dealing efficiently with large out-of-core AIA datacubes.
Nearly all of this is inspired by work done by Stuart Mumford, in particular how to handle
FITS files in Dask.
"""
import warnings
import dask.bytes
import dask.array as da
import dask
import distributed
from dask import delayed, compute
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from astropy.io.fits.hdu.base import BITPIX2DTYPE
import numpy as np
import sunpy.map
from sunpy.util.metadata import MetaDict

__all__ = ['DistributedAIACube',]


def validate_dtype_shape(head):
    naxes = head['NAXIS']
    dtype = BITPIX2DTYPE[head['BITPIX']]
    shape = [head[f'NAXIS{n}'] for n in range(naxes, 0, -1)]
    return dtype, shape


def get_header(fn, hdu=0):
    with fn as fi:
        return MetaDict(sunpy.io.fits.get_header(fi)[hdu])
    

class DelayedFITS:
    def __init__(self, file, shape, dtype, hdu=0):
        self.shape = shape
        self.dtype = dtype
        self.file = file
        self.hdu = hdu
    
    def __getitem__(self, item):
        with self.file as fi:
            with fits.open(fi) as hdul:
                return hdul[self.hdu].section[item]
            

class DistributedAIACube(object):
    """
    Load sequence of AIA images for a single channel and operate on them in a distributed and 
    parallelized way.

    #TODO: Refactor this to use ndcube instead
    """
    def __init__(self, maps, headers):
        if not all([m.data.shape == maps[0].data.shape for m in maps]):
            raise ValueError('All maps must have same dimensions')
        if not all([m.data.dtype == maps[0].data.dtype for m in maps]):
            raise ValueError('All maps must have same dtype')
        self.maps = maps
        self.headers = headers

    @classmethod
    def from_files(cls, read_template):
        openfiles = dask.bytes.open_files(read_template)
        headers = cls._get_headers(openfiles)
        dtype, shape = cls._get_dtype_and_shape(headers)
        maps = cls._get_maps(openfiles, headers, dtype, shape)
        return cls(maps, headers)

    @staticmethod
    def _get_maps(openfiles, headers, dtype, shape):
        arrays = [da.from_array(DelayedFITS(f, shape=shape, dtype=dtype, hdu=0),
                                chunks=shape) for f in openfiles]
        return [sunpy.map.Map(a, h) for a, h in zip(arrays, headers)]

    @staticmethod
    def _get_headers(openfiles):
        client = distributed.get_client()
        futures = client.map(get_header, openfiles)
        return client.gather(futures)

    @staticmethod
    def _get_dtype_and_shape(headers):
        dtypes = [validate_dtype_shape(h) for h in headers]
        if not all([d == dtypes[0] for d in dtypes]):
            raise ValueError('All maps must have same shape and dtype')
        return dtypes[0]

    @property
    def time(self,):
        """
        Note that this should really check for date-obs and then do something different if that
        exists
        """
        return u.Quantity([h['t_obs'] for h in self.headers],
                          self.headers[self.channels[0]][0]['tunit'])

    @property
    def shape(self,):
        return self.maps[0].data.shape

    @property
    def dtype(self,):
        return self.maps[0].data.dtype

    @property
    def unstacked_data(self,):
        return [da.from_delayed(m.data, dtype=m.data.dtype, shape=m.data.shape) for m in self.maps]

    @property
    def stacked_data(self,):
        return da.stack(self.unstacked_data)

    def rechunk(self, shape):
        return self.stacked_data.rechunk(shape)

    def prep(self,):
        """
        Lazily apply a prep operation to all maps and return a new distributed cube object
        """
        pass

    def derotate(self, reference_date):
        """
        Lazily apply a derotation to all maps and return a new distributed cube object
        """
        pass


class DistributedAIACollection(object):
    """
    A collection of distributed AIA images over multiple wavelengths

    #TODO: refactor this to use ndcube sequence
    #TODO: figure out useful methods to add here
    #TODO: Add a check for data being aligned?
    #NOTE: It is assumed that the data in this container are all aligned, i.e. prepped and derotated
    """

    def __init__(self, *args, **kwargs):
        # Check all spatial and time shapes the same
        if not all([a.shape == args[0].shape for a in args]):
            raise ValueError('All spatial dimensions must be the same')
        if not all([a.time.shape == args[0].time.shape for a in args]):
            # Not an error because may have missing timesteps in observations
            # Will interpolate later to account for this
            raise warnings.warn('Time dimensions are not all equal length')
        self._cubes = {f"{a.headers[0]['wavelnth']}": a for a in args}
        self.channels = self._cubes.keys()

    def __getitem__(self, channel):
        channel = f'{channel}'
        return self._cubes[channel]
