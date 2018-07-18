"""
Object for dealing efficiently with large out-of-core AIA datacubes.
Much of this was inspired by work done by Stuart Mumford, in particular how to handle
FITS files in Dask.
"""
import dask.bytes
import dask.array as da
import dask
import distributed
from dask import delayed,compute
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
    Load sequences of AIA images over multiple channels and operate on them
    in a distributed and parallelized way.
    """
    def __init__(self, read_template, channels=None):
        client = distributed.get_client()
        self.channels = ['94','131','171','193','211','335'] if channels is None else channels
        # Open files with Dask
        templates = {c: read_template.format(c) for c in self.channels}
        openfiles = {c: dask.bytes.open_files(temp) for c,temp in templates.items()}
        # Get headers from all FITS files
        self.headers = self._get_headers(openfiles)
        # Get dtype and shape
        self.dtype,self.shape = self._get_dtype_and_shape()
        # Construct maps
        self.maps = self._get_maps(openfiles)
        
    def _get_maps(self,openfiles,):
        arrays = {c: [da.from_array(DelayedFITS(f,shape=self.shape,dtype=self.dtype,hdu=0),chunks=self.shape) 
                      for f in files] 
                  for c,files in openfiles.items()}
        maps = {}
        for c in arrays.keys():
            maps[c] = [delayed(sunpy.map.sources.AIAMap)(a, h) for a,h in zip(arrays[c], self.headers[c])]
        return maps
        
    def _get_headers(self, openfiles):
        client = distributed.get_client()
        futures = {c: client.map(get_header, fn) for c,fn in openfiles.items()}
        return {c: client.gather(fut) for c,fut in futures.items()}
    
    def _get_dtype_and_shape(self,):
        dtypes = {c: [validate_dtype_shape(h) for h in head] for c, head in self.headers.items()}
        if not all([d == list(dtypes.values())[0] for d in dtypes.values()]):
            raise ValueError('All maps must have same shape and dtype')
        return list(dtypes.values())[0][0]
    
    @property
    def time(self,):
        """
        Note that this should really check for date-obs and then do something different if that exists
        """
        return u.Quantity([h['t_obs'] for h in self.headers[self.channels[0]]], 
                          self.headers[self.channels[0]][0]['tunit'])
    
    @property
    def unstacked_data(self,):
        return {c: [da.from_delayed(m.data, dtype=self.dtype, shape=self.shape) for m in maps] for c,maps in self.maps.items()}
    
    @property
    def stacked_data(self,):
        return {c: da.stack(array) for c,array in self.unstacked_data.items()}
        
    def rechunk(self,shape):
        return {c: d.rechunk(shape) for c,d in self.stacked_data.items()}
    