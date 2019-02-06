"""
Hack at computing DEM in a nice way
Move this somewhere else!
"""
import os

import numpy as np
import matplotlib.colors
import astropy.units as u
from sunpy.map import Map
import hissw

from synthesizAR.analysis.dem import EMCube

__all__ = ['GenericModel', 'IDLModel', 'HannahKontarModel', 'make_slope_map_tpeak']


class GenericModel(object):
    
    @u.quantity_input
    def __init__(self, maps, temperature: u.K, responses, **kwargs):
        self.temperature = temperature
        wvl = u.Quantity([m.meta['wavelnth']*u.Unit(m.meta['waveunit']) for m in maps])
        self.maps = [maps[i] for i in np.argsort(wvl)]
        self.response = [responses[i] for i in np.argsort(wvl)]
        self.wavelength = np.sort(wvl)
        
    def fit(self,):
        pass


class IDLModel(GenericModel):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if 'dem_path' not in kwargs:
            raise ValueError('Must specify path to DEM IDL code')
        self.ssw = hissw.ScriptMaker(ssw_packages=self.ssw_packages, ssw_paths=self.ssw_paths,
                                     extra_paths=[kwargs.get('dem_path')])
        self.mapcube = np.stack([m.data for m in self.maps], axis=2)
        self.response_matrix = np.stack(self.response)
        
    def _fit(self, **kwargs):
        input_args = self.input_args
        input_args.update(kwargs)
        return self.ssw.run(self._template, args=input_args, save_vars=self.save_vars,
                            verbose=kwargs.get('verbose',True),)
    
    def to_emcube(self, em):
        header = self.maps[0].meta.copy()
        del header['telescop']
        del header['detector']
        del header['waveunit']
        del header['instrume']
        del header['wavelnth']
        return EMCube(em, header, self.temperature,
                      plot_settings={'cmap': 'magma',
                                     'norm': matplotlib.colors.SymLogNorm(1, vmin=1e25, vmax=1e30)})


class HannahKontarModel(IDLModel):

    @property
    def ssw_packages(self,):
        return ['sdo/aia']

    @property
    def ssw_paths(self,):
        return ['aia', 'xrt']

    @property
    def input_args(self,):
        return {
            'log_temperature': np.log10(((self.temperature[1:]
                                          + self.temperature[:-1])/2.).to(u.K).value).tolist(),
            'temperature_bin_edges': self.temperature.to(u.K).value.tolist(),
            'saturation_level': 2e4,
            'edge_trim': 10,
            'aia_error_table_filename': os.path.join(self.ssw.ssw_home,
                                                     'sdo/aia/response/aia_V2_error_table.txt'),
            'n_sample': 1,
            'maps': self.mapcube.T.tolist(),
            'response_matrix': self.response_matrix.tolist(),
            'normalized': True,
            'max_iterations': 10,
            'alpha': 1.0,
            'increase_alpha': 1.5,
            'use_em_loci': False,
            'em_loci_indices': None,
            'initial_guess': None,
        }

    @property
    def save_vars(self,):
        return ['dem', 'dem_errors', 'logt_errors', 'chi_squared', 'dn_regularized']

    def fit(self, **kwargs):
        self.results = self._fit(**kwargs)
        em = np.transpose(self.results['dem'], axes=(2, 1, 0))/(u.cm**5)/u.K
        em = em*np.diff(self.temperature)
        return self.to_emcube(em)

    @property
    def _template(self,):
        return """
        data = {{maps}}
        ; Set some basic parameters
        nx=n_elements(data[*,0,0])
        ny=n_elements(data[0,*,0])
        nchannels=n_elements(data[0,0,*])
        
        ; May need exposure times (these are just taken from some FITS files)
        exptimes = [2.901051,2.901311,2.000197,1.999629,2.901278,2.900791]

        ; Filter out bad pixels
        sat_lvl={{ saturation_level }}
        id=where(data ge sat_lvl,nid)
        if (nid gt 1) then data[id]=0.0
        id=where(data le 0,nid)
        if (nid gt 1) then data[id]=0.0
        edg0={{edge_trim}}
        data[0:edg0-1,*,*]=0.0
        data[*,0:edg0-1,*]=0.0
        data[nx-edg0:nx-1,*,*]=0.0
        data[*,ny-edg0:ny-1,*]=0.0

        ; Calculate errors
        {% if percent_error is defined %}
        data_errors = {{ percent_error }}*data
        {% else %}
        data_errors=fltarr(nx,ny,nchannels)
        channels = [94,131,171,193,211,335]
        common aia_bp_error_common,common_errtable
        common_errtable=aia_bp_read_error_table('{{ aia_error_table_filename }}')
        for i=0,nchannels-1 do begin
            norm_factor=1
            {% if normalized %}norm_factor = exptimes[i]{% endif %}
            data_errors[*,*,i]=aia_bp_estimate_error(reform(data[*,*,i])*norm_factor,replicate(channels[i],nx,ny),n_sample={{ n_sample }})
            data_errors[*,*,i]=data_errors[*,*,i]/exptimes[i]
        endfor
        {% endif %}

        ; Get temperature bins
        response_logt = {{log_temperature}}
        temperature = {{temperature_bin_edges}}

        ; Calculate response functions
        response_matrix = {{ response_matrix }}

        ; DEM Calculation
        {% if not normalized %}
        for i=0,nchannels-1 do begin
            data[*,*,i] = data[*,*,i]/exptimes[i]
        endfor
        {% endif %}
        dn2dem_pos_nb,data,data_errors,response_matrix,response_logt,temperature,$
        dem,dem_errors,logt_errors,chi_squared,dn_regularized,{% if use_em_loci %}/gloci,{% endif %}$
        reg_tweak={{alpha}}, rgt_fact={{increase_alpha}}, max_iter={{max_iterations}}
        """


@u.quantity_input
def make_slope_map_tpeak(emcube, Tmin=1e6*u.K, rsquared_tolerance=0.9, safety=10**(0.4)):
    """
    Fit emission measure slope using the temperature at which the emission
    measure peaks as the upper bound.
    
    Parameters
    ----------
    emcube {[type]} -- [description]
    Tmin {[type]} -- [description] (default: {1e6*u.K})
    rsquared_tolerance {float} -- [description] (default: {0.9})
    safety {[type]} -- [description] (default: {10**(0.4)})
    
    Returns:
        [type] -- [description]
    """
    # Choose Tmax just below Tpeak to avoid fitting across the peak
    index = emcube.as_array().argmax(axis=2) - 1
    Tpeak = emcube.temperature_bin_centers[np.where(index < 0, 0, index)]  # No negative indices
    # Fit slopes for multiple Tmax values
    slopes, rsquared, Tmax = [], [], []
    for tp in np.unique(Tpeak):
        if tp < Tmin*safety:
            continue
        slope_map, _, r2 = emcube.make_slope_map(
            temperature_bounds=u.Quantity((Tmin, tp)),
            em_threshold=1*u.cm**(-5),
            rsquared_tolerance=0.0,
            full=True,
        )
        slopes.append(slope_map.data)
        rsquared.append(r2.reshape(slope_map.data.shape))
        Tmax.append(tp.value)
    Tmax = u.Quantity(Tmax, tp.unit)
    # Stack
    rsquared = np.stack(rsquared, axis=2)
    rsquared[np.isnan(rsquared)] = -np.inf  # Give NaNs a very small r^2
    slopes = np.stack(slopes, axis=2)
    # Find index of peak temperature
    i_peak = np.argmin(np.stack([np.fabs(Tpeak - tm) for tm in Tmax], axis=2), axis=2)
    i_row, i_col = np.indices(i_peak.shape)
    # Select slope and rsquared for correct peak temperature
    slopes = slopes[(i_row.flatten(), i_col.flatten(), i_peak.flatten())].reshape(i_peak.shape)
    rsquared = rsquared[(i_row.flatten(), i_col.flatten(), i_peak.flatten())].reshape(i_peak.shape)

    return Map(
        slopes,
        {**slope_map.meta, 'temp_b': 'temperature at peak of EM'},
        mask=rsquared < rsquared_tolerance,
        plot_settings=slope_map.plot_settings)
