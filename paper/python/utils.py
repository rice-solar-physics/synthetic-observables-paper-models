"""
Simple utilities for formatting figures nicely
"""
import warnings

import numpy as np
import seaborn
import astropy.units as u
from sunpy.map import Map


# Color palette for heating functions
def heating_palette(n_colors=3):
    return seaborn.color_palette(palette='colorblind', n_colors=n_colors,)


# General qualitative color palette
def qualitative_palette(n):
    return seaborn.color_palette('Dark2', n)


# Custom histogram plotting that respects rcparams
def plot_hist(ax, vals, bins, **kwargs):
    """
    Plot a histogram from bin values

    Parameters
    ----------
    ax : Matplotlib axis
    vals : value in each bin
    bins : Bin edges, including the rightmost edge
    kwargs : Plotting keyword arguments
    """
    ymin = ax.get_ylim()[0]
    ax.step(bins[:-1], vals, where='post', **kwargs)
    # Remove keywords we don't want duplicated
    for kw in ['label', 'marker']:
        if kw in kwargs:
            del kwargs[kw]
    ax.step(bins[-2:], [vals[-1], vals[-1]], where='pre', **kwargs)
    ax.vlines(bins[0], ymin, vals[0], **kwargs)
    ax.vlines(bins[-1], ymin, vals[-1], **kwargs)


def make_slope_map(emcube, temperature_lower_bound=None, em_threshold=None):
    """
    Fit emission measure distribution in every pixel
    
    The fit is computed between `temperature_lower_bound`
    and the temeperature at which the EM is maximum.
    
    Parameters
    ----------
    emcube: `EMCube`
        Emission measure map as a function space and temperature
    em_threshold: `~astropy.units.Quantity`, optional
        If the total EM in a pixel is below this, no slope is calculated
        
    Returns
    -------
    slope_map: `~sunpy.map.GenericMap`
    rsquared_map: `~sunpy.map.GenericMap`
    """
    if em_threshold is None:
        em_threshold = u.Quantity(1e25, u.cm**(-5))
    i_valid = np.where(
        u.Quantity(emcube.total_emission.data, emcube[0].meta['bunit']) > em_threshold)
    em_valid = np.log10(emcube.as_array()[i_valid])
    em_valid[np.logical_or(np.isinf(em_valid), np.isnan(em_valid))] = 0.0
    i_peak = em_valid.argmax(axis=1)
    log_temperature_bin_centers = np.log10(emcube.temperature_bin_centers.value)
    if temperature_lower_bound is None:
        i_lower = 0
    else:
        i_lower = np.fabs(emcube.temperature_bin_centers - temperature_lower_bound).argmin()
    slopes, rsquared = [], []
    for emv, ip in zip(em_valid, i_peak):
        t_fit = log_temperature_bin_centers[i_lower:ip]
        if t_fit.size < 3:
            warnings.warn('Fit should be over 3 or more bins in temperature.')
        if t_fit.size == 0:
            slopes.append(np.nan)
            rsquared.append(0.)
            continue
        em_fit = emv[i_lower:ip]
        w = np.where(em_fit > 0, 1, 0)
        coeff, rss, _, _, _ = np.polyfit(t_fit, em_fit, 1, full=True, w=w)
        rss = 1 if rss.size == 0 else rss[0]
        _, rss_flat, _, _, _ = np.polyfit(t_fit, em_fit, 0, full=True, w=w)
        rss_flat = 1 if rss_flat.size == 0 else rss_flat[0]
        slopes.append(coeff[0])
        rsquared.append(1-rss/rss_flat)
    slopes_data = np.zeros(emcube.total_emission.data.shape)
    slopes_data[i_valid] = slopes
    rsquared_data = np.zeros(emcube.total_emission.data.shape)
    rsquared_data[i_valid] = rsquared
    
    return Map(slopes_data, emcube[0].meta,), Map(rsquared_data, emcube[0].meta)