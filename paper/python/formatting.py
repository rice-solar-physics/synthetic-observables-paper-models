"""
Simple utilities for formatting figures nicely
"""
import seaborn


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
    if 'label' in kwargs:
        del kwargs['label']
    ax.step(bins[-2:], [vals[-1], vals[-1]], where='pre', **kwargs)
    ax.vlines(bins[0], ymin, vals[0], **kwargs)
    ax.vlines(bins[-1], ymin, vals[-1], **kwargs)
