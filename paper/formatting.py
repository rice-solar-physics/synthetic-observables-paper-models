"""
Few simple utilities for formatting figures nicely
"""
import seaborn
import matplotlib.colors


# Nice looking RWB diverging colormap
def rwb_cmap():
    p = seaborn.color_palette('deep')
    return matplotlib.colors.LinearSegmentedColormap.from_list(
        'rwb_nice', [p[0], (1, 1, 1), p[3],],N=1000,
    )


# Color palette for heating functions
def heating_palette():
    return seaborn.color_palette(palette='colorblind', n_colors=3,)


def heating_cmap():
    return matplotlib.colors.ListedColormap(heating_palette(),N=3),


# General qualitative color palette
def qualitative_palette(n):
    return seaborn.color_palette('Dark2', n)


# Function for consistent figure sizing in TeX
# Default columnwidth is hardcoded for this paper
def get_figsize(wf=0.5, hf=(5.**0.5-1.0)/2.0, columnwidth=513.11743):
    """
    Calculate figure sizes to get consistent font sizes between latex and 
    matplotlib. See `this SO post <https://stackoverflow.com/questions/29187618/matplotlib-and-latex-beamer-correct-size>`_.
    For this to work, your figures in LaTeX should NOT be scaled.
    
    Parameters
    ----------
    wf : `float`
        width fraction in columnwidth units
    hf : `float`
        height fraction in columnwidth units. Set by default to golden ratio.
    columnwidth : `float`
        width of the column in latex. Get this from LaTeX using \showthe\columnwidth
    Returns
    -------
    (fig_width, fig_height) : `tuple
        that should be given to matplotlib
    """
    fig_width_pt = columnwidth*wf 
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*hf      # height in inches
    return (fig_width, fig_height)


# Custom histogram plotting that respects rcparams
def hist_step(ax, vals, bins, **kwargs):
    """
    Alternative to matplotlib hist that respects matplotlibrc styling
    
    Parameters
    ----------
    ax : Matplotlib axis
    vals : value in each bin
    bins : Bin edges, including the rightmost edge
    kwargs : Plotting keyword arguments
    """
    ymin = ax.get_ylim()[0]
    ax.step(bins[:-1], vals, where='post', **kwargs)
    ax.step(bins[-2:],[vals[-1],vals[-1]], where='pre', **kwargs)
    ax.vlines(bins[0], ymin, vals[0], **kwargs)
    ax.vlines(bins[-1], ymin, vals[-1], **kwargs)
    