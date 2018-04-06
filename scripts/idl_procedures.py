"""
Miscellaneous IDL procedures in Python
"""
import numpy as np

def c_correlate(s_1, s_2, lags):
    """
    Numpy implementation of c_correlate.pro IDL routine
    """
    # ensure signals are of equal length
    assert s_1.shape == s_2.shape
    n_s = s_1.shape[0]
    # center both signals
    s_1_center = s_1 - s_1.mean()
    s_2_center = s_2 - s_2.mean()
    # allocate space for correlation
    correlation = np.zeros(lags.shape)
    # iterate over lags
    for i,l in enumerate(lags):
        if l >= 0:
            tmp = s_1_center[:(n_s - l)] * s_2_center[l:]
        else:
            tmp = s_1_center[-l:] * s_2_center[:(n_s + l)]
        correlation[i] = tmp.sum()
    # Divide by standard deviation of both
    correlation /= np.sqrt((s_1_center**2).sum() * (s_2_center**2).sum())
    
    return correlation