import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import correlate, correlation_lags
from astropy.io import fits
from astropy import units as u

def crosscor(spec_a, spec_b):
    """
    Normalised cross-correlation of two 1-D spectra.

    Returns
    -------
    lags   : pixel-lag array
    xcorr  : normalised cross-correlation array
    peak   : lag at the peak (integer pixels)
    """
    a = mask_nans(spec_a - np.nanmean(spec_a))
    b = mask_nans(spec_b - np.nanmean(spec_b))

    cc = correlate(a, b, mode="full")
    norm = np.sqrt(np.dot(a, a) * np.dot(b, b))
    if norm > 0:
        cc /= norm

    lags = correlation_lags(len(a), len(b), mode="full")
    peak = lags[np.argmax(cc)]
    return lags, cc, peak