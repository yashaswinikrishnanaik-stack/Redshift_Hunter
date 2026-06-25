import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import correlate, correlation_lags
from astropy.io import fits
from astropy import units as u

#functions to make both spectra into a uniform log-wavelength grid sinceredshift becomes an additive shift in the log-wavelength space
def _make_log_grid(wave_min: float, wave_max: float, dloglam: float) -> np.ndarray:
    n = int(np.ceil((np.log(wave_max) - np.log(wave_min)) / dloglam))
    n = max(n, 2)
    return np.log(wave_min) + np.arange(n) * dloglam

def _resample_to_log(wave: np.ndarray, flux: np.ndarray, loglam_grid: np.ndarray) -> np.ndarray:
    loglam = np.log(wave)
    order = np.argsort(loglam)
    return np.interp(loglam_grid, loglam[order], flux[order], left=np.nan, right=np.nan)


def _normalize(flux: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean-subtract / std-normalize, treating non-finite pixels as masked."""
    mask = np.isfinite(flux)
    out = np.zeros_like(flux, dtype=float)
    if mask.sum() < 2:
        return out, mask
    mean = flux[mask].mean()
    std = flux[mask].std()
    if std == 0 or not np.isfinite(std):
        return out, mask
    out[mask] = (flux[mask] - mean) / std
    return out, mask

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