import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import correlate, correlation_lags
from astropy.io import fits
from astropy import units as u
from spectra.py import observed_spectra, Eigenspectra

def find_redshift(obs_wave, obs_flux, tmpl_wave, tmpl_flux,
                   z_min=-0.01, z_max=1.0, dloglam=1e-4, min_overlap=50):

    # function to resample spectra into a log-wavelength grid    
        def log_grid(wave,flux):
            loglam0 = np.log(wave.min())
            n = int(np,ceil((np.log(wave.max)) - loglam0) /dloglam)
            grid = loglam0 + np.arange(max(n, 2)) * dloglam
            order = np.argsort(wave)
            resampled = np.interp(grid, np.log(wave[order]), flux[order],
                                left=np.nan, right=np.nan)
            return grid, resampled
        
    obs_grid, obs_resamp = to_log_grid(obs_wave, obs_flux)
    tmpl_grid, tmpl_resamp = to_log_grid(tmpl_wave, tmpl_flux)
    loglam_obs0, loglam_tmpl0 = obs_grid[0], tmpl_grid[0]  

    # function to normalizing flux
    def normalize(flux):
        mask = np.isfinite(flux)
        out = np.zeros_like(flux) 
        if mask.sum < 2:
            return out, mask
        std = flux[mask].std()
        if std == 0 or not np.isfinite(std):
            return out, mask
        out[mask] = (flux[mask] - flux[mask].mean())/std
        return out, mask
        
    obs_norm, obs_mask = normalize(obs_resamp)
    tmpl_norm, tmpl_mask = normalize(tmpl_resamp)
    n_obs, n_tmpl = len(obs_norm), len(tmpl_norm)

    # cross correlation by sliding the template across the spectrum, pixel by pixel
    lag_min = (np.log1p(z_min) - loglam_obs0 + loglam_tmpl0) / dloglam
    lag_max = (np.log1p(z_max) - loglam_obs0 + loglam_tmpl0) / dloglam
    lags = np.arange(int(np.floor(lag_min)), int(np.ceil(lag_max)) + 1)

    ccf = np.zeros(len(lags))
    for k, lag in enumerate(lags):
        lag = int(lag)
        i_start, i_end = max(0, lag), min(n_obs, n_tmpl + lag)
        if i_end - i_start < min_overlap:
            continue
        o = obs_norm[i_start:i_end]
        t = tmpl_norm[i_start - lag:i_end - lag]
        valid = obs_mask[i_start:i_end] & tmpl_mask[i_start - lag:i_end - lag]
        if valid.sum() < min_overlap:
            continue
        o, t = o[valid], t[valid]
        denom = np.sqrt(np.sum(o * o) * np.sum(t * t))
        if denom > 0:
            ccf[k] = np.sum(o * t) / denom

    # finding the peak, refining it to a sub-pixel precision and fiding redshift
    peak = int(np.argmax(ccf))
    best_lag = float(lags[peak])
    if 0 < peak < len(ccf) - 1:
        y0, y1, y2 = ccf[peak - 1], ccf[peak], ccf[peak + 1]
        denom = y0 - 2 * y1 + y2
        if denom != 0:
            delta = np.clip(0.5 * (y0 - y2) / denom, -1.0, 1.0)
            best_lag += delta

    z = float(np.exp(loglam_obs0 - loglam_tmpl0 + best_lag * dloglam) - 1.0)
    return z