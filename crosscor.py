import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import correlate, correlation_lags
from astropy.io import fits
from astropy import units as u


@dataclass
class ZCrossResult:
    """Result of a cross-correlation redshift measurement.

    Attributes
    ----------
    z : float
        Best-fit redshift.
    z_err : float
        1-sigma uncertainty on `z`, estimated via Monte Carlo resampling.
        ``nan`` if `n_mc` was 0.
    z_grid : np.ndarray
        Redshift values corresponding to each searched pixel lag.
    ccf : np.ndarray
        Normalized cross-correlation coefficient at each `z_grid` value.
    peak_index : int
        Index into `z_grid` / `ccf` of the (integer-pixel) peak.
    z_mc : np.ndarray
        Best-fit redshift recovered in each Monte Carlo trial (empty if
        `n_mc` was 0).
    """

    z: float
    z_err: float
    z_grid: np.ndarray
    ccf: np.ndarray
    peak_index: int
    z_mc: np.ndarray = field(default_factory=lambda: np.array([]))

    @property
    def peak_ccf(self) -> float:
        """Cross-correlation coefficient at the (integer-pixel) peak."""
        return float(self.ccf[self.peak_index])

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

# The peak of the cross-correlation function (CCF) is located and refined to sub-pixel precision with a parabolic fit, then converted to redshift.

def _best_z(
    obs_norm: np.ndarray,
    obs_mask: np.ndarray,
    tmpl_norm: np.ndarray,
    tmpl_mask: np.ndarray,
    lags: np.ndarray,
    dloglam: float,
    loglam_obs0: float,
    loglam_tmpl0: float,
    min_overlap: int,
) -> tuple[np.ndarray, int, float]:
    ccf = np.array(
        [
            _ccf_at_lag(obs_norm, obs_mask, tmpl_norm, tmpl_mask, int(lag), min_overlap)
            for lag in lags
        ]
    )
    peak_idx = int(np.argmax(ccf))

    # Parabolic sub-pixel refinement around the peak (if not at an edge).
    refined_lag = float(lags[peak_idx])
    if 0 < peak_idx < len(ccf) - 1:
        y0, y1, y2 = ccf[peak_idx - 1], ccf[peak_idx], ccf[peak_idx + 1]
        denom = y0 - 2 * y1 + y2
        if denom != 0:
            delta = 0.5 * (y0 - y2) / denom
            delta = float(np.clip(delta, -1.0, 1.0))
            refined_lag = float(lags[peak_idx]) + delta

    z = float(np.exp(loglam_obs0 - loglam_tmpl0 + refined_lag * dloglam) - 1.0)
    return ccf, peak_idx, z

def crosscor(
    obs_wave: np.ndarray,
    obs_flux: np.ndarray,
    tmpl_wave: np.ndarray,
    tmpl_flux: np.ndarray,
    z_min: float = -0.01,
    z_max: float = 1.0,
    dloglam: float = 1e-4,
    obs_err: np.ndarray | None = None,
    min_overlap: int = 50,
    n_mc: int = 25,
    random_state: int | np.random.Generator | None = None,
) -> ZCrossResult:
    """Estimate redshift via cross-correlation against a rest-frame template.

    Parameters
    ----------
    obs_wave, obs_flux : array_like
        Observed-frame wavelength (any units consistent with `tmpl_wave`,
        e.g. Angstrom) and flux. Must be the same length.
    tmpl_wave, tmpl_flux : array_like
        Rest-frame template wavelength and flux.
    z_min, z_max : float
        Redshift search range.
    dloglam : float
        Natural-log-wavelength pixel spacing used for resampling both
        spectra. Smaller values give finer redshift resolution at the cost
        of speed. ``1e-4`` corresponds to roughly 70 km/s pixels.
    obs_err : array_like, optional
        1-sigma flux uncertainty per pixel of the observed spectrum, same
        length as `obs_flux`. Used to scale the Monte Carlo noise. If not
        given, the RMS of `obs_flux` is used as a rough proxy.
    min_overlap : int
        Minimum number of valid overlapping pixels required to evaluate
        the CCF at a given lag; lags with less overlap are assigned ccf=0.
    n_mc : int
        Number of Monte Carlo noise realizations used to estimate `z_err`.
        Set to 0 to skip uncertainty estimation (faster).
    random_state : int or np.random.Generator, optional
        Seed or generator for reproducible Monte Carlo noise.

    Returns
    -------
    ZCrossResult
    """
    obs_wave = np.asarray(obs_wave, dtype=float)
    obs_flux = np.asarray(obs_flux, dtype=float)
    tmpl_wave = np.asarray(tmpl_wave, dtype=float)
    tmpl_flux = np.asarray(tmpl_flux, dtype=float)

    if obs_wave.shape != obs_flux.shape:
        raise ValueError("obs_wave and obs_flux must have the same shape")
    if tmpl_wave.shape != tmpl_flux.shape:
        raise ValueError("tmpl_wave and tmpl_flux must have the same shape")
    if z_max <= z_min:
        raise ValueError("z_max must be greater than z_min")

    loglam_obs0 = float(np.log(obs_wave.min()))
    loglam_tmpl0 = float(np.log(tmpl_wave.min()))

    obs_grid = _make_log_grid(obs_wave.min(), obs_wave.max(), dloglam)
    tmpl_grid = _make_log_grid(tmpl_wave.min(), tmpl_wave.max(), dloglam)

    obs_resamp = _resample_to_log(obs_wave, obs_flux, obs_grid)
    tmpl_resamp = _resample_to_log(tmpl_wave, tmpl_flux, tmpl_grid)

    obs_norm, obs_mask = _normalize(obs_resamp)
    tmpl_norm, tmpl_mask = _normalize(tmpl_resamp)

    # Lag range corresponding to the requested redshift search range.
    lag_min = (np.log1p(z_min) - loglam_obs0 + loglam_tmpl0) / dloglam
    lag_max = (np.log1p(z_max) - loglam_obs0 + loglam_tmpl0) / dloglam
    lags = np.arange(int(np.floor(lag_min)), int(np.ceil(lag_max)) + 1)
    if len(lags) == 0:
        raise ValueError("Redshift search range produced an empty lag grid")

    ccf, peak_idx, z_best = _best_z(
        obs_norm, obs_mask, tmpl_norm, tmpl_mask, lags, dloglam,
        loglam_obs0, loglam_tmpl0, min_overlap,
    )
    z_grid = np.exp(loglam_obs0 - loglam_tmpl0 + lags * dloglam) - 1.0

    z_mc = np.array([])
    z_err = float("nan")
    if n_mc > 0:
        rng = np.random.default_rng(random_state)
        if obs_err is not None:
            err = np.asarray(obs_err, dtype=float)
            err_resamp = _resample_to_log(obs_wave, err, obs_grid)
            err_resamp = np.where(np.isfinite(err_resamp), err_resamp, 0.0)
        else:
            rms = np.nanstd(obs_resamp[np.isfinite(obs_resamp)])
            err_resamp = np.full_like(obs_resamp, rms)

        z_trials = np.empty(n_mc)
        for k in range(n_mc):
            noisy = obs_resamp + rng.normal(scale=err_resamp)
            noisy_norm, noisy_mask = _normalize(noisy)
            _, _, z_trial = _best_z(
                noisy_norm, noisy_mask, tmpl_norm, tmpl_mask, lags, dloglam,
                loglam_obs0, loglam_tmpl0, min_overlap,
            )
            z_trials[k] = z_trial
        z_mc = z_trials
        z_err = float(np.std(z_trials))

    return ZCrossResult(
        z=z_best, z_err=z_err, z_grid=z_grid, ccf=ccf,
        peak_index=peak_idx, z_mc=z_mc,
    )