import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit

LINES = {
    'H_beta':  4861.3,
    'OIII_5007': 5006.8,
    'H_alpha': 6562.8,
    'NII_6584': 6583.4
}

def gaussian(x, amplitude, mean, sigma, continuum):
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * sigma ** 2)) + continuum

def fit_emission_line(wave, flux, rest_wave, window=30.0):
    mask = (wave >= rest_wave - window) & (wave <= rest_wave + window)
    x = wave[mask]
    y = flux[mask]
    if len(x) < 5:
        return 0.0  

    cont_guess = np.median(y)
    amp_guess = np.max(y) - cont_guess
    p0 = [amp_guess, rest_wave, 3.0, cont_guess]
    bounds = ([0, rest_wave - 5, 0.5, -np.inf], [np.inf, rest_wave + 5, 15.0, np.inf])

    try:
        popt, _ = curve_fit(gaussian, x, y, p0=p0, bounds=bounds, maxfev=5000)
        amp, mean, sigma, cont = popt
        return amp * np.abs(sigma) * np.sqrt(2 * np.pi)
    except RuntimeError:
        return 0.0

def classify_galaxy(fits_path, redshift):
    with fits.open(fits_path) as hdul:
        data = hdul[1].data
        flux = data['flux']
        wave_obs = 10**data['loglam'] if 'loglam' in data.names else data['wavelength']

    wave_rest = wave_obs / (1.0 + redshift)
    fluxes = {line_name: fit_emission_line(wave_rest, flux, rest_w) for line_name, rest_w in LINES.items()}

    hb, o3, ha, n2 = fluxes['H_beta'], fluxes['OIII_5007'], fluxes['H_alpha'], fluxes['NII_6584']
    if min(hb, o3, ha, n2) <= 0:
        return "Unclassified / Weak Emission Lines (Possible Passive Galaxy)", None, None

    log_n2_ha = np.log10(n2 / ha)
    log_o3_hb = np.log10(o3 / hb)

    k03_boundary = 0.61 / (log_n2_ha - 0.05) + 1.3 if log_n2_ha < 0.05 else -np.inf
    k01_boundary = 0.61 / (log_n2_ha - 0.47) + 1.19 if log_n2_ha < 0.47 else -np.inf

    if log_n2_ha < 0.05 and log_o3_hb < k03_boundary:
        classification = "Star-Forming Galaxy"
    elif log_n2_ha >= 0.47 or log_o3_hb > k01_boundary:
        classification = "Active Galactic Nucleus (AGN)"
    else:
        classification = "Composite Galaxy (Mixed SF/AGN Activity)"

    return classification, log_n2_ha, log_o3_hb
