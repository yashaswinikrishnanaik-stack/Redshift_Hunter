import numpy as np
from astropy.io import fits
from scipy.signal import medfilt

def Eigenspectra(filename):
    """Extract the template spectra from SDSS eigenspectra."""
    with fits.open(filename) as hdul:
        E1 = hdul[0].header["EIGEN0"]
        E2 = hdul[0].header["EIGEN1"]
        E3 = hdul[0].header["EIGEN2"]
        E4 = hdul[0].header["EIGEN3"]
        data_full = (hdul[0].data[0, :] * E1 + 
                     hdul[0].data[1, :] * E2 + 
                     hdul[0].data[2, :] * E3 + 
                     hdul[0].data[3, :] * E4)
        spec = data_full.astype(np.float64)
        data_continuum = medfilt(spec, kernel_size=301)
        eig_lines = data_full - data_continuum
        
        coeff0 = hdul[0].header["COEFF0"]
        coeff1 = hdul[0].header["COEFF1"]
        pix = np.arange(len(eig_lines))
        loglam = coeff0 + coeff1 * pix
        wavelength = 10**loglam   
    return wavelength, eig_lines

def observed_spectra(filename):
    """Extract the observed flux versus wavelength data."""
    with fits.open(filename) as hdul: # Fixed hardcoded string bug here
        spec = hdul[1].data
    flux = spec["flux"]               
    wavelength = 10**spec["loglam"]   
    return wavelength, flux
