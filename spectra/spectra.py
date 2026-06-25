import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import medfilt
#################################################################


def Eigenspectra(filename):
    """Eigenspectra
    
    extract the template spectra from sdss eigenspectra .

    Args:
        filename(string): fits file containing eigenspectra of sdss galaxy 
    Returns:
        list: wavelength, continumm subtracted template spectra
    """
    with fits.open(filename) as hdul:
        E1 = hdul[0].header["EIGEN0"]
        E2 = hdul[0].header["EIGEN1"]
        E3 = hdul[0].header["EIGEN2"]
        E4 = hdul[0].header["EIGEN3"]
        data_full = hdul[0].data[0,:]*E1 + hdul[0].data[1,:]*E2+hdul[0].data[2,:]*E3+hdul[0].data[3,:]*E4
        spec = data_full.astype(np.float64)
        data_continumm = medfilt(spec, kernel_size=301)
        eig_lines = data_full - data_continumm #eigenspectra flux
    coeff0 = hdul[0].header["COEFF0"]
    coeff1 = hdul[0].header["COEFF1"]
    pix = np.arange(7012)
    loglam = coeff0 + coeff1 * pix
    wavelength = 10**loglam   # eigenspectra wavelength
    return([wavelength,eig_lines])
#####################################################################
def observed_spectra(filename):
    """Observed spectra
    
    Extract the observed flux versus wavelength data

    Args:
        filename(string): fits file containing the observed spectra of a galaxy

    Returns:
        list: wavelength, flux 

    """

    with fits.open("spec-1678-53433-0425.fits") as hdul:
        spec = hdul[1].data

    flux = spec["flux"]               #flux of the observed galaxy
    wavelength = 10**spec["loglam"]   #wavelength axis of observed galaxy
    return([wavelength,flux])
