"""
classify.py

Main interface for galaxy classification.
"""

from .identify_lines import identify_emission_lines
from .bpt import classify_bpt
from .confidence import compute_confidence


def classify_spectrum(peaks, redshift):
    """
    Main classification function.

    Parameters
    ----------
    peaks : list
        List of tuples

        (observed_wavelength, flux)

    redshift : float

    Returns
    -------
    dict
    """

    # Step 1
    identified_lines = identify_emission_lines(peaks, redshift)

    # Step 2
    classification, x, y = classify_bpt(identified_lines)

    # Step 3
    confidence = compute_confidence(x, y, classification)

    return {
        "classification": classification,
        "confidence": confidence,
        "identified_lines": identified_lines,
        "bpt_x": x,
        "bpt_y": y,
    }