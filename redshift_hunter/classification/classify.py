"""
Main classification interface.
"""

from .identify_lines import identify_emission_lines
from .bpt import classify_bpt
from .agn import classify_agn
from .confidence import compute_confidence


def classify_spectrum(peaks, redshift):

    identified = identify_emission_lines(peaks, redshift)

    galaxy_type, x, y = classify_bpt(identified)

    if galaxy_type == "AGN":
        galaxy_type = classify_agn(identified, y)

    confidence = compute_confidence(x, y, galaxy_type)

    return {
        "classification": galaxy_type,
        "confidence": confidence,
        "identified_lines": identified,
        "bpt_x": x,
        "bpt_y": y,
    }