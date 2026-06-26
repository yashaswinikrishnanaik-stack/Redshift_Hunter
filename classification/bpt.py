"""
bpt.py

Implements the standard BPT classification using the
Kauffmann (2003) and Kewley (2001) diagnostic curves.
"""

from .ratios import compute_bpt_ratios


def kauffmann_curve(x):
    """
    Kauffmann et al. (2003)
    """

    return (0.61 / (x - 0.05)) + 1.30


def kewley_curve(x):
    """
    Kewley et al. (2001)
    """

    return (0.61 / (x - 0.47)) + 1.19


def classify_bpt(lines):
    """
    Classify a galaxy using the [NII] BPT diagram.

    Parameters
    ----------
    lines : dict

    Returns
    -------
    tuple

        classification

        x

        y
    """

    x, y = compute_bpt_ratios(lines)

    if y < kauffmann_curve(x):

        classification = "Star-forming"

    elif y < kewley_curve(x):

        classification = "Composite"

    else:

        classification = "AGN"

    return classification, x, y