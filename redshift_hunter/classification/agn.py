"""
AGN subtype classification.

Uses the [SII] diagnostic diagram to distinguish
between Seyfert galaxies and LINERs.
"""

from .ratios import compute_sii_ratio


def seyfert_liner_boundary(x):
    """
    Kewley et al. (2006)

    Boundary separating Seyfert and LINER.
    """

    return 1.89 * x + 0.76


def classify_agn(lines, oiii_hbeta):
    """
    Parameters
    ----------
    lines : dict

    oiii_hbeta : float
        log10([OIII]/Hbeta)

    Returns
    -------
    str

        Seyfert
        or
        LINER
    """

    sii_ratio = compute_sii_ratio(lines)

    boundary = seyfert_liner_boundary(sii_ratio)

    if oiii_hbeta >= boundary:
        return "Seyfert"

    return "LINER"