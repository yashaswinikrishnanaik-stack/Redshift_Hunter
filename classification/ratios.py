"""
ratios.py

Functions to compute emission-line ratios used in BPT diagrams.
"""

import math


def compute_bpt_ratios(lines):
    """
    Compute the BPT diagnostic ratios.

    Parameters
    ----------
    lines : dict
        Dictionary containing identified emission line fluxes.

        Example
        -------
        {
            "Halpha": 120,
            "Hbeta": 40,
            "OIII5007": 95,
            "NII6584": 55
        }

    Returns
    -------
    tuple
        (x, y)

        x = log10([NII]6584 / Halpha)

        y = log10([OIII]5007 / Hbeta)
    """

    required = [
        "Halpha",
        "Hbeta",
        "OIII5007",
        "NII6584",
    ]

    for line in required:
        if line not in lines:
            raise ValueError(f"Missing required emission line: {line}")

    x = math.log10(lines["NII6584"] / lines["Halpha"])
    y = math.log10(lines["OIII5007"] / lines["Hbeta"])

    return x, y