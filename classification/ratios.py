"""
ratios.py

Compute emission-line ratios used in BPT diagnostics.
"""

import math


def compute_bpt_ratios(lines):
    """
    Compute the [NII]-based BPT ratios.

    Returns
    -------
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
            raise ValueError(f"Missing emission line: {line}")

    x = math.log10(lines["NII6584"] / lines["Halpha"])
    y = math.log10(lines["OIII5007"] / lines["Hbeta"])

    return x, y


def compute_sii_ratio(lines):
    """
    Compute log10(([SII]6716+[SII]6731)/Halpha)
    """

    required = [
        "Halpha",
        "SII6716",
        "SII6731",
    ]

    for line in required:
        if line not in lines:
            raise ValueError(f"Missing emission line: {line}")

    sii = lines["SII6716"] + lines["SII6731"]

    return math.log10(sii / lines["Halpha"])


def compute_oi_ratio(lines):
    """
    Compute log10([OI]6300/Halpha)
    """

    required = [
        "Halpha",
        "OI6300",
    ]

    for line in required:
        if line not in lines:
            raise ValueError(f"Missing emission line: {line}")

    return math.log10(lines["OI6300"] / lines["Halpha"])