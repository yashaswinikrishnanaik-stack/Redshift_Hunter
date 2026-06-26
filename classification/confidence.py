"""
confidence.py

Estimate confidence of the BPT classification.
"""

from .bpt import kauffmann_curve, kewley_curve


def compute_confidence(x, y, classification):
    """
    Estimate confidence based on distance from
    the nearest classification boundary.

    Returns
    -------
    float

        Between 0 and 1
    """

    if classification == "Star-forming":

        boundary = kauffmann_curve(x)

    elif classification == "Composite":

        lower = kauffmann_curve(x)
        upper = kewley_curve(x)

        boundary = min(abs(y - lower), abs(y - upper))

        return min(1.0, 0.60 + boundary)

    else:

        boundary = kewley_curve(x)

    distance = abs(y - boundary)

    confidence = min(1.0, 0.70 + distance)

    return confidence