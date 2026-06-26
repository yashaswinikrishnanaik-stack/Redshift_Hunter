"""Functions for identifying emission lines from observed peaks."""

from .constants import EMISSION_LINES, DEFAULT_TOLERANCE


def observed_to_rest(observed_wavelength, redshift):
    """Convert observed wavelength to rest-frame wavelength."""
    return observed_wavelength / (1 + redshift)


def identify_emission_lines(peaks, redshift, tolerance=DEFAULT_TOLERANCE):
    """Identify emission lines from detected peaks."""
    identified = {}

    for observed_wave, flux in peaks:
        rest_wave = observed_to_rest(observed_wave, redshift)
        best_line = None
        smallest_difference = float("inf")

        for line_name, line_wave in EMISSION_LINES.items():
            difference = abs(rest_wave - line_wave)
            if difference < smallest_difference:
                smallest_difference = difference
                best_line = line_name

        if smallest_difference <= tolerance:
            identified[best_line] = flux

    return identified
