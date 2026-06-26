import os
import sys

# Crucial step: Ensure Python looks at the current directory for 'speczfind'
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from speczfind.spectra import Eigenspectra, observed_spectra
from speczfind.crosscorrelation import find_redshift
from speczfind.galaxyclassifier import classify_galaxy

def run_manual_test():
    # ---------------------------------------------------------
    # DEFINE YOUR INPUT FILES HERE MANUALLY
    # ---------------------------------------------------------
    observed_fits = "spec-1678-53433-0425.fits"
    eigenspectra_fits = "template_galaxy.fits"
    
    print(f"Checking for files...")
    if not os.path.exists(observed_fits) or not os.path.exists(eigenspectra_fits):
        print(f"Error: Make sure '{observed_fits}' and '{eigenspectra_fits}' are in this same folder!")
        return

    try:
        print("-> Step 1: Extracting spectra data...")
        obs_wave, obs_flux = observed_spectra(observed_fits)
        tmpl_wave, tmpl_flux = Eigenspectra(eigenspectra_fits)

        print("-> Step 2: Running cross-correlation to determine redshift...")
        calculated_z = find_redshift(
            obs_wave, obs_flux, tmpl_wave, tmpl_flux,
            z_min=-0.01, z_max=1.0, dloglam=1e-4, min_overlap=50
        )
        print(f"   Calculated Redshift (z): {calculated_z:.5f}\n")

        print("-> Step 3: Performing line fitting and BPT diagnostics...")
        classification, log_n2_ha, log_o3_hb = classify_galaxy(observed_fits, calculated_z)
        
        print("\n" + "="*45)
        print(" MANUAL ANALYSIS RESULTS")
        print("="*45)
        print(f"Observed File:   {observed_fits}")
        print(f"Template File:   {eigenspectra_fits}")
        print(f"Redshift (z):    {calculated_z:.5f}")
        if log_n2_ha is not None:
            print(f"log([NII]/Ha):   {log_n2_ha:.3f}")
            print(f"log([OIII]/Hb):  {log_o3_hb:.3f}")
        print(f"Galaxy Class:    {classification}")
        print("="*45)

    except Exception as e:
        print(f"\n[Execution Failed]: {e}")

if __name__ == "__main__":
    run_manual_test()
