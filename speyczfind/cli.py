import argparse
import sys
from .spectra import Eigenspectra, observed_spectra
from .crosscorrelation import find_redshift
from .galaxyclassifier import classify_galaxy

def main():
    parser = argparse.ArgumentParser(
        description="Galaxy Analyzer: Compute redshift and classify spectral type."
    )
    
    # Required inputs
    parser.add_argument("--obs", required=True, help="Path to observed spectrum FITS file.")
    parser.add_argument("--eigen", required=True, help="Path to eigenspectra template FITS file.")
    
    # Optional tuning configurations
    parser.add_argument("--zmin", type=float, default=-0.01, help="Minimum redshift bound (default: -0.01).")
    parser.add_argument("--zmax", type=float, default=1.0, help="Maximum redshift bound (default: 1.0).")
    parser.add_argument("--dloglam", type=float, default=1e-4, help="Log-lam step sizing (default: 1e-4).")
    parser.add_argument("--overlap", type=int, default=50, help="Minimum cross-correlation pixel overlap (default: 50).")

    args = parser.parse_args()

    try:
        print(" -> Extracting spectra data...")
        obs_wave, obs_flux = observed_spectra(args.obs)
        tmpl_wave, tmpl_flux = Eigenspectra(args.eigen)

        print(" -> Running cross-correlation to determine redshift...")
        calculated_z = find_redshift(
            obs_wave, obs_flux, tmpl_wave, tmpl_flux,
            z_min=args.zmin, z_max=args.zmax, dloglam=args.dloglam, min_overlap=args.overlap
        )
        print(f" -> Calculated Redshift (z): {calculated_z:.5f}\n")

        print(" -> Performing line fitting and BPT diagnostics classification...")
        classification, log_n2_ha, log_o3_hb = classify_galaxy(args.obs, calculated_z)
        
        print("="*45)
        print(" ANALYSIS RESULTS")
        print("="*45)
        print(f"Calculated Redshift: {calculated_z:.5f}")
        if log_n2_ha is not None:
            print(f"log([NII]/Ha):       {log_n2_ha:.3f}")
            print(f"log([OIII]/Hb):      {log_o3_hb:.3f}")
        print(f"Galaxy Class:        {classification}")
        print("="*45)

    except Exception as e:
        print(f"\n[Error executing pipeline]: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
