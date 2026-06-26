# specZfind

`specZfind` is a modular Python package designed to extract astronomical spectra, calculate cosmological redshifts using cross-correlation against template eigenspectra, and classify galaxies via BPT line-ratio diagnostics.

---

## 🚀 Features

* **Spectra Extraction**: Automatically parses observed galaxy spectra and template eigenspectra from FITS files.
* **Redshift Computation**: Implements a robust cross-correlation pipeline shifting across pixel grids to determine the optimal redshift ($z$).
* **BPT Classification**: Automatically fits Gaussian profiles to major emission lines ($\text{H}\beta$, $[\text{OIII}]$, $\text{H}\alpha$, $[\text{NII}]$) and maps the results to standard empirical boundaries (Kauffmann vs. Kewley) to classify the galaxy type (Star-Forming, AGN, or Composite).

---

## 🛠️ Installation

Clone the repository, navigate into your directory, and install it locally in editable mode:

```bash
git clone [https://github.com/yashaswinikrishnanaik-stack/specZfind.git](https://github.com/yashaswinikrishnanaik-stack/specZfind.git)
cd specZfind
python3 -m pip install -e .

```

### Dependencies

The package handles installation of the following prerequisites automatically:

* `numpy`
* `scipy`
* `astropy`
* `matplotlib`

---

## 💻 Usage

The package comes equipped with a built-in Command Line Interface (CLI) pipeline that connects your data seamlessly.

### Option 1: Standard Command

If your environment's path variables are configured correctly, run the workflow globally via your terminal:

```bash
speczfind --obs 5837856538216388608.fits --eigen input.fits

```

### Option 2: Explicit Python Executable (Fallback)

If your terminal returns a `command not found` error for `speczfind`, bypass it entirely by targeting the module framework directly using your active Python interpreter:

```bash
python3 -m speczfind.cli --obs 5837856538216388608.fits --eigen input.fits

```

### Optional Arguments

You can fine-tune the cross-correlation search thresholds via the CLI:

* `--zmin`: Minimum redshift boundary search constraint (Default: `-0.01`).
* `--zmax`: Maximum redshift boundary search constraint (Default: `1.0`).
* `--dloglam`: Log-wavelength step resizing accuracy parameter (Default: `1e-4`).
* `--overlap`: Minimum required pixel overlap threshold (Default: `50`).

Example with custom configurations:

```bash
python3 -m speczfind.cli --obs 5837856538216388608.fits --eigen input.fits --zmin 0.0 --zmax 0.5

```

---

## 📁 Package Structure

```text
specZfind/
│
├── setup.py                  # Package installer script
├── README.md                 # This documentation file
└── speczfind/                # Core lowercase package module
    ├── __init__.py           # Package initializer
    ├── cli.py                # Command Line interface orchestrator
    ├── spectra.py            # FITS file parser tools
    ├── crosscorrelation.py   # Redshift calculation algorithm
    └── galaxyclassifier.py   # BPT diagnostic & line profile fitting tools

```

---

## Acknowledgements

This project is built by **Kanan**, **Yashaswini**, **Vaishnavi** and **Mywish** as part of the Code/Astro in-person workshop at Raman Research Institute, Bengaluru. We thank **Sonith LS** (the local TA) and all the organizers of the Code/Astro fraternity for this learning and building opportunity.

```
