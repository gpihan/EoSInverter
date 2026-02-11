# EoSInverter

**EoSInverter** is a small, self-contained tool for inverting tabulated equations of state (EoS).  
It computes thermodynamic variables (e.g. temperature and chemical potentials) as functions of
hydrodynamic quantities such as energy density and conserved charges.

This repository is a simplified, standalone version of a larger framework that will be available at:

- [gpihan/EoS-TrENsMUTher](https://github.com/gpihan/EoS-TrENsMUTher)

---

## Requirements

- Python 3.x  
- Python packages listed in `requirements.txt`:
  - `numpy`
  - `pandas`
  - `scipy`
  - `matplotlib`

Install them via:

```bash
pip install -r requirements.txt
```

---

## Installation

1. **Clone the repository**

```bash
git clone https://github.com/gpihan/EoSInverter.git
cd EoSInverter
```

2. **(Recommended) Create and activate a virtual environment**

```bash
python3 -m venv --clear venv
source ./venv/bin/activate       # on Linux / macOS
# .\venv\Scripts\activate        # on Windows (PowerShell / CMD)

pip3 install --upgrade pip
pip3 install -r requirements.txt
```

---

## Quick Start

1. **Place your input EoS table** (e.g. `EoS.dat` or `EoS2DISINGnew.dat`) in the working directory.

2. **Edit `parameters.py`** and adjust the entries in the `general_parameters` dictionary.

3. **Run the inversion**

```bash
python3 EoSInverter.py parameters.py
```

4. **Inspect the output directory** specified by `OutputFolder` in `parameters.py`.

---

## Configuration

All runtime options are stored in `parameters.py` inside the `general_parameters` dictionary.

### Example parameters (2D inversion, Ising-type EoS)

```python
general_parameters = {
    "RunMode": 0,  # Run Mode: 0 -> run locally
    "Dimension": 2,  # 1 -> T(e)
                     # 2 -> T(e, nB), muB(e, nB)
                     # 3 -> T(e, nB, nQ), muB(e, nB, nQ), muQ(e, nB, nQ)
                     # 4 -> T(e, nB, nQ, nS), muB, muQ, muS
    "AutoSetBoundaries": True,  # automatically set inversion boundaries from input EoS
    "NT": 200,  # grid points in tildeT direction
    "NB": 200,  # grid points in tildemuB direction
    "NQ": 2,
    "NS": 2,
    "Ttilde": [0.009777504201390019, 0.31403115649771435, 10],
    "muBtilde": [-2.2732322158914076, 2.2732322158914076, 10],
    "muQtilde": [-0.6292193138512704, 0.6292193138512704, 10],
    "muStilde": [-0.6880440931493973, 0.6880440931493973, 10],
    "Accuracy": 1e-6,
    "MAXITER": 100,
    "Number_of_cores": 12,
    "hydro_model": "vHLLE",  # or "MUSIC"
    "EoS_table": "EoS.dat",
    "OutputFolder": "invEoS",
    "TransitionLine": "TransitionlineFilePar.dat",  # (muB, T) transition line
    "RegionS": "RegionS.dat",                       # transition region S data
    "OutputMergedEoS": "invIsing-2DTExS.dat",       # merged EoS output
    "NTildeT": 100,     # points in energy density direction for merged EoS
    "NTildemuB": 100,   # points in baryon density direction for merged EoS
}
```

### Key options

- **`RunMode`**  
  Run mode selector. In this simplified version:
  - `0` – run locally.

- **`Dimension`**  
  Controls which variables are inverted; see [Inversion modes by dimension](#inversion-modes-by-dimension).

- **`AutoSetBoundaries`**  
  If `True`, boundaries for the inversion grid are determined from the input EoS table.

- **Grid resolution and ranges**
  - `NT`, `NB`, `NQ`, `NS` – number of grid points for tilde-variables.
  - `Ttilde`, `muBtilde`, `muQtilde`, `muStilde` – ranges and sampling parameters for tilde-variables.

- **Input and output files**
  - `EoS_table` – path to the input EoS table.
  - `OutputFolder` – directory where inverted tables and helper files are written.

- **Ising / merging related**
  - `TransitionLine` – path to transition line data (two-column `muB`, `T` format).
  - `RegionS` – path for the transition region S data.
  - `OutputMergedEoS` – output file name for the merged EoS table.
  - `NTildeT`, `NTildemuB` – resolution in $(e, n_B)$ space for merged EoS.

---

## Running the Code

To perform an inversion:

```bash
python3 EoSInverter.py parameters.py
```

This will:

1. Read settings from `parameters.py`.
2. Invert the provided EoS table.
3. Write results to the directory given by `OutputFolder`.

---

## Post-processing

After the main inversion finishes, you can post-process results with:

```bash
python3 PostProcessor.py converted_thermodynamics parameters.py
```

Replace `converted_thermodynamics` with the name of the `OutputFolder` you configured in `parameters.py`.

---

## Running with an Ising-model EoS

To use Ising-model based EoS tables:

1. Place the required Ising EoS files (e.g. `EoS2DISINGnew.dat`) in the working directory.
2. Update the paths in `parameters.py`.
3. Run the inversion and post-processing as described above.

### Example Ising configuration

```python
general_parameters = {
    "RunMode": 0,
    "Dimension": 2,
    "AutoSetBoundaries": True,
    "NT": 200,
    "NB": 200,
    "NQ": 2,
    "NS": 2,
    "Ttilde": [0.009777504201390019, 0.31403115649771435, 10],
    "muBtilde": [-2.2732322158914076, 2.2732322158914076, 10],
    "muQtilde": [-0.6292193138512704, 0.6292193138512704, 10],
    "muStilde": [-0.6880440931493973, 0.6880440931493973, 10],
    "Accuracy": 1e-6,
    "MAXITER": 100,
    "Number_of_cores": 12,
    "hydro_model": "vHLLE",
    "EoS_table": "EoS2DISINGnew.dat",
    "OutputFolder": "EoS2DISINGnew",
    "TransitionLine": "TransitionlineFilePar.dat",
    "RegionS": "RegionS.dat",
    "OutputMergedEoS": "invIsing-2DTExS.dat",
    "NTildeT": 100,
    "NTildemuB": 100,
}
```

After configuring, run the helper scripts **in this order**:

1. **Find the transition region S**

   ```bash
   python3 FindVals.py parameters.py
   ```

   This produces the transition region description (`RegionS.dat`).

2. **Merge the inverted EoS tables**

   ```bash
   python3 merger_eos.py parameters.py
   ```

   The merged output is written to the path set in `OutputMergedEoS`.  
   It contains the inverted thermodynamic variables as functions of energy density and
   baryon density.

   Example column format (comment header):

   ```text
   # Ttilde(GeV) muBtilde(GeV) e(GeV^4) nB(GeV^3) T(GeV) muB(GeV) P(GeV^4) s(GeV^3) cs2 chi2(GeV^2) chi
   ```

---

## Directory Structure

A typical layout of the repository:

```text
EoSInverter/
├── EoSInverter.py        # Main script (inversion)
├── PostProcessor.py      # Post-processing of inverted tables
├── FindVals.py           # Find transition region S (Ising case)
├── merger_eos.py         # Merge inverted Ising tables
├── parameters.py         # User-defined settings (general_parameters)
├── requirements.txt      # Python dependencies
├── scripts/              # Helper scripts for analysis and plotting
└── data/                 # Example or user-supplied EoS tables
```

---

## Inversion Modes by Dimension

The `Dimension` parameter controls which set of variables is inverted:

- `1` →  
  - $T(e)$
- `2` →  
  - $T(e, n_B)$  
  - $\mu_B(e, n_B)$
- `3` →  
  - $T(e, n_B, n_Q)$  
  - $\mu_B(e, n_B, n_Q)$  
  - $\mu_Q(e, n_B, n_Q)$
- `4` →  
  - $T(e, n_B, n_Q, n_S)$  
  - $\mu_B(e, n_B, n_Q, n_S)$  
  - $\mu_Q(e, n_B, n_Q, n_S)$  
  - $\mu_S(e, n_B, n_Q, n_S)$

---

## License and Citation

The full version of this project will be released as open source at:

- [gpihan/EoS-TrENsMUTher](https://github.com/gpihan/EoS-TrENsMUTher)

License and citation information for this simplified version will follow the full framework.  
If you use this code in a scientific publication before the full release, please contact the authors.
