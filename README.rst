EoSInverter
===========

**EoSInverter** is a small, self-contained tool for inverting tabulated
equations of state (EoS). It computes thermodynamic variables (such as
temperature and chemical potentials) as functions of hydrodynamic
quantities like energy density and conserved charges.

This repository is a simplified, standalone version of the full
framework that will be available at:

-  https://github.com/gpihan/EoS-TrENsMUTher

-----------------------
Requirements
-----------------------

- Python 3.x
- Python packages listed in ``requirements.txt``:

  - numpy
  - pandas
  - scipy
  - matplotlib

Install the dependencies with:

.. code-block:: bash

   pip install -r requirements.txt

-----------------------
Installation
-----------------------

1. Clone the repository

   .. code-block:: bash

      git clone https://github.com/gpihan/EoSInverter.git
      cd EoSInverter

2. (Recommended) Create and activate a virtual environment

   .. code-block:: bash

      python3 -m venv --clear venv
      source ./venv/bin/activate        # Linux / macOS
      # .\venv\Scripts\activate         # Windows (PowerShell / CMD)

      pip3 install --upgrade pip
      pip3 install -r requirements.txt

-----------------------
Quick start
-----------------------

1. Place your input EoS table (for example ``EoS.dat`` or
   ``EoS2DISINGnew.dat``) in the working directory.

2. Edit ``parameters.py`` and adjust the entries in the
   ``general_parameters`` dictionary.

3. Run the inversion:

   .. code-block:: bash

      python3 EoSInverter.py parameters.py

4. Inspect the output directory given by ``OutputFolder`` in
   ``parameters.py``.

-----------------------
Configuration
-----------------------

Edit ``parameters.py`` to configure the run. Parameters are stored in
the dictionary ``general_parameters``.

Example of default parameters (2D inversion, Ising-type EoS):

.. code-block:: python

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

Some important options:

- ``RunMode`` – run mode selector. In this simplified version:

  - ``0`` – run locally.

- ``Dimension`` – controls which variables are inverted; see
  section :ref:`Inversion modes by dimension <dim-modes>`.

- ``AutoSetBoundaries`` – if ``True``, boundaries for the inversion
  grid are determined from the input EoS table.

- ``NT``, ``NB``, ``NQ``, ``NS`` – number of grid points in the
  corresponding tilde-variables.

- ``Ttilde``, ``muBtilde``, ``muQtilde``, ``muStilde`` – ranges and
  sampling parameters for the tilde-variables.

- ``EoS_table`` – path to the input EoS table.

- ``OutputFolder`` – directory where the inverted tables and helper
  files are written.

- ``TransitionLine``, ``RegionS``, ``OutputMergedEoS``,
  ``NTildeT``, ``NTildemuB`` – parameters used when working with the
  Ising-model EoS and when merging several inverted tables into a
  single one.

-----------------------
Running the code
-----------------------

To run the inversion:

.. code-block:: bash

   python3 EoSInverter.py parameters.py

This will:

1. Read settings from ``parameters.py``.
2. Invert the provided EoS table.
3. Write results to the directory specified by ``OutputFolder``.

-----------------------
Post-processing
-----------------------

After the main script finishes you can post-process results with:

.. code-block:: bash

   python3 PostProcessor.py converted_thermodynamics parameters.py

Replace ``converted_thermodynamics`` with the ``OutputFolder`` you set
in ``parameters.py``.

-----------------------
Running with an Ising-model EoS
-----------------------

To run with Ising-model based tables, place the appropriate EoS files in
the working directory and update the paths in ``parameters.py``. Run
the inversion and post-processing as above.

Example Ising configuration:

.. code-block:: python

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

After configuring, run the helper scripts in this order:

1. Find the transition region S

   .. code-block:: bash

      python3 FindVals.py parameters.py

   This produces data describing the transition region
   (``RegionS.dat``).

2. Merge the inverted EoS tables

   .. code-block:: bash

      python3 merger_eos.py parameters.py

   The merged output will be written to the path set in
   ``OutputMergedEoS`` and contains the inverted thermodynamic
   variables as functions of energy density and baryon density.

   Example column format (comment header):

   .. code-block:: text

      # Ttilde(GeV) muBtilde(GeV) e(GeV^4) nB(GeV^3) T(GeV) muB(GeV) P(GeV^4) s(GeV^3) cs2 chi2(GeV^2) chi

-----------------------
Directory structure
-----------------------

Typical repository layout:

.. code-block:: text

   EoSInverter/
   ├── EoSInverter.py        # Main script (inversion)
   ├── PostProcessor.py      # Post-processing of inverted tables
   ├── FindVals.py           # Find transition region S (Ising case)
   ├── merger_eos.py         # Merge inverted Ising tables
   ├── parameters.py         # User-defined settings (general_parameters)
   ├── requirements.txt      # Python dependencies
   ├── scripts/              # Helper scripts for analysis and plotting
   └── data/                 # Example or user-supplied EoS tables

.. _dim-modes:

-----------------------
Inversion modes by dimension
-----------------------

Depending on the value of ``Dimension``, the inversion computes different
sets of variables:

- ``1`` → T(e)
- ``2`` → T(e, nB), muB(e, nB)
- ``3`` → T(e, nB, nQ), muB(e, nB, nQ), muQ(e, nB, nQ)
- ``4`` → T(e, nB, nQ, nS), muB(e, nB, nQ, nS),
  muQ(e, nB, nQ, nS), muS(e, nB, nQ, nS)

-----------------------
License and citation
-----------------------

The full version of this project will be released under an open-source
license at:

- https://github.com/gpihan/EoS-TrENsMUTher

License and citation information for this simplified version will follow
the full framework. If you use this code in a scientific publication
before the full release, please contact the authors.
