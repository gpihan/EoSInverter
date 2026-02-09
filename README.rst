
=======================
EoSInverter
=======================

**EoSInverter** is a small, self-contained tool for inverting
tabulated equations of state (EoS). It expresses thermodynamic
variables such as temperature and chemical potentials as functions of
hydrodynamic quantities like energy density and conserved charges.

This repository is a simplified, standalone version of the full
framework that will be available at:

   https://github.com/gpihan/EoS-TrENsMUTher

-----------------------
Requirements
-----------------------

- Python 3 ("python3" on most systems)
- Python packages listed in ``requirements.txt``:

  - numpy
  - pandas
  - scipy
  - matplotlib

-----------------------
Installation
-----------------------

1. Clone the repository

.. code-block:: bash

   git clone https://github.com/gpihan/EoSInverter.git
   cd EoSInverter

2. Create and activate a virtual environment (recommended)

.. code-block:: bash

   python3 -m venv --clear venv
   source ./venv/bin/activate
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

Configuration
-----------------------

Edit ``parameters.py`` to configure the run. Parameters are stored in
the dictionary ``general_parameters``.

Example of default parameters (2D inversion, Ising-type EoS):

.. code-block:: python

   general_parameters = {
       "RunMode": 0, # Run Mode: 0 -> Run locally
       "Dimension": 2, # Dimension of inversion: 1 -> T(e), 2 -> T(e, nb), muB(e, nb), 3 -> T(e, nb, nq), muB(e, nb, nq), muQ(e, nb, nq), 4 -> T(e, nb, nq, ns), muB(e, nb, nq, ns), muQ(e, nb, nq, ns), muS(e, nb, nq, ns)
       "AutoSetBoundaries": True, # Automatically set boundaries for inversion based on input EoS table
       "NT": 200, # Number of grid points in tildeT direction for inversion
       "NB": 200, # Number of grid points in tildemuB direction for inversion
       "NQ": 2,
       "NS": 2,
       "Ttilde": [0.009777504201390019, 0.31403115649771435, 10],
       "muBtilde": [-2.2732322158914076, 2.2732322158914076, 10],
       "muQtilde": [-0.6292193138512704, 0.6292193138512704, 10],
       "muStilde": [-0.6880440931493973, 0.6880440931493973, 10],
       "Accuracy": 1e-6, 
       "MAXITER": 100,
       "Number_of_cores": 12,
       "hydro_model": "vHLLE", # MUSIC
       "EoS_table": "EoS.dat",
       "OutputFolder": "invEoS",
       "EoS_table": "EoS2DISINGnew.dat",
       "OutputFolder": "EoS2DISINGnew",
       "TransitionLine": "TransitionlineFilePar.dat",  # Path to transition line data (muB, T format)
       "RegionS": "RegionS.dat",  # path to transition region S data
       "OutputMergedEoS": "invIsing-2DTExS.dat",  # output merged EoS table
       "NTildeT": 100,  # number of points in energy density direction for merged EoS
       "NTildemuB": 100,  # number of points in baryon density direction for merged ES
      }

   Some important options:

   - ``RunMode`` – run mode selector. In this simplified version,
     ``0`` means running locally.
   - ``Dimension`` – controls which variables are inverted; see the
     section "Inversion modes by dimension" below.
   - ``AutoSetBoundaries`` – if ``True``, boundaries for the inversion
     grid are determined from the input EoS table.
   - ``NT``, ``NB``, ``NQ``, ``NS`` – number of grid points in the
     corresponding tilde-variables used for the inversion.
   - ``Ttilde``, ``muBtilde``, ``muQtilde``, ``muStilde`` – ranges and
     sampling parameters for the tilde-variables.
   - ``EoS_table`` – path to the input EoS table.
   - ``OutputFolder`` – directory where the inverted tables and helper
     files are written.
   - ``TransitionLine``, ``RegionS``, ``OutputMergedEoS``, ``NTildeT``,
     ``NTildemuB`` – parameters used when working with the Ising-model
     EoS and when merging several inverted tables into a single one.

   -----------------------
   Running the code
   -----------------------

To run the inversion:

.. code-block:: bash

   python3 EoSInverter.py parameters.py

This reads settings from ``parameters.py``, inverts the provided EoS
table and writes results to the directory specified by ``OutputFolder``.

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

To run with Ising-model based tables, place the appropriate EoS files
in the working directory and update the paths in ``parameters.py``.
Run the inversion and post-processing as above.

To merge the three inverted tables into a single merged EoS, set the
merger-related keys in ``general_parameters``. Example:

.. code-block:: python

   general_parameters = {
       "RunMode": 0, # Run Mode: 0 -> Run locally
       "Dimension": 2, # Dimension of inversion: 1 -> T(e), 2 -> T(e, nb), muB(e, nb), 3 -> T(e, nb, nq), muB(e, nb, nq), muQ(e, nb, nq), 4 -> T(e, nb, nq, ns), muB(e, nb, nq, ns), muQ(e, nb, nq, ns), muS(e, nb, nq, ns)
       "AutoSetBoundaries": True, # Automatically set boundaries for inversion based on input EoS table
       "NT": 200, # Number of grid points in tildeT direction for inversion
       "NB": 200, # Number of grid points in tildemuB direction for inversion
       "NQ": 2,
       "NS": 2,
       "Ttilde": [0.009777504201390019, 0.31403115649771435, 10],
       "muBtilde": [-2.2732322158914076, 2.2732322158914076, 10],
       "muQtilde": [-0.6292193138512704, 0.6292193138512704, 10],
       "muStilde": [-0.6880440931493973, 0.6880440931493973, 10],
       "Accuracy": 1e-6, 
       "MAXITER": 100,
       "Number_of_cores": 12,
       "hydro_model": "vHLLE", # MUSIC
       "EoS_table": "EoS.dat",
       "OutputFolder": "invEoS",
       "EoS_table": "EoS2DISINGnew.dat",
       "OutputFolder": "EoS2DISINGnew",
       "TransitionLine": "TransitionlineFilePar.dat",  # Path to transition line data (muB, T format)
       "RegionS": "RegionS.dat",  # path to transition region S data
       "OutputMergedEoS": "invIsing-2DTExS.dat",  # output merged EoS table
       "NTildeT": 100,  # number of points in energy density direction for merged EoS
       "NTildemuB": 100,  # number of points in baryon density direction for merged ES
   }

After configuring, run the helper scripts in this order:

1. Find the transition region S

.. code-block:: bash

   python3 FindVals.py parameters.py

This produces data describing the transition region (``RegionS.dat``).

2. Merge the inverted EoS tables

.. code-block:: bash

   python3 merger_eos.py parameters.py

The merged output will be written to the path set in ``OutputMergedEoS``
and contains the inverted thermodynamic variables as functions of energy
density and baryon density.

Example column format (comment header):

.. code-block:: text

   # Ttilde(GeV) muBtilde(GeV) e(GeV^4) nB(GeV^3) T(GeV) muB(GeV) P(GeV^4) s(GeV^3) cs2 chi2(GeV^2) chi


Directory structure
-----------------------

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

-----------------------
Inversion modes by dimension
-----------------------

Depending on the value of ``Dimension``, the inversion computes different
sets of variables:

- ``1`` → T(e)
- ``2`` → T(e, nB), muB(e, nB)
- ``3`` → T(e, nB, nQ), muB(e, nB, nQ), muQ(e, nB, nQ)
- ``4`` → T(e, nB, nQ, nS), muB, muQ, muS as functions of
  (e, nB, nQ, nS)

-----------------------
License and citation
-----------------------

The full version of this project will be released under an open-source
license at:

    https://github.com/gpihan/EoS-TrENsMUTher

