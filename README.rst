
=======================
EoSInverter
=======================

**EoSInverter** is a tool for inverting tabulated equations of state (EoS) to express thermodynamic variables—such as temperature and chemical potentials—as functions of hydrodynamic quantities like energy density and conserved charges.

This repository is a simplified version of the full framework that will be available at:

   https://github.com/gpihan/EoS-TrENsMUTher

-----------------------
Installation
-----------------------

1. Clone the repository

.. code-block:: bash

   git clone https://github.com/gpihan/EoSInverter.git
   cd EoSInverter

2. Create and activate a virtual environment

.. code-block:: bash

   python3 -m venv --clear venv
   source ./venv/bin/activate
   pip3 install --upgrade pip
   pip3 install -r requirements.txt

-----------------------
Configuration
-----------------------

Edit `parameters.py` to configure the run. Parameters are stored in the dictionary
`general_parameters`.

Example of default parameters:

.. code-block:: python

   general_parameters = {
       "RunMode": 0,  # 0 = run locally, 1 = run on Slurm cluster
       "Dimension": 2,  # 2D inversion: T(e, nb), muB(T, nb)
       "AutoSetBoundaries": True,
       "NT": 100,
       "NB": 100,
       "NQ": 2,
       "NS": 2,
       "Ttilde": [0.0098, 0.3140, 10],
       "muBtilde": [-2.2732, 2.2732, 10],
       "muQtilde": [-0.6292, 0.6292, 10],
       "muStilde": [-0.6880, 0.6880, 10],
       "Accuracy": 1e-6,
       "MAXITER": 100,
       "Number_of_cores": 12,
       "hydro_model": "vHLLE",  # or "MUSIC"
       "EoS_table": "converted_thermodynamics.dat",
       "OutputFolder": "converted_thermodynamics",
   }

-----------------------
Running the code
-----------------------

To run the inversion:

.. code-block:: bash

   python3 EoSInverter.py parameters.py

This reads settings from `parameters.py`, inverts the provided EoS table and writes results to
the directory specified by `OutputFolder`.

-----------------------
Post-processing
-----------------------

After the main script finishes you can post-process results with:

.. code-block:: bash

   python3 PostProcessor.py converted_thermodynamics parameters.py

Replace `converted_thermodynamics` with the `OutputFolder` you set in `parameters.py`.

-----------------------
Running with an Ising-model EoS
-----------------------

To run with Ising-model based tables, place the appropriate EoS files (for example
`EoS_Above.dat`, `EoS_Below.dat`, `EoS_Cross.dat`) in the working directory and update the paths
in `parameters.py`. Run the inversion and post-processing as above.

To merge the three inverted tables into a single merged EoS, set the merger-related keys in
`general_parameters`. Example:

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
       "EoS_table": "EoSCrossGeV.dat",
       "OutputFolder": "EoSCrossGeV_test",
       # Paths to EoS files for merger
       "EoS_above_input": "EoSAboveGeV.dat",
       "EoS_below_input": "EoSBelowGeV.dat",
       "EoS_cross_input": "EoSCrossGeV.dat",
       "premerger_eos": True,
       "premerger_output": "output",
       "EoS_above": "EoS_inv_above.dat",
       "EoS_below": "EoS_inv_below.dat",
       "EoS_cross": "EoS_inv_cross.dat",
       "TransitionLine": "TransitionLine.dat",
       "RegionS": "RegionS.dat",
       "OutputMergedEoS": "invMergedEoS_new.dat",
       "Ne": 400,
       "Nn": 400,
   }

After configuring, run the helper scripts in this order:

1. Find the transition region S:

.. code-block:: bash

   python3 FindVals.py parameters.py

This produces data describing the transition region (`RegionS.dat`).

2. Preprocess (clean) data before merging:

.. code-block:: bash

   python3 premerger_eos.py parameters.py

3. Merge the inverted EoS tables:

.. code-block:: bash

   python3 merger_eos.py parameters.py

The merged output will be written to the path set in `OutputMergedEoS` and contains the
inverted thermodynamic variables as functions of energy density and baryon density.

Example column format (comment header):

.. code-block:: text

   # e nb Ttilde muBtilde T muB P S chi chi_e chi_n

-----------------------
Directory structure
-----------------------

.. code-block:: text

   EoSInverter/
   ├── EoSInverter.py        # Main script
   ├── PostProcessor.py      # Post-processing
   ├── parameters.py         # User-defined settings
   ├── requirements.txt      # Python dependencies
   └── ...

-----------------------
Inversion modes by dimension
-----------------------

Depending on the value of `Dimension`, the inversion computes different sets of variables:

- `1` → T(e)
- `2` → T(e, nb), muB(e, nb)
- `3` → T(e, nb, nq), muB(e, nb, nq), muQ(e, nb, nq)
- `4` → T(e, nb, nq, ns), muB, muQ, muS as functions of (e, nb, nq, ns)

-----------------------
License and citation
-----------------------

The full version of this project will be released under an open-source license at:

    https://github.com/gpihan/EoS-TrENsMUTher

