=======================
EoSInverter
=======================

**EoSInverter** is a tool for inverting tabulated equations of state (EoS) to express thermodynamic variables—such as temperature and chemical potentials—as functions of hydrodynamic quantities like energy density and conserved charges.

This repository is a simplified version of the full framework that will be available at:

    https://github.com/gpihan/EoS-TrENsMUTher

-----------------------
Installation
-----------------------

1. **Clone the repository**:

.. code-block:: bash

   git clone https://github.com/gpihan/EoSInverter.git
   cd EoSInverter

2. **Create and activate a virtual environment**:

.. code-block:: bash

   python3 -m venv --clear venv
   source ./venv/bin/activate
   pip3 install --upgrade pip
   pip3 install -r requirements.txt    

-----------------------
Configuration
-----------------------

Edit the file ``parameters.py`` to configure your run. The parameters are stored in a dictionary called ``general_parameters``.

Here’s a breakdown of the default parameters provided:

.. code-block:: python

   general_parameters = {
       "RunMode": 0,  # 0 = run locally, 1 = run on Slurm cluster
       "Dimension": 2,  # 2D inversion: T(e, nb), μB(T, nb)
       "AutoSetBoundaries": True,  # Automatically infer tilde variable ranges from input EoS table
       "NT": 100,  # Number of temperature grid points
       "NB": 100,  # Number of baryon density grid points
       "NQ": 2,    # Grid points in electric charge direction (ignored if Dimension < 3)
       "NS": 2,    # Grid points in strangeness direction (ignored if Dimension < 4)
       "Ttilde": [0.0098, 0.3140, 10],  # Manual override for temperature tilde range
       "muBtilde": [-2.2732, 2.2732, 10],  # Manual override for μB tilde range
       "muQtilde": [-0.6292, 0.6292, 10],  # Manual override for μQ tilde range
       "muStilde": [-0.6880, 0.6880, 10],  # Manual override for μS tilde range
       "Accuracy": 1e-6,  # Convergence criterion for root-finding
       "MAXITER": 100,  # Max iterations allowed for solvers
       "Number_of_cores": 12,  # Number of CPU cores for parallel cluster execution
       "hydro_model": "vHLLE",  # Choose between "vHLLE" and "MUSIC" hydro models -> different output formats
       "EoS_table": "converted_thermodynamics.dat",  # Path to input EoS table
       "OutputFolder": "converted_thermodynamics",   # Folder to save results
   }

-----------------------
Running the Code
-----------------------

To start the inversion process, use:

.. code-block:: bash

   python3 EoSInverter.py parameters.py

This reads the settings from `parameters.py`, inverts the EoS table, and saves the output to the specified folder.

-----------------------
Post-Processing
-----------------------

After running the main script, you can post-process the results using:

.. code-block:: bash

   python3 PostProcessor.py converted_thermodynamics parameters.py

Replace ``converted_thermodynamics`` with the output folder specified in your config if different.

-----------------------
Directory Structure
-----------------------

.. code-block:: text

   EoSInverter/
   ├── EoSInverter.py        # Main script
   ├── PostProcessor.py      # Post-processing
   ├── parameters.py         # User-defined settings
   ├── requirements.txt      # Python dependencies
   └── ...

-----------------------
Inversion Modes by Dimension
-----------------------

Depending on the ``Dimension`` value, the inversion behavior changes:

- ``1`` → T(e)
- ``2`` → T(e, nb), μB(T, nb)
- ``3`` → T(e, nb, nq), μB(e, nb, nq), μQ(e, nb, nq)
- ``4`` → T(e, nb, nq, ns), μB, μQ, μS as functions of (e, nb, nq, ns)

-----------------------
License and Citation
-----------------------

The full version of this project will be released under an open-source license at:

    https://github.com/gpihan/EoS-TrENsMUTher

