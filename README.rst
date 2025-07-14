=======================
EoSInverter
=======================

**EoSInverter** is a simplified implementation of the forthcoming
`EoS-TrENsMUTher <https://github.com/gpihan/EoS-TrENsMUTher>`_, designed for
inverting tabulated Equations of State (EoS) in high-energy nuclear physics.

The code reconstructs thermodynamic variables—such as temperature and
chemical potentials—as functions of hydrodynamic quantities like energy
density and conserved charges.

-----------------------
Installation
-----------------------

1. **Clone the repository**:

.. code-block:: bash

   git clone https://github.com/gpihan/EoSInverter.git
   cd EoSInverter

2. **Create a virtual environment and install dependencies**:

.. code-block:: bash

   python3 -m venv --clear venv
   ./venv/bin/pip install --upgrade pip
   ./venv/bin/pip install -r requirements.txt

-----------------------
Configuration
-----------------------

Edit the ``parameters.py`` file to set up your run. All parameters are defined in a Python dictionary format.

**Key Parameters Overview**:

- **RunMode**:
  
  - ``0``: Run locally using one CPU
  - ``1``: Run on a Slurm cluster

- **Dimension**:

  - ``1``: T(e)
  - ``2``: T(e, nb), μB(T, nb)
  - ``3``: T(e, nb, nq), μB(e, nb, nq), μQ(e, nb, nq)
  - ``4``: T(e, nb, nq, ns), μB, μQ, μS as functions of (e, nb, nq, ns)

- **Grid Settings**:

  - ``NT, NB, NQ, NS``: Number of points in each direction
  - ``Ttilde, muBtilde, muQtilde, muStilde``: Ranges as ``[MIN, MAX, N]``
  - ``AutoSetBoundaries``: If true, automatically sets variable boundaries based on the input table

- **Numerical Parameters**:

  - ``Accuracy``: Convergence threshold for root finding
  - ``MAXITER``: Maximum iterations allowed
  - ``Number_of_cores``: Used for parallel execution on a cluster

- **Paths**:

  - ``EoS_table``: Name of the input EoS table
  - ``OutputFolder``: Folder to store the results

-----------------------
Running the Code
-----------------------

Once you've configured ``parameters.py``, you can run the main inversion script:

.. code-block:: bash

   python3 EoSInverter.py parameters.py

-----------------------
Post-Processing
-----------------------

After the inversion finishes, run the post-processing step to analyze the results:

.. code-block:: bash

   python3 PostProcessor.py [OutputFolder] parameters.py

Replace ``[OutputFolder]`` with the actual folder name specified in your configuration.

-----------------------
Project Structure
-----------------------

.. code-block:: text

   EoSInverter/
   ├── EoSInverter.py        # Main inversion script
   ├── PostProcessor.py      # Post-processing utility
   ├── parameters.py         # Configuration file
   ├── requirements.txt      # Python dependencies
   └── ...

-----------------------
Related Projects
-----------------------

This repository is a minimal version of the full framework:
`EoS-TrENsMUTher <https://github.com/gpihan/EoS-TrENsMUTher>`_ (coming soon)


