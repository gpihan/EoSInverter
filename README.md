# EoSInverter


This repository is a simplified version of the code 
that will be public at https://github.com/gpihan/EoS-TrENsMUTher


# Step 1: Download the code

Open a terminal and write: 
git clone https://github.com/gpihan/EoSInverter.git

# Step 2: Change the relevant parameters and run the code


The parameters are in "parameters.py" file in the form of a python dictionary. 

"RunMode"   # Run Mode: 0 -> Run locally on one CPU\
            # Run Mode: 1 -> Run on slurm cluster

"Dimension":4, # Dimensions = 1 gives T(e)\
               # Dimensions = 2 gives T(e, nb), mub(T, nb)\
               # Dimensions = 3 gives T(e, nb, nq), mub(e, nb, nq), muq(e, nb, nq)\
               # Dimensions = 4 gives T(e, nb, nq, ns), mub(e, nb, nq, ns), muq(e, nb, nq, ns), mus(e, nb, nq, ns)

"AutoSetBoundaries" # Sets the boundary of the tilde variables automatically considering the input table boundaries.\
"NT" # number of points in the Temperature direction\
"NB" # number of points in the baryon direction\
"NQ" # number of points in the electric charge direction\
"NS" # number of points in the strangeness direction
 
"Ttilde" # [MIN, MAX, N] # set the range of tilde temperature\
"muBtilde" # [MIN, MAX, N] # set the range of mubtilde temperature\
"muQtilde" # [MIN, MAX, N] # set the range of muqtilde temperature\
"muStilde" # [MIN, MAX, N] # set the range of mustilde temperature

"Accuracy" # accuracy of the root finding algorithm\
"MAXITER"  # maximum number of iterations.\
"Number_of_cores" # Number of cores for cluster parallellized calculations.\
"EoS_table" # The name of the EoS table giving the hydro variables as a function of the thermo variables.\
"OutputFolder", the output folder path.

To run the code: 

python3 EoSInverter.py parameters.py

# Step 3: Proceed to post process

In the EoSInverter folder write \
python3 PostProcessor.py [OutputFolder] parameters.py

