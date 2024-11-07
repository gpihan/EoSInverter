# EoSInverter


This repository is a simplified version of the code 
that will be public at https://github.com/gpihan/EoS-TrENsMUTher


# Step 1: Download the code

Open a terminal and write: \

git clone https://github.com/gpihan/EoSInverter.git

# Step 2: Change the relevant parameters

Essentially, modify the dimension under "Dimension"\
The location of the EoS table under "EoSTable"\
The output folder where you want the inverted EoS table to\
be output.

# Step 3: Proceed to post process

In the EoSInverter folder write \
python3 PostProcessor.py [OutputFolder] parameters.py

# Step 4: Analyze the inverted table

You get several files all shaped in the\
same way:\

a header of size 3*Dimension (index 3*Dimension - 1)\
and long line that should be read using the\
"key" function at the begining of each map_...py scripts.\
(We can see together)
