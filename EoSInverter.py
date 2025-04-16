import sys
import os
import importlib
from utils import *
import numpy as np
import pickle

try:
    param_path = sys.argv[1]
except FileNotFoundError:
    print("Parameter file is not found or specified.")
    sys.exit()

Param = read_parameters(param_path)
folder_name = create_folder(Param["OutputFolder"])

dimension = Param["Dimension"]

if Param["AutoSetBoundaries"]:
    Bound = SetBoundaries(Param["EoS_table"], Param["Dimension"])
    Names = ["Ttilde", "muBtilde", "muQtilde", "muStilde"]
    Npoints = ["NT", "NB", "NQ", "NS"]
    Boundaries = {}
    for i, name, Npoint in zip(range(len(Bound)), Names, Npoints):
        Boundaries[name] = list(Bound[i]) + [Param[Npoint]]
    with open('boundaries_temp.dat', 'wb') as f:
        pickle.dump(Boundaries, f)
    TTILDE_MIN, TTILDE_MAX, TTILDE_N = Boundaries["Ttilde"][0], Boundaries["Ttilde"][1], Boundaries["Ttilde"][2]
else:
     TTILDE_MIN, TTILDE_MAX, TTILDE_N = Param["Ttilde"][0], Param["Ttilde"][1], Param["Ttilde"][2]


if Param["RunMode"] == 0:
    # Run locally
    TArr = np.linspace(TTILDE_MIN, TTILDE_MAX, TTILDE_N)
    for i, T in enumerate(TArr):
        print("Running temperature", T)
        os.system("python3 mapEOS_"+str(dimension)+"D.py "+param_path+" "+ str(T) +" "+ str(i) +" "+ folder_name)

elif Param["RunMode"] == 1:
    # Run on Slurm cluster
    queue = "primary"
    Jobname = "InvertEoS_"+folder_name
    Script_name = "JobArray_"+folder_name+".sh"

    dTtilde = (TTILDE_MAX - TTILDE_MIN)/(TTILDE_N-1)
    create_job_script(Jobname, queue, Param["Number_of_cores"], TTILDE_N-1, 
                      TTILDE_MIN, TTILDE_MAX, dTtilde, param_path, folder_name, output_filename=Script_name)

    os.system("sbatch "+Script_name)


