import sys
import os
import importlib
from utils import *
import numpy as np

try:
    param_path = sys.argv[1]
except FileNotFoundError:
    print("Parameter file is not found or specified.")
    sys.exit()

Param = read_parameters(param_path)
folder_name = create_folder(Param["OutputFolder"])

dimension = Param["dimension"]

if Param["RunMode"] == 0:
    # Run locally
    TArr = np.linspace(Param["Ttilde"][0], Param["Ttilde"][1], Param["Ttilde"][2])
    for i, T in enumerate(TArr):
        print("Running temperature", T)
        os.system("python3 mapEOS_"+str(dimension)+"D.py "+param_path+" "+ str(T) +" "+ str(i) +" "+ folder_name)

elif Param["RunMode"] == 1:
    # Run on Slurm cluster
    queue = "primary"
    Jobname = "InvertEoS_"+folder_name
    Script_name = "JobArray_"+folder_name+".sh"

    dTtilde = (Param["Ttilde"][1] - Param["Ttilde"][0])/(Param["Ttilde"][2]-1)
    create_job_script(Jobname, queue, Param["Number_of_cores"], Param["Ttilde"][2]-1, 
                      Param["Ttilde"][0], Param["Ttilde"][1], dTtilde, param_path, folder_name, output_filename=Script_name)

    os.system("sbatch "+Script_name)


