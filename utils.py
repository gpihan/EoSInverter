import os
import sys
import importlib
import random
import string

def read_parameters(path):
        try:
           spec = importlib.util.spec_from_file_location("general_parameters", path)
           module = importlib.util.module_from_spec(spec)
           sys.modules["general_parameters"] = module
           spec.loader.exec_module(module)
           return module.general_parameters
        except FileNotFoundError:
            print(f"File '{file_path}' not found.")
            sys.exit()
        except Exception as e:
            print(f"An error occurred: {e}")
            sys.exit()

def create_folder(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        return folder_path
    else:
        random_char = random.choice(string.ascii_letters)
        folder_path_ = folder_path+"_"+str(random_char)
        return create_folder(folder_path_)

def create_job_script(job_name, queue, n_cores, nx, xm, xM, dx, ParametersPath, folder_name, output_filename='JobArray.sh') -> None:
    script_template = """#!/bin/bash

#SBATCH --job-name="{JOBNAMES}"
#SBATCH -q {QUEUE}
#SBATCH -N 1
#SBATCH -n {NCORES}
#SBATCH --mem=8G
#SBATCH --constraint=intel
#SBATCH -o output_%A_%a.out
#SBATCH -e errors_%A_%a.err
#SBATCH -t 20-0:0:0
#SBATCH --array=0-{NX}

TArr=($(seq {Xm} {dX} {XM}))
index=$SLURM_ARRAY_TASK_ID
T_value=${{TArr[index]}}
python3 mapEOS_4D.py {Parameters_path} $T_value $index {folder_name}  
"""

    script_content = script_template.format(
        JOBNAMES=job_name,
        QUEUE=queue,
        NCORES=n_cores,
        NX=nx,
        Xm=xm,
        XM=xM,
        dX=dx,
        Parameters_path=ParametersPath,
        folder_name=folder_name
    )
    with open(output_filename, 'w') as f:
        f.write(script_content)

    print(f"Job script written to {output_filename}")


