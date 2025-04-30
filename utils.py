import os
import sys
import importlib
import random
import string
import numpy as np
import pandas as pd

# TArr=($(seq {Xm} {dX} {XM}))


def read_parameters(path):
    try:
        spec = importlib.util.spec_from_file_location("general_parameters", path)
        module = importlib.util.module_from_spec(spec)
        sys.modules["general_parameters"] = module
        spec.loader.exec_module(module)
        return module.general_parameters
    except FileNotFoundError:
        print(f"File '{path}' not found.")
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
        folder_path_ = folder_path + "_" + str(random_char)
        return create_folder(folder_path_)


def create_job_script(
    job_name,
    queue,
    n_cores,
    nx,
    xm,
    xM,
    dx,
    ParametersPath,
    folder_name,
    output_filename="JobArray.sh",
) -> None:
    script_template = r"""#!/bin/bash

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

TArr=($(awk -v min={Xm} -v dx={dX} -v nx={NX} 'BEGIN {{ for(i=0;i<=nx;i++) printf "%.6f ", min + i*dx }}'))
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
        folder_name=folder_name,
    )
    with open(output_filename, "w") as f:
        f.write(script_content)

    print(f"Job script written to {output_filename}")


def Get2DTilde(vec):  # Input GeV^4, GeV^3
    e, nb = vec.T
    Ttilde = (12 / (19 * np.pi**2) * e) ** 0.25  # GeV
    muBtilde = 5 * nb / Ttilde**2  # GeV
    return np.column_stack((Ttilde, muBtilde))


def Get4DTilde(vec):  # Input GeV^4, GeV^3
    e, nb, nq, ns = vec.T
    Ttilde = (12 / (19 * np.pi**2) * e) ** 0.25  # GeV
    muBtilde = (5 * nb - nq + 2 * ns) / Ttilde**2  # GeV
    muQtilde = (-1 * nb + 2 * nq - ns) / Ttilde**2  # GeV
    muStilde = (2 * nb - nq + 2 * ns) / Ttilde**2  # GeV
    return np.column_stack((Ttilde, muBtilde, muQtilde, muStilde))


def GetTilde(vec, dim):
    if dim == 4:
        return Get4DTilde(vec)
    elif dim == 2:
        return Get2DTilde(vec)
    else:
        print("AutoBoundaries are not implemented for 1D and 3D EoS.")
        sys.exit()


def readTable(EoS_path):
    df = pd.read_csv(EoS_path, sep="\s+", comment="#", header=None)
    return df.to_numpy()


def SetBoundaries(EoS_path, dim):
    EOS = readTable(EoS_path)
    T = EOS[:, 0]
    # To do: Check if user wants normalization ------
    exponents = [4] + [3] * (dim - 1)
    print(f"exponents: {exponents}")
    scaling_factors = np.column_stack([T**p for p in exponents])
    print(f"scaling_factors: {scaling_factors}")
    dynVars = EOS[:, dim : 2 * dim] * scaling_factors
    Tildes = GetTilde(dynVars, dim)
    return np.array([[np.min(Tildes[:, i]), np.max(Tildes[:, i])] for i in range(dim)])
