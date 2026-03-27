import os
import pickle
import sys

import numpy as np
from utils import (
    SetBoundaries,
    SetENBoundaries,
    create_folder,
    create_job_script,
    read_parameters,
)

try:
    param_path = sys.argv[1]
except FileNotFoundError:
    print("Parameter file is not found or specified.")
    sys.exit()

Param = read_parameters(param_path)
folder_name = create_folder(Param["OutputFolder"])

dimension = Param["Dimension"]
tilde_switch = Param.get("TildeSwitch", True)

if Param["AutoSetBoundaries"]:
    Bound = SetBoundaries(Param["EoS_table"], Param["Dimension"])
    Names = ["Ttilde", "muBtilde", "muQtilde", "muStilde"]
    Npoints = ["NT", "NB", "NQ", "NS"]
    OutputFolder = Param["OutputFolder"]
    Boundaries = {}
    for i, name, Npoint in zip(range(len(Bound)), Names, Npoints):
        Boundaries[name] = list(Bound[i]) + [Param[Npoint]]

    # Optional EN-space boundaries used when TildeSwitch is off.
    if not tilde_switch and dimension == 2:
        BoundEN = SetENBoundaries(Param["EoS_table"], Param["Dimension"])
        Boundaries["e"] = [BoundEN[0][0], BoundEN[0][1], Param["NT"]]
        Boundaries["nB"] = [BoundEN[1][0], BoundEN[1][1], Param["NB"]]

    with open(OutputFolder + "_boundaries_temp.dat", "wb") as f:
        pickle.dump(Boundaries, f)
    if tilde_switch:
        TTILDE_MIN, TTILDE_MAX, TTILDE_N = (
            Boundaries["Ttilde"][0],
            Boundaries["Ttilde"][1],
            Boundaries["Ttilde"][2],
        )
    elif dimension == 2:
        TTILDE_MIN, TTILDE_MAX, TTILDE_N = (
            Boundaries["e"][0],
            Boundaries["e"][1],
            Boundaries["e"][2],
        )
    else:
        raise ValueError(
            "TildeSwitch=False is currently supported only for Dimension=2."
        )
else:
    if tilde_switch:
        TTILDE_MIN, TTILDE_MAX, TTILDE_N = (
            Param["Ttilde"][0],
            Param["Ttilde"][1],
            Param["Ttilde"][2],
        )
    elif dimension == 2:
        TTILDE_MIN, TTILDE_MAX, TTILDE_N = (
            Param["e"][0],
            Param["e"][1],
            Param["e"][2],
        )
    else:
        raise ValueError(
            "TildeSwitch=False is currently supported only for Dimension=2."
        )


if Param["RunMode"] == 0:
    # Run locally
    TArr = np.linspace(TTILDE_MIN, TTILDE_MAX, TTILDE_N)
    for i, T in enumerate(TArr):
        print("Running temperature", T)
        os.system(
            "python3 mapEOS_"
            + str(dimension)
            + "D.py "
            + param_path
            + " "
            + str(T)
            + " "
            + str(i)
            + " "
            + folder_name
        )
elif Param["RunMode"] == 1:
    # Run on Slurm cluster
    queue = "primary"
    Jobname = "InvertEoS_" + folder_name
    Script_name = "JobArray_" + folder_name + ".sh"

    dTtilde = (TTILDE_MAX - TTILDE_MIN) / (TTILDE_N - 1)
    create_job_script(
        Jobname,
        queue,
        Param["Number_of_cores"],
        TTILDE_N - 1,
        TTILDE_MIN,
        TTILDE_MAX,
        dTtilde,
        param_path,
        folder_name,
        output_filename=Script_name,
    )

    os.system("sbatch " + Script_name)
elif Param["RunMode"] == 2:
    # Generate a bash script to run locally later
    TArr = np.linspace(TTILDE_MIN, TTILDE_MAX, TTILDE_N)
    # Create a directory to store individual job scripts
    scripts_dir = f"jobs_{dimension}D_{folder_name}"
    try:
        os.makedirs(scripts_dir, exist_ok=True)
    except Exception:
        pass

    # Create a master script that runs all jobs sequentially
    master_script = os.path.join(scripts_dir, f"run_all_{dimension}D_{folder_name}.sh")
    master_lines = [
        "#!/usr/bin/env bash\n",
        "set -euo pipefail\n\n",
        'cd "$(dirname "$0")"\n\n',
    ]

    submit_failures = []

    for i, T in enumerate(TArr):
        job_script = os.path.join(
            scripts_dir, f"job_{dimension}D_{folder_name}_Tidx{i}.sh"
        )
        job_lines = [
            "#!/bin/bash\n",
            "#PBS -l walltime=28:00:00\n",
            "#PBS -l select=2:mem=10gb\n",
            'cd "/storage/brno12-cerit/home/poledto1/EoSInverter"\n',
            "module add python\n",
            "source venv/bin/activate\n",
            "python3 mapEOS_"
            + str(dimension)
            + "D.py "
            + param_path
            + " "
            + str(T)
            + " "
            + str(i)
            + " "
            + folder_name
            + "\n",
        ]
        with open(job_script, "w") as jf:
            jf.writelines(job_lines)
        try:
            os.chmod(job_script, 0o755)
        except Exception:
            pass
        rc = os.system("qsub " + job_script)
        if rc != 0:
            submit_failures.append((i, job_script, rc))

        # Add to master script
        master_lines.append(f"./{os.path.basename(job_script)}\n")

    with open(master_script, "w") as mf:
        mf.writelines(master_lines)
    try:
        os.chmod(master_script, 0o755)
    except Exception:
        pass

    if submit_failures:
        failure_log = os.path.join(scripts_dir, "submit_failures.txt")
        with open(failure_log, "w") as lf:
            for i, job_script, rc in submit_failures:
                lf.write(f"Tidx{i}\t{job_script}\treturn_code={rc}\n")
        print(
            f"Warning: {len(submit_failures)} job(s) failed to submit. "
            f"Details in {failure_log}."
        )
