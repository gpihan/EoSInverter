import os
import numpy as np
import struct
import sys
import pickle
import shutil
from utils import read_parameters

def fill_TILDE(DATA):
    d = {}
    d["n"] = DATA[2]
    d["MIN"] = DATA[0]
    d["MAX"] = DATA[1]
    d["Arr"] = np.linspace(DATA[0], DATA[1], DATA[2])
    d["d"] = d["Arr"][1] - d["Arr"][0]
    return d


def write_header(tables, quantities):
    header = []
    for q in quantities:
        d = tables[q]
        header.extend([d["MIN"], d["d"], d["n"]])
    return header


try:
    target_folder = sys.argv[1]
    param_path = sys.argv[2]
except IndexError:
    print("Usage: PostProcessor.py output_folder param_path")
    sys.exit(1)

Param = read_parameters(param_path)
dimension = Param["Dimension"]
solver = Param["hydro_model"]

dim_map = {1: ["T"], 2: ["T", "MUB"], 3: ["T", "MUB", "MUQ"], 4: ["T", "MUB", "MUQ", "MUS"]}
Thermo = dim_map.get(dimension, ["T", "MUB"])
quantity_keys = {"T": "Ttilde", "MUB": "muBtilde", "MUQ": "muQtilde", "MUS": "muStilde"}

if Param.get("AutoSetBoundaries", False):
    with open("boundaries_temp.dat", "rb") as f:
        BOUNDS = pickle.load(f)
else:
    BOUNDS = Param


tilde_bounds = [BOUNDS[quantity_keys[q]] for q in Thermo]
TILDE = {q: fill_TILDE(data) for q, data in zip(Thermo, tilde_bounds)}

if solver == "MUSIC":
    fields = ["t", "mub", "p", "s"]
    COL = {"t": 0, "mub": 1, "p": 2, "s": 3}
elif solver == "vHLLE":
    fields = ["e", "nb", "t", "mub", "p", "s"]
    if dimension >= 3:
        fields.append("muq")
    if dimension == 4:
        fields.append("mus")
    COL = {f: i for i, f in enumerate(fields)}
else:
    raise ValueError(f"Unknown solver: {solver}")

NT = TILDE["T"]["n"] if "T" in TILDE else 1

if solver == "MUSIC":
    header = write_header(TILDE, Thermo)
    for i in range(NT):
        print(f"Treating case MUSIC: {i}")
        data = np.loadtxt(os.path.join(target_folder, f"TEMP_unordered_inversion_{i}.dat"))
        data = data[np.argsort(data[:, -1])]
        mode = "wb" if i == 0 else "ab"
        for fname in fields:
            with open(f"EoS_{fname}_b.dat", mode) as f:
                if i == 0:
                    for h in header:
                        f.write(struct.pack('f', h))
                for val in data[:, COL[fname]]:
                    f.write(struct.pack('f', val))
    for f in [f"EoS_{f}_b.dat" for f in fields]:
        shutil.move(f, target_folder)

elif solver == "vHLLE":
    out_file = os.path.join(target_folder, "EoS_all.dat")
    with open(out_file, "w") as f_out:
        for i in range(NT):
            print(f"Treating case vHLLE: {i}")
            file_path = os.path.join(target_folder, f"TEMP_unordered_inversion_{i}.dat")
            data = np.loadtxt(file_path)
            data = data[np.argsort(data[:, -1])]
            
            num_cols = data.shape[1]
            used_fields = fields[:num_cols] 

            for row in data:
                values = [str(row[COL[f]]) for f in used_fields]
                f_out.write(" ".join(values) + "\n")
    print(f"Output saved to {out_file}")
else:
    raise ValueError(f"Unknown solver: {solver}")

raw_dir = os.path.join(target_folder, "RAW_DATA")
os.makedirs(raw_dir, exist_ok=True)
for file in os.listdir(target_folder):
    if file.startswith("TEMP_unordered_inversion_"):
        shutil.move(os.path.join(target_folder, file), raw_dir)
if os.path.exists("boundaries_temp.dat"):
    os.remove("boundaries_temp.dat")
for pat in ["JobArray*", "output*", "errors*"]:
    os.system(f"rm -rf {pat}")
