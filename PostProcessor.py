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


def write_header(TILDE_TABLES, dimension):
    keys = ["T", "MUB", "MUQ", "MUS"]
    header = []
    for key in keys:
        if key in TILDE_TABLES:
            d = TILDE_TABLES[key]
            header.extend([d["MIN"], d["d"], d["n"]])
        else:
            header.extend([0.0, 0.0, 1])
    return header


try:
    output_folder = sys.argv[1]
    param_path = sys.argv[2]
except IndexError:
    print("Usage: PostProcessor.py output_folder param_path")
    sys.exit(1)

Param = read_parameters(param_path)
dimension = Param.get("Dimension", 1)

dim_map = {
    1: ["T"],
    2: ["T", "MUB"],
    3: ["T", "MUB", "MUQ"],
    4: ["T", "MUB", "MUQ", "MUS"],
}
Thermodynamic_quantities = dim_map.get(dimension, ["T", "MUB"])


if Param.get("AutoSetBoundaries", False):
    with open("boundaries_temp.dat", "rb") as f:
        BOUNDS = pickle.load(f)
else:
    BOUNDS = Param

quantity_keys = {"T": "Ttilde", "MUB": "muBtilde", "MUQ": "muQtilde", "MUS": "muStilde"}
TILDE_BOUNDARIES = [BOUNDS[quantity_keys[q]] for q in Thermodynamic_quantities]
TILDE_TABLES = {
    q: fill_TILDE(data) for q, data in zip(Thermodynamic_quantities, TILDE_BOUNDARIES)
}

header = write_header(TILDE_TABLES, dimension)
NT = TILDE_TABLES["T"]["n"] if "T" in TILDE_TABLES else 1
FNAMES = ["t", "mub", "p", "s"]

for i in range(NT):
    print(f"Treating case: {i}")
    file_path = os.path.join(output_folder, f"TEMP_unordered_inversion_{i}.dat")
    DATA = np.loadtxt(file_path)
    DATA = DATA[np.argsort(DATA[:, -1])]

    mode = "wb" if i == 0 else "ab"
    for j, fname in enumerate(FNAMES):
        with open(f"EoS_{fname}_b.dat", mode) as f:
            if i == 0:
                for h in header:
                    f.write(struct.pack("f", h))
            for val in DATA[:, j]:
                f.write(struct.pack("f", val))


for fn in FNAMES:
    shutil.move(f"EoS_{fn}_b.dat", output_folder)


raw_data_dir = os.path.join(output_folder, "RAW_DATA")
os.makedirs(raw_data_dir, exist_ok=True)

for file in os.listdir(output_folder):
    if file.startswith("TEMP_unordered_inversion_"):
        shutil.move(os.path.join(output_folder, file), raw_data_dir)

if os.path.exists("boundaries_temp.dat"):
    os.remove("boundaries_temp.dat")

for pattern in ["JobArray*", "output*", "errors*"]:
    os.system(f"rm -rf {pattern}")
