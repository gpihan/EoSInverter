import os
import numpy as np
import struct
import sys
from utils import read_parameters
import pickle

def fill_TILDE(DATA):
    d = {}
    d["n"] = DATA[2] 
    d["MIN"] = DATA[0]
    d["MAX"] = DATA[1]
    d["Arr"] = np.linspace(DATA[0], DATA[1], DATA[2])
    d["d"] = d["Arr"][1] - d["Arr"][0]
    return d

def write_header(TILDE_TABLES):
    dmubtilde = TILDE_TABLES["MUB"]["d"]
    dmuqtilde = TILDE_TABLES["MUQ"]["d"]
    dmustilde = TILDE_TABLES["MUS"]["d"]
    dTtilde   = TILDE_TABLES["T"]["d"] 
    Btm = TILDE_TABLES["MUB"]["MIN"] 
    Qtm = TILDE_TABLES["MUQ"]["MIN"] 
    Stm = TILDE_TABLES["MUS"]["MIN"] 
    Ttm = TILDE_TABLES["T"]["MIN"] 
    NBtilde = TILDE_TABLES["MUB"]["n"] 
    NQtilde = TILDE_TABLES["MUQ"]["n"] 
    NStilde = TILDE_TABLES["MUS"]["n"] 
    NTtilde = TILDE_TABLES["T"]["n"] 
    return [Btm, Qtm, Stm, Ttm, dmubtilde, dmuqtilde, dmustilde, dTtilde, NBtilde, NQtilde, NStilde, NTtilde]

Thermodynamic_quantities = ["T", "MUB", "MUQ", "MUS"]
try:
    output_folder = sys.argv[1] 
    param_path = sys.argv[2]
except:
    print("PostProcessor.py output_folder param_path")
    sys.exit()

Param = read_parameters(param_path)

if Param["AutoSetBoundaries"]:
    with open('boundaries_temp.dat', 'rb') as f:
        Boundaries = pickle.load(f)
    BOUNDS = Boundaries
else:
    BOUNDS = Param

TILDE_BOUNDARIES = [BOUNDS["Ttilde"], BOUNDS["muBtilde"], BOUNDS["muQtilde"], BOUNDS["muStilde"]]
TILDE_TABLES = {Quantity:fill_TILDE(DATA) for Quantity, DATA in zip(Thermodynamic_quantities, TILDE_BOUNDARIES)}

here = output_folder
h_, F, f = next(os.walk(here))
FNAMES = ["t", "mub", "muq", "mus", "p", "s"]
header = write_header(TILDE_TABLES) 

NT = BOUNDS["Ttilde"][2] 
for i in range(NT): 
    print("Treating case : ", i)
    DATA = np.loadtxt(here+"/TEMP_unordered_inversion_"+str(i)+".dat")
    DATA[:] = DATA[np.argsort(DATA[:, -1])]
    if i==0:
        files = {fname:open("EoS_"+fname+"_b.dat", 'wb') for fname in FNAMES}
        for j, f in enumerate(files.values()):
            for h in header:
                f.write(struct.pack('f', h))
            for data in DATA[:,j]:
                f.write(struct.pack('f',data))
            f.close()
    else:
        files = {fname:open("EoS_"+fname+"_b.dat", 'ab') for fname in FNAMES}
        for j, (fn, f) in enumerate(files.items()):
            for data in DATA[:,j]:
                f.write(struct.pack('f',data))
            f.close()

for fn in FNAMES:
    os.system("mv "+"EoS_"+fn+"_b.dat "+ output_folder)
os.system("mkdir "+output_folder+"/RAW_DATA")
os.system("mv "+output_folder+"/TEMP_unordered_inversion_* "+output_folder+"/RAW_DATA")
os.system("rm -rf JobArray*")
os.system("rm -rf output*")
os.system("rm -rf errors*")
os.system("rm -rf boundaries_temp.dat")
