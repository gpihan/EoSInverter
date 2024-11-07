#/usr/bin/env python3
# Copyright Gregoire Pihan @ 2024

import numpy as np
import sys
from os import path
from scipy import interpolate
import array
from itertools import product
import time
import multiprocessing
import os
import pickle
import random
import struct
from concurrent.futures import ProcessPoolExecutor, as_completed
from utils import read_parameters

def key(iT):
    return  iT

def ToTilde(E): # Input GeV^4, GeV^3
    Ttilde = (12/(19*np.pi**2) * E)**0.25# GeV
    return Ttilde, muBtilde

def ToEN(Tt):
    e = 19*np.pi**2/12 * Tt**4 
    return np.array([e])

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
    dTtilde   = TILDE_TABLES["T"]["d"] 
    Btm = TILDE_TABLES["MUB"]["MIN"] 
    Ttm = TILDE_TABLES["T"]["MIN"] 
    NBtilde = TILDE_TABLES["MUB"]["n"] 
    NTtilde = TILDE_TABLES["T"]["n"] 
    return [Ttm, dTtilde, NTtilde]

def extract_DATA_from(Arr):
    L = list(set(Arr))
    d = {}
    d["n"] = len(L)
    d["MIN"] = np.min(L)
    d["MAX"] = np.max(L)
    d["Table"] = Arr
    d["Arr"] = np.linspace(d["MIN"], d["MAX"], d["n"])
    return d 

def check_lim(X, XMIN, XMAX, tol=1e-8):
    Inside = [(x>=xmin and x<=xmax) for x,xmin,xmax in zip(X, XMIN, XMAX)]
    OnLow = [np.abs(x-xmin) < tol for x,xmin in zip(X,XMIN)]
    OnHigh = [np.abs(x-xmax) < tol for x,xmax in zip(X,XMAX)]
    return all([(I or L or H) for I,L,H in zip(Inside, OnLow, OnHigh)])

def get_EN_interpolation_at(TMUS, INTERP):
    return np.array([f(np.array(TMUS))[0] for f in INTERP.values()])

def getJacobian(local_TMUS, INTERP, absol_deriv=1e-3, deriv_percent=0.02):
    Jacobian = []
    dts = np.array([max(absol_deriv, abs(local_quantity*deriv_percent)) for local_quantity in local_TMUS])
    L_TMUS = list(local_TMUS)
    for f in INTERP.values():
        Jacobian_row = []
        for i, (lv, hv, dx) in enumerate(zip(local_TMUS - dts, local_TMUS + dts, dts)):
            LV, HV = np.array(L_TMUS[:i] + [lv] + L_TMUS[i+1:]), np.array(L_TMUS[:i] + [hv] + L_TMUS[i+1:])
            Jacobian_row.append((f(HV)[0]-f(LV)[0])/(2*dx))
        Jacobian.append(Jacobian_row)
    return np.array(Jacobian)

def Newton(EN, guessSol, INTERP, MTP):
    iterations, status = 0, True
    TMUs_sol = guessSol
    DeltaEN = get_EN_interpolation_at(TMUs_sol, INTERP) - EN 
    Newton_criterion_not_meet = any([(np.abs(Delta)>MTP["ACCURACY"]) for Delta in DeltaEN])
    while Newton_criterion_not_meet:
        Inverse_jacobian = np.linalg.inv(getJacobian(TMUs_sol, INTERP))
        TMUs_sol -= np.dot(Inverse_jacobian, DeltaEN)
        DeltaEN = get_EN_interpolation_at(TMUs_sol, INTERP) - EN 
        Newton_criterion_not_meet = any([(np.abs(Delta)>MTP["ACCURACY"]) for Delta in DeltaEN])
        iterations += 1
        if iterations > MTP["MAXITER"]:
            status = False
            break
    status = check_lim(TMUs_sol, MTP["MINS"], MTP["MAXS"]) and status
    return status, TMUs_sol[0], TMUs_sol[1]

def binary_search_1d(ed_local, MUs, INTERP, TMU_TABLES, MTP):
    iteration, LMUs = 0, list(MUs)
    T_min = TMU_TABLES["T"]["MIN"]; T_max = TMU_TABLES["T"]["MAX"]
    (e_low, e_up) = (INTERP["e"](np.array([T_]+ LMUs))[0] for T_ in [T_min, T_max])
    if (ed_local < e_low):
        return(T_min)
    elif (ed_local > e_up):
        return(T_max)
    else:
        while iteration < MTP["MAXITER"]:
            T_mid = (T_max + T_min)/2.
            e_mid = INTERP["e"](np.array([T_mid] + LMUs))[0]
            abs_err = abs(e_mid - ed_local)
            rel_err = abs_err/abs(e_mid + ed_local + 1e-15)
            if rel_err <= MTP["ACCURACY"] and abs_err <= MTP["ACCURACY"] * 1e-2:
                return T_mid
            if (ed_local < e_mid):
                T_max = T_mid
            else:
                T_min = T_mid
            iteration += 1
        return T_mid

def invert_EOS_tables(T_ID, Ttilde, TILDE_TABLES, TMU_TABLES, INTERP, MTP, guessSol=[0.05, 0.0]):
    hyper_index = key(T_ID)
    EN = ToEN(Ttilde)
    success = False
    try:
        success, T_local = Newton(EN, guessSol, INTERP, MTP)
    except:
        T_local = binary_search_1d(EN, INTERP, TMU_TABLES, MTP)
        success = True
    if not success:
        T_local = binary_search_1d(EN, INTERP, TMU_TABLES, MTP)

    (P_local, s_local) = (f([T_local])[0] for f in MTP["INTERP_PS"].values())
    return [T_local, P_local, s_local, int(hyper_index)]


if __name__ == "__main__":

    try:
        ParamPath = str(sys.argv[1])
        Ttilde = float(sys.argv[2])
        T_ID = float(sys.argv[3])
        OUTPUT_FOLDER = str(sys.argv[4])
    except:
        print("Usage: python3 {} EOS_table_file Ttilde ID output folder".format(sys.argv[0]))
        exit(1)
    
    Param = read_parameters(ParamPath)
    EoS_table_file = Param["EoS_table"] 
    ACCURACY = Param["Accuracy"] 
    MAXITER  = Param["MAXITER"] 
    N_CORES = Param["Number_of_cores"] 
    RunMode = Param["RunMode"]

    hbarC = 0.19733

    eos_table = np.loadtxt(EoS_table_file)
    Thermodynamic_quantities = ["T"]
    Dynamic_quantities = ["e"]
    Press_Entro = ["P", "S"]
    TMU_TABLES = {Quantity_name:extract_DATA_from(eos_table[:, i]) for i, Quantity_name  in enumerate(Thermodynamic_quantities)}
    MINS = [TMU_TABLES[Quantity]["MIN"] for Quantity in Thermodynamic_quantities]
    MAXS = [TMU_TABLES[Quantity]["MAX"] for Quantity in Thermodynamic_quantities]

    # EN in GeV**powers
    EN_powers = [4.0]
    # reshape Temperature
    
    NT, NB, NQ, NS = (TMU_TABLES[Q]["n"] for Q in Thermodynamic_quantities)
    reshaped_T = TMU_TABLES["T"]["Table"].reshape(NT, NB, NQ, NS)
    EN_TABLES = {Quantity_name:eos_table[:, i].reshape(NT, NB, NQ, NS)*reshaped_T**j for i, j, Quantity_name in zip(range(4,7), EN_powers, Dynamic_quantities)} 

    PS_powers = [4.0]
    PS_TABLES = {Quantity_name:eos_table[:, i].reshape(NT, NB, NQ, NS)*reshaped_T**j/hbarC**3 for i, j, Quantity_name in zip(range(7,9), PS_powers, Press_Entro)} 

    GRID = tuple(TMU_TABLES[Q]["Arr"] for Q in Thermodynamic_quantities)
    INTERP = {Quantity_name:interpolate.RegularGridInterpolator(GRID, EN_TABLES[Quantity_name], bounds_error=False, fill_value=None) for Quantity_name in Dynamic_quantities}          
    INTERP_PS = {Quantity_name:interpolate.RegularGridInterpolator(GRID, PS_TABLES[Quantity_name], bounds_error=False, fill_value=None) for Quantity_name in Press_Entro}

    META_PARAMS = {"ACCURACY":ACCURACY, "MAXITER":MAXITER, "MINS":MINS, "MAXS":MAXS, "INTERP_PS":INTERP_PS}

    TILDE_BOUNDARIES = [Param["Ttilde"], Param["muBtilde"], Param["muQtilde"], Param["muStilde"]]

    TILDE_TABLES = {Quantity:fill_TILDE(DATA) for Quantity, DATA in zip(Thermodynamic_quantities, TILDE_BOUNDARIES)}

    lists = [list(range(TILDE_TABLES["MUQ"]["n"])), 
             list(range(TILDE_TABLES["MUB"]["n"]))]

    BQS_combinations = list(product(*lists))
    results = []
    with ProcessPoolExecutor(max_workers=N_CORES) as executor:
        future_to_x = {executor.submit(invert_EOS_tables, T_ID, Ttilde, iQt, iBt, TILDE_TABLES, TMU_TABLES, INTERP, META_PARAMS): (iQt, iBt) 
                       for iQt, iBt in BQS_combinations}
        for future in as_completed(future_to_x):  
            x = future_to_x[future]
            results = future.result()
            with open(OUTPUT_FOLDER+"/TEMP_unordered_inversion_"+str(int(T_ID))+".dat", "a") as file:
                #if RunMode == 0:
                #    print("Writting ", results)
                file.write(" ".join(map(str, results)) + "\n")
        file.close()
