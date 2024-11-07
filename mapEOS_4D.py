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

def key(iT, ib, iq, i_s, Nb=5, Nq=5, Ns=5):
    return  ((iT*Ns + i_s)*Nq + iq)*Nb + ib

def ToTilde(E, Nb, Nq, Ns): # Input GeV^4, GeV^3
    Ttilde = (12/(19*np.pi**2) * E)**0.25# GeV
    muBtilde = (5 * Nb - Nq + 2 * Ns)/ Ttilde**2 # GeV
    muQtilde = (-1 * Nb + 2 * Nq - Ns)/ Ttilde**2 # GeV
    muStilde = (2 * Nb - Nq + 2 * Ns)/ Ttilde**2 # GeV
    return Ttilde, muBtilde, muQtilde, muStilde

def ToEN(Tt, mbt, mqt, mst):
    e = 19*np.pi**2/12 * Tt**4 
    nb = 1/3*(mbt - mst)*Tt**2
    nq = 1/3*(2*mqt + mst)*Tt**2
    ns = 1/3*(mqt - mbt + 3 * mst)*Tt**2
    return np.array([e, nb, nq, ns])

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

def Newton4D(EN, guessSol, INTERP, MTP):
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
    return status, TMUs_sol[0], TMUs_sol[1], TMUs_sol[2], TMUs_sol[3],

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

def binary_search_2d(ed_local, nB_local, MUs, INTERP, TMU_TABLES, MTP):
    iteration, LMUs = 0, list(MUs)
    muB_min = TMU_TABLES["MUB"]["MIN"]; muB_max = TMU_TABLES["MUB"]["MAX"]
    (T_min, T_max) = (binary_search_1d(ed_local, [muB_]+LMUs, INTERP, TMU_TABLES, MTP) for muB_ in [muB_min, muB_max])
    (nB_min, nB_max) = (INTERP["nB"](np.array([T_, muB_]+ LMUs))[0] for T_, muB_ in zip([T_min, T_max], [muB_min, muB_max]))
    if (nB_local < nB_min):
        return T_min, muB_min
    elif (nB_local > nB_max):
        return T_max, muB_max
    else:
        while iteration < MTP["MAXITER"]:
            muB_mid = (muB_min + muB_max)/2.
            T_mid = binary_search_1d(ed_local, [muB_mid] + LMUs, INTERP, TMU_TABLES, MTP)
            nB_mid = INTERP["nB"](np.array([T_mid, muB_mid] + LMUs))[0]

            abs_err = abs(nB_mid - nB_local)
            rel_err = abs_err/max(1e-15, abs(nB_mid + nB_local))

            if rel_err <= MTP["ACCURACY"] and abs_err <= MTP["ACCURACY"] * 1e-2:
                return T_mid, muB_mid

            if (nB_local < nB_mid):
                muB_max = muB_mid
            else:
                muB_min = muB_mid
            iteration += 1
        return T_mid, muB_mid

def binary_search_3d(ed_local, nB_local, nQ_local, MUs, INTERP, TMU_TABLES, MTP):
    iteration, LMUs = 0, [MUs]
    muQ_min = TMU_TABLES["MUQ"]["MIN"]; muQ_max = TMU_TABLES["MUQ"]["MAX"]
    ((T_min, muB_min), (T_max, muB_max)) = (binary_search_2d(ed_local, nB_local, [muQ_] + LMUs, INTERP, TMU_TABLES, MTP) for muQ_ in [muQ_min, muQ_max])
    (nQ_min, nQ_max) = (INTERP["nQ"](np.array([T_, muB_, muQ_]+ LMUs))[0] for T_, muB_, muQ_ in zip([T_min, T_max], [muB_min, muB_max], [muQ_min, muQ_max]))
    if (nQ_local < nQ_min):
        return(T_min, muB_min, muQ_min)
    elif (nQ_local > nQ_max):
        return(T_max, muB_max, muQ_max)
    else:
        while iteration < MTP["MAXITER"]:
            muQ_mid = (muQ_min + muQ_max)/2.
            T_mid, muB_mid = binary_search_2d(ed_local, nB_local, [muQ_mid] + LMUs, INTERP, TMU_TABLES, MTP)
            nQ_mid = INTERP["nQ"](np.array([T_mid, muB_mid, muQ_mid] + LMUs))[0]
            abs_err = abs(nQ_mid - nQ_local)
            rel_err = abs_err/max(1e-15, abs(nQ_mid + nQ_local))

            if rel_err <= MTP["ACCURACY"] and abs_err <= MTP["ACCURACY"] * 1e-2:
                return T_mid, muB_mid, muQ_mid

            if (nQ_local < nQ_mid):
                muQ_max = muQ_mid
            else:
                muQ_min = muQ_mid
            iteration += 1
        return T_mid, muB_mid, muQ_mid

def binary_search_4d(EN, INTERP, TMU_TABLES, MTP):
    iteration = 0
    ed_local, nB_local, nQ_local, nS_local = EN[0], EN[1], EN[2], EN[3]
    muS_min = TMU_TABLES["MUS"]["MIN"]; muS_max = TMU_TABLES["MUS"]["MAX"]

    ((T_min, muB_min, muQ_min), (T_max, muB_max, muQ_max)) = (binary_search_3d(ed_local, nB_local, nQ_local, muS_, INTERP, TMU_TABLES, MTP) for muS_ in [muS_min, muS_max])
    (nS_min, nS_max) = (INTERP["nS"](np.array([T_, muB_, muQ_, muS_]))[0] 
                        for T_, muB_, muQ_, muS_ in zip([T_min, T_max], [muB_min, muB_max], [muQ_min, muQ_max], [muS_min, muS_max]))
    if (nS_local < nS_min):
        return(T_min, muB_min, muQ_min, muS_min)
    elif (nS_local > nS_max):
        return(T_max, muB_max, muQ_max, muS_max)
    else:
        while iteration < MTP["MAXITER"]:
            muS_mid = (muS_min + muS_max)/2.
            T_mid, muB_mid, muQ_mid = binary_search_3d(ed_local, nB_local, nQ_local, muS_mid, INTERP, TMU_TABLES, MTP)
            nS_mid = INTERP["nS"](np.array([T_mid, muB_mid, muQ_mid, muS_mid]))[0]
            abs_err = abs(nS_mid - nS_local)
            rel_err = abs_err/max(1e-15, abs(nS_mid + nS_local))

            if rel_err <= MTP["ACCURACY"] and abs_err <= MTP["ACCURACY"] * 1e-2:
                return T_mid, muB_mid, muQ_mid, muS_mid

            if (nS_local < nS_mid):
                muS_max = muS_mid
            else:
                muS_min = muS_mid
            iteration += 1
        return T_mid, muB_mid, muQ_mid, muS_mid

def invert_EOS_tables(T_ID, Ttilde, iSt, iQt, iBt, TILDE_TABLES, TMU_TABLES, INTERP, MTP, guessSol=[0.05, 0.0, 0.0, 0.0]):
    mbtilde, mqtilde, mstilde = TILDE_TABLES["MUB"]["Arr"][iBt], TILDE_TABLES["MUQ"]["Arr"][iQt], TILDE_TABLES["MUS"]["Arr"][iSt]
    hyper_index = key(T_ID, iBt, iQt, iSt, Nb=TILDE_TABLES["MUB"]["n"], Nq=TILDE_TABLES["MUQ"]["n"], Ns=TILDE_TABLES["MUS"]["n"])
    EN = ToEN(Ttilde, mbtilde, mqtilde, mstilde)
    success = False
    try:
        success, T_local, muB_local, muQ_local, muS_local = Newton4D(EN, guessSol, INTERP, MTP)
    except:
        T_local, muB_local, muQ_local, muS_local = binary_search_4d(EN, INTERP, TMU_TABLES, MTP)
        success = True
    if not success:
        T_local, muB_local, muQ_local, muS_local = binary_search_4d(EN, INTERP, TMU_TABLES, MTP)

    (P_local, s_local) = (f([T_local, muB_local, muQ_local, muS_local])[0] for f in MTP["INTERP_PS"].values())
    return [T_local, muB_local, muQ_local, muS_local, P_local, s_local, int(hyper_index)]


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
    Thermodynamic_quantities = ["T", "MUB", "MUQ", "MUS"]
    Dynamic_quantities = ["e", "nB", "nQ", "nS"]
    Press_Entro = ["P", "S"]
    TMU_TABLES = {Quantity_name:extract_DATA_from(eos_table[:, i]) for i, Quantity_name  in enumerate(Thermodynamic_quantities)}
    MINS = [TMU_TABLES[Quantity]["MIN"] for Quantity in Thermodynamic_quantities]
    MAXS = [TMU_TABLES[Quantity]["MAX"] for Quantity in Thermodynamic_quantities]

    # EN in GeV**powers
    EN_powers = [4.0, 3.0, 3.0, 3.0]
    # reshape Temperature
    
    NT, NB, NQ, NS = (TMU_TABLES[Q]["n"] for Q in Thermodynamic_quantities)
    reshaped_T = TMU_TABLES["T"]["Table"].reshape(NT, NB, NQ, NS)
    EN_TABLES = {Quantity_name:eos_table[:, i].reshape(NT, NB, NQ, NS)*reshaped_T**j for i, j, Quantity_name in zip(range(4,8), EN_powers, Dynamic_quantities)} 

    PS_powers = [4.0, 3.0]
    PS_TABLES = {Quantity_name:eos_table[:, i].reshape(NT, NB, NQ, NS)*reshaped_T**j/hbarC**3 for i, j, Quantity_name in zip(range(8,10), PS_powers, Press_Entro)} 

    GRID = tuple(TMU_TABLES[Q]["Arr"] for Q in Thermodynamic_quantities)
    INTERP = {Quantity_name:interpolate.RegularGridInterpolator(GRID, EN_TABLES[Quantity_name], bounds_error=False, fill_value=None) for Quantity_name in Dynamic_quantities}          
    INTERP_PS = {Quantity_name:interpolate.RegularGridInterpolator(GRID, PS_TABLES[Quantity_name], bounds_error=False, fill_value=None) for Quantity_name in Press_Entro}

    META_PARAMS = {"ACCURACY":ACCURACY, "MAXITER":MAXITER, "MINS":MINS, "MAXS":MAXS, "INTERP_PS":INTERP_PS}

    TILDE_BOUNDARIES = [Param["Ttilde"], Param["muBtilde"], Param["muQtilde"], Param["muStilde"]]

    TILDE_TABLES = {Quantity:fill_TILDE(DATA) for Quantity, DATA in zip(Thermodynamic_quantities, TILDE_BOUNDARIES)}

    lists = [list(range(TILDE_TABLES["MUS"]["n"])), 
             list(range(TILDE_TABLES["MUQ"]["n"])), 
             list(range(TILDE_TABLES["MUB"]["n"]))]

    BQS_combinations = list(product(*lists))
    results = []
    with ProcessPoolExecutor(max_workers=N_CORES) as executor:
        future_to_x = {executor.submit(invert_EOS_tables, T_ID, Ttilde, iSt, iQt, iBt, TILDE_TABLES, TMU_TABLES, INTERP, META_PARAMS): (iSt, iQt, iBt) 
                       for iSt, iQt, iBt in BQS_combinations}
        for future in as_completed(future_to_x):  
            x = future_to_x[future]
            results = future.result()
            with open(OUTPUT_FOLDER+"/TEMP_unordered_inversion_"+str(int(T_ID))+".dat", "a") as file:
                #if RunMode == 0:
                #    print("Writting ", results)
                file.write(" ".join(map(str, results)) + "\n")
        file.close()
