# /usr/bin/env python3
# Copyright Gregoire Pihan @ 2024

import pickle
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy import interpolate
from utils import read_parameters, readTable


def key(iT, ib, Nb=5):
    return iT * Nb + ib


def ToEN(Tt, mbt):
    e = 19 * np.pi**2 / 12 * Tt**4
    nb = 1 / 5 * mbt * Tt**2
    return np.array([e, nb])


def fill_TILDE(DATA):
    d = {}
    d["n"] = DATA[2]
    d["MIN"] = DATA[0]
    d["MAX"] = DATA[1]
    d["Arr"] = np.linspace(DATA[0], DATA[1], DATA[2])
    d["d"] = d["Arr"][1] - d["Arr"][0]
    return d


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
    Inside = [(x >= xmin and x <= xmax) for x, xmin, xmax in zip(X, XMIN, XMAX)]
    OnLow = [np.abs(x - xmin) < tol for x, xmin in zip(X, XMIN)]
    OnHigh = [np.abs(x - xmax) < tol for x, xmax in zip(X, XMAX)]
    return all([(ins or L or H) for ins, L, H in zip(Inside, OnLow, OnHigh)])


def get_EN_interpolation_at(TMUS, INTERP):
    return np.array([f(np.array(TMUS))[0] for f in INTERP.values()])


def getJacobian(local_TMUS, INTERP, absol_deriv=1e-3, deriv_percent=0.02):
    Jacobian = []
    dts = np.array(
        [
            max(absol_deriv, abs(local_quantity * deriv_percent))
            for local_quantity in local_TMUS
        ]
    )
    L_TMUS = list(local_TMUS)
    for f in INTERP.values():
        Jacobian_row = []
        for i, (lv, hv, dx) in enumerate(zip(local_TMUS - dts, local_TMUS + dts, dts)):
            LV, HV = (
                np.array(L_TMUS[:i] + [lv] + L_TMUS[i + 1 :]),
                np.array(L_TMUS[:i] + [hv] + L_TMUS[i + 1 :]),
            )
            Jacobian_row.append((f(HV)[0] - f(LV)[0]) / (2 * dx))
        Jacobian.append(Jacobian_row)
    return np.array(Jacobian)


def Newton(EN, guessSol, INTERP, MTP):
    iterations, status = 0, True
    TMUs_sol = np.array(guessSol, dtype=float)
    DeltaEN = get_EN_interpolation_at(TMUs_sol, INTERP) - EN
    Newton_criterion_not_meet = any(
        [(np.abs(Delta) > MTP["ACCURACY"]) for Delta in DeltaEN]
    )
    while Newton_criterion_not_meet:
        try:
            J = getJacobian(TMUs_sol, INTERP)
            # if the matrix is very ill-conditioned, use the pseudo-inverse
            cond = np.linalg.cond(J)
            if not np.isfinite(cond) or cond > 1e12:
                invJ = np.linalg.pinv(J)
            else:
                invJ = np.linalg.inv(J)
        except np.linalg.LinAlgError:
            # cannot invert -> abort Newton (fallback will be called outside)
            status = False
            break

        TMUs_sol -= np.dot(invJ, DeltaEN)
        DeltaEN = get_EN_interpolation_at(TMUs_sol, INTERP) - EN
        Newton_criterion_not_meet = any(
            [(np.abs(Delta) > MTP["ACCURACY"]) for Delta in DeltaEN]
        )
        iterations += 1
        if iterations > MTP["MAXITER"]:
            status = False
            break
    status = check_lim(TMUs_sol, MTP["MINS"], MTP["MAXS"]) and status
    return status, TMUs_sol[0], TMUs_sol[1]


def binary_search_1d(ed_local, MUs, INTERP, TMU_TABLES, MTP):
    iteration, LMUs = 0, list(MUs)
    T_min = TMU_TABLES["T"]["MIN"]
    T_max = TMU_TABLES["T"]["MAX"]
    (e_low, e_up) = (INTERP["e"](np.array([T_] + LMUs))[0] for T_ in [T_min, T_max])
    if ed_local < e_low:
        return T_min
    elif ed_local > e_up:
        return T_max
    else:
        while iteration < MTP["MAXITER"]:
            T_mid = (T_max + T_min) / 2.0
            e_mid = INTERP["e"](np.array([T_mid] + LMUs))[0]
            abs_err = abs(e_mid - ed_local)
            rel_err = abs_err / abs(e_mid + ed_local + 1e-15)
            if rel_err <= MTP["ACCURACY"] and abs_err <= MTP["ACCURACY"] * 1e-2:
                return T_mid
            if ed_local < e_mid:
                T_max = T_mid
            else:
                T_min = T_mid
            iteration += 1
        return T_mid


def binary_search_2d(EN, INTERP, TMU_TABLES, MTP):
    iteration = 0
    ed_local, nB_local = EN[0], EN[1]
    muB_min = TMU_TABLES["MUB"]["MIN"]
    muB_max = TMU_TABLES["MUB"]["MAX"]
    (T_min, T_max) = (
        binary_search_1d(ed_local, [muB_], INTERP, TMU_TABLES, MTP)
        for muB_ in [muB_min, muB_max]
    )
    (nB_min, nB_max) = (
        INTERP["nB"](np.array([T_, muB_]))[0]
        for T_, muB_ in zip([T_min, T_max], [muB_min, muB_max])
    )
    if nB_local < nB_min:
        return T_min, muB_min
    elif nB_local > nB_max:
        return T_max, muB_max
    else:
        while iteration < MTP["MAXITER"]:
            muB_mid = (muB_min + muB_max) / 2.0
            T_mid = binary_search_1d(ed_local, [muB_mid], INTERP, TMU_TABLES, MTP)
            nB_mid = INTERP["nB"](np.array([T_mid, muB_mid]))[0]

            abs_err = abs(nB_mid - nB_local)
            rel_err = abs_err / max(1e-15, abs(nB_mid + nB_local))

            if rel_err <= MTP["ACCURACY"] and abs_err <= MTP["ACCURACY"] * 1e-2:
                return T_mid, muB_mid

            if nB_local < nB_mid:
                muB_max = muB_mid
            else:
                muB_min = muB_mid
            iteration += 1
        return T_mid, muB_mid


def invert_EOS_tables(
    T_ID,
    Ttilde,
    iBt,
    TILDE_TABLES,
    TMU_TABLES,
    INTERP,
    MTP,
    hydro_model,
    guessSol=[0.05, 0.0],
):
    mbtilde = TILDE_TABLES["MUB"]["Arr"][iBt]
    hyper_index = key(T_ID, iBt, Nb=TILDE_TABLES["MUB"]["n"])
    EN = ToEN(Ttilde, mbtilde)
    success = False
    try:
        success, T_local, muB_local = Newton(EN, guessSol, INTERP, MTP)
    except Exception:
        # unexpected error in Newton — fall back to backup
        success = False

    if not success or T_local is None:
        T_local, muB_local = binary_search_2d(EN, INTERP, TMU_TABLES, MTP)
        success = True

    (P_local, s_local) = (f([T_local, muB_local])[0] for f in MTP["INTERP_PS"].values())

    # Optional: speed of sound squared and chi2 from forward EoS
    cs2_local = None
    chi2_local = None
    cs2_interp = MTP.get("INTERP_CS2")
    chi2_interp = MTP.get("INTERP_CHI2")
    if cs2_interp is not None:
        try:
            cs2_local = cs2_interp([T_local, muB_local])[0]
        except Exception:
            cs2_local = None
    if chi2_interp is not None:
        try:
            chi2_local = chi2_interp([T_local, muB_local])[0]
        except Exception:
            chi2_local = None
    if "vHLLE" == hydro_model:
        # output for vHLLE: e, nB, T, muB, P, S, cs2, chi2, index
        return [
            EN[0],
            EN[1],
            T_local,
            muB_local,
            P_local,
            s_local,
            cs2_local if cs2_local is not None else np.nan,
            chi2_local if chi2_local is not None else np.nan,
            int(hyper_index),
        ]
    elif "MUSIC" == hydro_model:
        return [
            T_local,
            muB_local,
            P_local,
            s_local,
            int(hyper_index),
        ]  # output for MUSIC: T, muB, P, S, index
    else:
        raise ValueError("Unknown hydro model: {}".format(hydro_model))


if __name__ == "__main__":
    try:
        ParamPath = str(sys.argv[1])
        Ttilde = float(sys.argv[2])
        T_ID = float(sys.argv[3])
        OUTPUT_FOLDER = str(sys.argv[4])
    except IndexError:
        print(
            "Usage: python3 {} EOS_table_file Ttilde ID output folder".format(
                sys.argv[0]
            )
        )
        exit(1)

    import os as _os

    OUTPUT_FOLDER = _os.path.abspath(OUTPUT_FOLDER)
    _os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    output_file = _os.path.join(
        OUTPUT_FOLDER, "TEMP_unordered_inversion_" + str(int(T_ID)) + ".dat"
    )

    Param = read_parameters(ParamPath)
    EoS_table_file = Param["EoS_table"]
    ACCURACY = Param["Accuracy"]
    MAXITER = Param["MAXITER"]
    N_CORES = Param["Number_of_cores"]
    RunMode = Param["RunMode"]
    hydro_model = Param["hydro_model"]
    OutputFolder = Param["OutputFolder"]

    if Param["AutoSetBoundaries"]:
        f = open(OutputFolder + "_boundaries_temp.dat", "rb")
        Boundaries = pickle.load(f)
        TILDE_BOUNDARIES = [Boundaries["Ttilde"], Boundaries["muBtilde"]]
    else:
        TILDE_BOUNDARIES = [
            Param["Ttilde"],
            Param["muBtilde"],
            Param["muQtilde"],
            Param["muStilde"],
        ]

    hbarC = 0.19733

    eos_table = readTable(EoS_table_file)
    Thermodynamic_quantities = ["T", "MUB"]
    Dynamic_quantities = ["e", "nB"]
    Press_Entro = ["P", "S"]
    TMU_TABLES = {
        Quantity_name: extract_DATA_from(eos_table[:, i])
        for i, Quantity_name in enumerate(Thermodynamic_quantities)
    }
    MINS = [TMU_TABLES[Quantity]["MIN"] for Quantity in Thermodynamic_quantities]
    MAXS = [TMU_TABLES[Quantity]["MAX"] for Quantity in Thermodynamic_quantities]

    # EN in GeV**powers
    EN_powers = [4.0, 3.0]

    NT, NB = (
        TMU_TABLES[Q]["n"] for Q in Thermodynamic_quantities
    )  # reshape Temperature
    reshaped_T = TMU_TABLES["T"]["Table"].reshape(NT, NB)
    EN_TABLES = {
        Quantity_name: eos_table[:, i].reshape(NT, NB) * reshaped_T**j
        for i, j, Quantity_name in zip(range(2, 4), EN_powers, Dynamic_quantities)
    }

    PS_powers = [4.0, 3.0]
    PS_TABLES = {
        Quantity_name: eos_table[:, i].reshape(NT, NB) * reshaped_T**j
        for i, j, Quantity_name in zip(range(4, 6), PS_powers, Press_Entro)
    }

    # Grid v (T, muB) prostoru pro všechny interpolátory
    GRID = tuple(TMU_TABLES[Q]["Arr"] for Q in Thermodynamic_quantities)

    # Optional extra tables: cs2 (dimensionless) and chi2(GeV^2)
    INTERP_CS2 = None
    INTERP_CHI2 = None
    if eos_table.shape[1] >= 8:
        CS2_TABLE = eos_table[:, 6].reshape(NT, NB)
        CHI2_over_T2_TABLE = eos_table[:, 7].reshape(NT, NB)
        CHI2_TABLE = CHI2_over_T2_TABLE * reshaped_T**2

        INTERP_CS2 = interpolate.RegularGridInterpolator(
            GRID, CS2_TABLE, bounds_error=False, fill_value=None
        )
        INTERP_CHI2 = interpolate.RegularGridInterpolator(
            GRID, CHI2_TABLE, bounds_error=False, fill_value=None
        )
    INTERP = {
        Quantity_name: interpolate.RegularGridInterpolator(
            GRID, EN_TABLES[Quantity_name], bounds_error=False, fill_value=None
        )
        for Quantity_name in Dynamic_quantities
    }
    INTERP_PS = {
        Quantity_name: interpolate.RegularGridInterpolator(
            GRID, PS_TABLES[Quantity_name], bounds_error=False, fill_value=None
        )
        for Quantity_name in Press_Entro
    }

    META_PARAMS = {
        "ACCURACY": ACCURACY,
        "MAXITER": MAXITER,
        "MINS": MINS,
        "MAXS": MAXS,
        "INTERP_PS": INTERP_PS,
    }

    if INTERP_CS2 is not None:
        META_PARAMS["INTERP_CS2"] = INTERP_CS2
    if INTERP_CHI2 is not None:
        META_PARAMS["INTERP_CHI2"] = INTERP_CHI2

    TILDE_TABLES = {
        Quantity: fill_TILDE(DATA)
        for Quantity, DATA in zip(Thermodynamic_quantities, TILDE_BOUNDARIES)
    }

    B_combinations = list(range(TILDE_TABLES["MUB"]["n"]))

    results = []
    with ProcessPoolExecutor(max_workers=N_CORES) as executor:
        future_to_x = {
            executor.submit(
                invert_EOS_tables,
                T_ID,
                Ttilde,
                iBt,
                TILDE_TABLES,
                TMU_TABLES,
                INTERP,
                META_PARAMS,
                hydro_model,
            ): (iBt)
            for iBt in B_combinations
        }
        for future in as_completed(future_to_x):
            x = future_to_x[future]
            results = future.result()
            with open(output_file, "a") as file:
                file.write(" ".join(map(str, results)) + "\n")
        file.close()
