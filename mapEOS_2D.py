# /usr/bin/env python3
# Copyright Gregoire Pihan @ 2024

import os
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
    return np.array([e, nb], dtype=float)


def fill_LIN_AXIS(DATA):
    n = int(DATA[2])
    vmin = float(DATA[0])
    vmax = float(DATA[1])
    arr = np.linspace(vmin, vmax, n)
    d = {
        "n": n,
        "MIN": vmin,
        "MAX": vmax,
        "Arr": arr,
        "d": (arr[1] - arr[0]) if n > 1 else 0.0,
    }
    return d


def fill_LOG_AXIS(DATA, *, include_zero=True, min_pos=1e-12):
    vmin, vmax, n = float(DATA[0]), float(DATA[1]), int(DATA[2])

    if n < 2:
        arr = np.array([vmin], dtype=float)
        return {
            "n": 1,
            "MIN": float(arr[0]),
            "MAX": float(arr[-1]),
            "Arr": arr,
            "d": 0.0,
        }

    if vmax <= 0 or not np.isfinite(vmax):
        return fill_LIN_AXIS(DATA)

    min_pos = float(min_pos)
    if not np.isfinite(min_pos) or min_pos <= 0:
        min_pos = vmax * 1e-12
    min_pos = min(min_pos, vmax * 0.999999)

    if include_zero:
        if n == 2:
            arr = np.array([0.0, vmax], dtype=float)
        else:
            arr_log = np.logspace(np.log10(min_pos), np.log10(vmax), n - 1)
            arr = np.concatenate(([0.0], arr_log))
    else:
        arr = np.logspace(np.log10(min_pos), np.log10(vmax), n)

    return {
        "n": int(len(arr)),
        "MIN": float(arr[0]),
        "MAX": float(arr[-1]),
        "Arr": arr,
        "d": np.nan,
    }


def extract_axis(col):
    arr = np.asarray(col, dtype=float)
    uniq = np.unique(arr)
    return {
        "n": int(len(uniq)),
        "MIN": float(uniq[0]),
        "MAX": float(uniq[-1]),
        "Table": arr,
        "Arr": uniq,
    }


def check_lim(X, XMIN, XMAX, tol=1e-10):
    Inside = [(x >= xmin and x <= xmax) for x, xmin, xmax in zip(X, XMIN, XMAX)]
    OnLow = [np.abs(x - xmin) < tol for x, xmin in zip(X, XMIN)]
    OnHigh = [np.abs(x - xmax) < tol for x, xmax in zip(X, XMAX)]
    return all([(ins or L or H) for ins, L, H in zip(Inside, OnLow, OnHigh)])


def get_EN_interpolation_at(TMUS, INTERP):
    TMUS = np.asarray(TMUS, dtype=float)
    return np.array([f(TMUS)[0] for f in INTERP.values()], dtype=float)


def getJacobian(local_TMUS, INTERP, absol_deriv=1e-4, deriv_percent=0.02):
    Jacobian = []
    local_TMUS = np.asarray(local_TMUS, dtype=float)
    dts = np.array(
        [max(absol_deriv, abs(q * deriv_percent)) for q in local_TMUS], dtype=float
    )
    L_TMUS = list(local_TMUS)
    for f in INTERP.values():
        Jacobian_row = []
        for i, (lv, hv, dx) in enumerate(zip(local_TMUS - dts, local_TMUS + dts, dts)):
            LV = np.array(L_TMUS[:i] + [lv] + L_TMUS[i + 1 :], dtype=float)
            HV = np.array(L_TMUS[:i] + [hv] + L_TMUS[i + 1 :], dtype=float)
            Jacobian_row.append((f(HV)[0] - f(LV)[0]) / (2.0 * dx))
        Jacobian.append(Jacobian_row)
    return np.asarray(Jacobian, dtype=float)


def Newton(EN, guessSol, INTERP, MTP):
    iterations, status = 0, True
    TMUs_sol = np.array(guessSol, dtype=float)
    DeltaEN = get_EN_interpolation_at(TMUs_sol, INTERP) - EN
    while np.linalg.norm(DeltaEN) > MTP["ACCURACY"]:
        try:
            J = getJacobian(TMUs_sol, INTERP)
            cond = np.linalg.cond(J)
            if not np.isfinite(cond):
                status = False
                break
            if cond > 1e12:
                step = np.dot(np.linalg.pinv(J), DeltaEN)
            else:
                step = np.linalg.solve(J, DeltaEN)
        except np.linalg.LinAlgError:
            status = False
            break

        if np.any(~np.isfinite(step)):
            status = False
            break

        TMUs_sol -= step
        DeltaEN = get_EN_interpolation_at(TMUs_sol, INTERP) - EN
        iterations += 1
        if iterations > MTP["MAXITER"]:
            status = False
            break

    status = status and check_lim(TMUs_sol, MTP["MINS"], MTP["MAXS"])
    return status, TMUs_sol[0], TMUs_sol[1]


def binary_search_1d(e_target, muB, INTERP, TMU_TABLES, MTP):
    T_min = TMU_TABLES["T"]["MIN"]
    T_max = TMU_TABLES["T"]["MAX"]

    e_low = INTERP["e"]([T_min, muB])[0]
    e_up = INTERP["e"]([T_max, muB])[0]

    if e_target <= e_low:
        return T_min
    if e_target >= e_up:
        return T_max

    T_mid = 0.5 * (T_min + T_max)
    for _ in range(MTP["MAXITER"]):
        T_mid = 0.5 * (T_min + T_max)
        e_mid = INTERP["e"]([T_mid, muB])[0]
        abs_err = abs(e_mid - e_target)
        rel_err = abs_err / (abs(e_mid + e_target) + 1e-15)
        if rel_err <= MTP["ACCURACY"] and abs_err <= MTP["ACCURACY"] * 1e-2:
            return T_mid
        if e_target < e_mid:
            T_max = T_mid
        else:
            T_min = T_mid
    return T_mid


def nB_reachable_range_for_e(e_target, INTERP, TMU_TABLES, MTP):
    muB_min = TMU_TABLES["MUB"]["MIN"]
    muB_max = TMU_TABLES["MUB"]["MAX"]
    T_at_min = binary_search_1d(e_target, muB_min, INTERP, TMU_TABLES, MTP)
    T_at_max = binary_search_1d(e_target, muB_max, INTERP, TMU_TABLES, MTP)
    nB_min = float(INTERP["nB"]([T_at_min, muB_min])[0])
    nB_max = float(INTERP["nB"]([T_at_max, muB_max])[0])
    return nB_min, nB_max, float(T_at_min), float(T_at_max)


def binary_search_2d(EN, INTERP, TMU_TABLES, MTP):
    e_target, nB_target = float(EN[0]), float(EN[1])
    muB_min = TMU_TABLES["MUB"]["MIN"]
    muB_max = TMU_TABLES["MUB"]["MAX"]

    nB_min, nB_max, T_at_min, T_at_max = nB_reachable_range_for_e(
        e_target, INTERP, TMU_TABLES, MTP
    )
    if nB_target <= nB_min:
        return T_at_min, muB_min
    if nB_target >= nB_max:
        return T_at_max, muB_max

    mu_lo, mu_hi = muB_min, muB_max
    T_mid, mu_mid = T_at_min, 0.5 * (mu_lo + mu_hi)
    for _ in range(MTP["MAXITER"]):
        mu_mid = 0.5 * (mu_lo + mu_hi)
        T_mid = binary_search_1d(e_target, mu_mid, INTERP, TMU_TABLES, MTP)
        nB_mid = float(INTERP["nB"]([T_mid, mu_mid])[0])

        abs_err = abs(nB_mid - nB_target)
        rel_err = abs_err / (abs(nB_mid + nB_target) + 1e-15)
        if rel_err <= MTP["ACCURACY"] and abs_err <= MTP["ACCURACY"] * 1e-2:
            return float(T_mid), float(mu_mid)

        if nB_target < nB_mid:
            mu_hi = mu_mid
        else:
            mu_lo = mu_mid

    return float(T_mid), float(mu_mid)


def invert_from_EN(
    E_ID,
    e_value,
    iBn,
    EN_AXES,
    TMU_TABLES,
    INTERP,
    MTP,
    hydro_model,
    guessSol=(0.05, 0.0),
):
    nB_value = float(EN_AXES["NB"]["Arr"][iBn])
    hyper_index = key(E_ID, iBn, Nb=EN_AXES["NB"]["n"])
    EN = np.array([float(e_value), nB_value], dtype=float)

    # Optional realizability clamp
    if MTP.get("ENClamp", True):
        nB_min, nB_max, _, _ = nB_reachable_range_for_e(EN[0], INTERP, TMU_TABLES, MTP)
        if EN[1] < nB_min:
            EN[1] = nB_min
        elif EN[1] > nB_max:
            EN[1] = nB_max

    success = False
    try:
        success, T_local, muB_local = Newton(EN, guessSol, INTERP, MTP)
    except Exception:
        success = False

    if not success or T_local is None or muB_local is None:
        T_local, muB_local = binary_search_2d(EN, INTERP, TMU_TABLES, MTP)

    P_local = float(MTP["INTERP_PS"]["P"]([T_local, muB_local])[0])
    s_local = float(MTP["INTERP_PS"]["S"]([T_local, muB_local])[0])

    cs2_local = np.nan
    chi2_local = np.nan
    if MTP.get("INTERP_CS2") is not None:
        try:
            cs2_local = float(MTP["INTERP_CS2"]([T_local, muB_local])[0])
        except Exception:
            cs2_local = np.nan
    if MTP.get("INTERP_CHI2") is not None:
        try:
            chi2_local = float(MTP["INTERP_CHI2"]([T_local, muB_local])[0])
        except Exception:
            chi2_local = np.nan

    if "vHLLE" == hydro_model:
        return [
            EN[0],
            EN[1],
            T_local,
            muB_local,
            P_local,
            s_local,
            cs2_local,
            chi2_local,
            int(hyper_index),
        ]
    if "MUSIC" == hydro_model:
        return [T_local, muB_local, P_local, s_local, int(hyper_index)]
    raise ValueError("Unknown hydro model: {}".format(hydro_model))


def diag_once(TMU_TABLES, EN_TABLES_FORWARD, INTERP, EN_AXES=None):
    T_arr = TMU_TABLES["T"]["Arr"]
    mu_arr = TMU_TABLES["MUB"]["Arr"]
    Tmin, Tmax = float(T_arr[0]), float(T_arr[-1])
    mumin, mumax = float(mu_arr[0]), float(mu_arr[-1])
    print("=== DIAG: axis ranges from table (unique) ===")
    print(f"T:   n={TMU_TABLES['T']['n']}   min={Tmin}   max={Tmax}")
    print(f"muB: n={TMU_TABLES['MUB']['n']} min={mumin}  max={mumax}")
    print("=== DIAG: forward EN min/max ===")
    print(
        f"e:  min={np.nanmin(EN_TABLES_FORWARD['e'])} max={np.nanmax(EN_TABLES_FORWARD['e'])}"
    )
    print(
        f"nB: min={np.nanmin(EN_TABLES_FORWARD['nB'])} max={np.nanmax(EN_TABLES_FORWARD['nB'])}"
    )
    print("=== DIAG: corners ===")
    corners = [(Tmin, mumin), (Tmin, mumax), (Tmax, mumin), (Tmax, mumax)]
    for T, mu in corners:
        e = float(INTERP["e"]([T, mu])[0])
        nB = float(INTERP["nB"]([T, mu])[0])
        print(f"(T={T}, muB={mu}) -> e={e}, nB={nB}")

    if EN_AXES is not None:
        nb = EN_AXES["NB"]["Arr"]
        print("=== DIAG: EN NB axis ===")
        print("NB first 10:", nb[:10])
        print("NB last 10 :", nb[-10:])


if __name__ == "__main__":
    try:
        ParamPath = str(sys.argv[1])
        X_value = float(sys.argv[2])  # Ttilde if tilde_switch else e_value
        X_ID = int(float(sys.argv[3]))
        OUTPUT_FOLDER = str(sys.argv[4])
    except IndexError:
        print(
            "Usage: python3 {} parameters.py X_value ID output_folder".format(
                sys.argv[0]
            )
        )
        sys.exit(1)

    Param = read_parameters(ParamPath)
    EoS_table_file = Param["EoS_table"]
    ACCURACY = float(Param["Accuracy"])
    MAXITER = int(Param["MAXITER"])
    N_CORES = int(Param["Number_of_cores"])
    hydro_model = Param["hydro_model"]
    OutputFolder = Param["OutputFolder"]
    tilde_switch = Param.get("TildeSwitch", True)

    # control clamp in EN-mode
    ENClamp = bool(Param.get("ENClamp", True))

    OUTPUT_FOLDER = os.path.abspath(OUTPUT_FOLDER)
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    out_path = os.path.join(
        OUTPUT_FOLDER, "TEMP_unordered_inversion_" + str(int(X_ID)) + ".dat"
    )

    # boundaries
    if Param.get("AutoSetBoundaries", False):
        candidates = [OutputFolder + "_boundaries_temp.dat", "boundaries_temp.dat"]
        Boundaries = None
        for c in candidates:
            try:
                with open(c, "rb") as f:
                    Boundaries = pickle.load(f)
                break
            except Exception:
                continue
        if Boundaries is None:
            raise FileNotFoundError(
                "Cannot read boundaries file. Tried: " + ", ".join(candidates)
            )

        if tilde_switch:
            TILDE_BOUNDARIES = [Boundaries["Ttilde"], Boundaries["muBtilde"]]
        else:
            EN_BOUNDARIES = [Boundaries["e"], Boundaries["nB"]]
    else:
        if tilde_switch:
            TILDE_BOUNDARIES = [Param["Ttilde"], Param["muBtilde"]]
        else:
            EN_BOUNDARIES = [Param["e"], Param["nB"]]

    eos_table = readTable(EoS_table_file)

    Thermodynamic_quantities = ["T", "MUB"]
    Dynamic_quantities = ["e", "nB"]
    Press_Entro = ["P", "S"]

    TMU_TABLES = {
        q: extract_axis(eos_table[:, i]) for i, q in enumerate(Thermodynamic_quantities)
    }
    MINS = [TMU_TABLES[q]["MIN"] for q in Thermodynamic_quantities]
    MAXS = [TMU_TABLES[q]["MAX"] for q in Thermodynamic_quantities]

    NT, NB = (TMU_TABLES[q]["n"] for q in Thermodynamic_quantities)
    if eos_table.shape[0] != NT * NB:
        raise ValueError(
            f"Table rows={eos_table.shape[0]} but NT*NB={NT * NB}. "
            "Table must be rectangular and ordered for reshape."
        )

    reshaped_T = eos_table[:, 0].reshape(NT, NB)

    EN_powers = [4.0, 3.0]
    EN_TABLES_FORWARD = {
        q: eos_table[:, i].reshape(NT, NB) * reshaped_T**p
        for i, p, q in zip(range(2, 4), EN_powers, Dynamic_quantities)
    }

    PS_powers = [4.0, 3.0]
    PS_TABLES = {
        q: eos_table[:, i].reshape(NT, NB) * reshaped_T**p
        for i, p, q in zip(range(4, 6), PS_powers, Press_Entro)
    }

    GRID = (TMU_TABLES["T"]["Arr"], TMU_TABLES["MUB"]["Arr"])
    INTERP = {
        q: interpolate.RegularGridInterpolator(
            GRID, EN_TABLES_FORWARD[q], bounds_error=False, fill_value=None
        )
        for q in Dynamic_quantities
    }
    INTERP_PS = {
        q: interpolate.RegularGridInterpolator(
            GRID, PS_TABLES[q], bounds_error=False, fill_value=None
        )
        for q in Press_Entro
    }

    # optional extra
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

    META_PARAMS = {
        "ACCURACY": ACCURACY,
        "MAXITER": MAXITER,
        "MINS": MINS,
        "MAXS": MAXS,
        "INTERP_PS": {"P": INTERP_PS["P"], "S": INTERP_PS["S"]},
        "ENClamp": ENClamp,
    }
    if INTERP_CS2 is not None:
        META_PARAMS["INTERP_CS2"] = INTERP_CS2
    if INTERP_CHI2 is not None:
        META_PARAMS["INTERP_CHI2"] = INTERP_CHI2

    if tilde_switch:
        raise RuntimeError(
            "This file version focuses on EN-mode log nB. Set TildeSwitch=False for this test."
        )

    # EN axes: E linear, NB log
    e_value = float(X_value)

    nB_forward = EN_TABLES_FORWARD["nB"].ravel()
    nB_pos = nB_forward[np.isfinite(nB_forward) & (nB_forward > 0)]
    # small quantile to hit low-T tiny densities
    nB_min_pos = float(np.quantile(nB_pos, 0.001)) if nB_pos.size else 1e-12

    EN_AXES = {
        "E": fill_LIN_AXIS(EN_BOUNDARIES[0]),
        "NB": fill_LOG_AXIS(EN_BOUNDARIES[1], include_zero=True, min_pos=nB_min_pos),
    }

    if X_ID == 0:
        diag_once(TMU_TABLES, EN_TABLES_FORWARD, INTERP, EN_AXES=EN_AXES)
        print("Using ENClamp =", ENClamp)
        print("Using nB_min_pos =", nB_min_pos)
        print("EN e boundary:", EN_BOUNDARIES[0])
        print("EN nB boundary:", EN_BOUNDARIES[1])

    nB_idxs = list(range(EN_AXES["NB"]["n"]))
    with ProcessPoolExecutor(max_workers=N_CORES) as executor:
        futures = {
            executor.submit(
                invert_from_EN,
                X_ID,
                e_value,
                iBn,
                EN_AXES,
                TMU_TABLES,
                INTERP,
                META_PARAMS,
                hydro_model,
            ): iBn
            for iBn in nB_idxs
        }
        for fut in as_completed(futures):
            res = fut.result()
            with open(out_path, "a") as f:
                f.write(" ".join(map(str, res)) + "\n")
