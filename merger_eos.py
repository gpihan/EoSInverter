# Copyright Tomáš Poledníček, Gregoire Pihan @2026
import os
import sys
from typing import Tuple

import numpy as np
from scipy.interpolate import LinearNDInterpolator, interp1d

from utils import read_parameters, readTable


def ToEN(Tt: float, mbt: float) -> Tuple[float, float]:
    e = 19 * np.pi**2 / 12 * Tt**4
    nb = 1 / 5 * (mbt * Tt**2)
    return e, nb


def Get2DTilde(e: float, nb: float) -> Tuple[float, float]:
    Ttilde = (12 / (19 * np.pi**2) * e) ** 0.25
    muBtilde = 5 * (nb / Ttilde**2)
    return Ttilde, muBtilde


def GetGoodTmuB(
    T_interp,
    muB_interp,
    P_interp,
    S_interp,
    cs2_interp,
    chi2_interp,
    TrLine,
    e: float,
    nB: float,
):

    T = float(T_interp(e, nB))
    muB = float(muB_interp(e, nB))
    P = float(P_interp(e, nB)) if P_interp is not None else None
    S = float(S_interp(e, nB)) if S_interp is not None else None
    cs2 = float(cs2_interp(e, nB)) if cs2_interp is not None else None
    chi2 = float(chi2_interp(e, nB)) if chi2_interp is not None else None

    region = determine_region(T, muB, TrLine)
    return T, muB, P, S, cs2, chi2, region


def determine_region(T, muB, TrLine):

    if not np.isfinite(T) or not np.isfinite(muB):
        return None

    if muB <= 0.4:
        return -1

    try:
        Tc = float(TrLine(muB))
    except Exception:
        return None

    if not np.isfinite(Tc):
        return None

    if T < Tc:
        return 0
    else:
        return 1


def project_param_and_distance(px, py, x1, y1, x2, y2):
    vx, vy = (x2 - x1), (y2 - y1)
    wx, wy = (px - x1), (py - y1)
    seg_len2 = vx * vx + vy * vy
    if seg_len2 == 0.0:
        return 0.0, float(np.hypot(wx, wy)), 0.0
    t = (wx * vx + wy * vy) / seg_len2
    t = float(np.clip(t, 0.0, 1.0))
    projx = x1 + t * vx
    projy = y1 + t * vy
    d = float(np.hypot(px - projx, py - projy))
    L = float(np.sqrt(seg_len2))
    return t, d, L


def find_closest_transition(e, nB, table, seg_norm_dist_max=10):
    rows = []
    for i, row in enumerate(table):
        try:
            # datS format: T muB eH nBH eQ nBQ pH sH pQ sQ
            Tc, muBc, eH, nBH, eQ, nBQ = row[:6]
            pH = sH = pQ = sQ = None
            if len(row) >= 10:
                pH, sH, pQ, sQ = row[6:10]
            eH_f, nBH_f, eQ_f, nBQ_f = float(eH), float(nBH), float(eQ), float(nBQ)
        except Exception:
            continue
        rows.append(
            (
                i,
                float(Tc),
                float(muBc),
                eH_f,
                nBH_f,
                eQ_f,
                nBQ_f,
                None if pH is None else float(pH),
                None if sH is None else float(sH),
                None if pQ is None else float(pQ),
                None if sQ is None else float(sQ),
            )
        )

    if len(rows) == 0:
        return False, False, False, False, False, False, None

    dists = []
    for i, Tc_f, muBc_f, eH_f, nBH_f, eQ_f, nBQ_f, pH_f, sH_f, pQ_f, sQ_f in rows:
        t, d, L = project_param_and_distance(e, nB, eH_f, nBH_f, eQ_f, nBQ_f)
        norm_d = d / max(L, 1e-12)
        dists.append(
            (
                d,
                norm_d,
                t,
                i,
                Tc_f,
                muBc_f,
                eH_f,
                nBH_f,
                eQ_f,
                nBQ_f,
                pH_f,
                sH_f,
                pQ_f,
                sQ_f,
            )
        )

    dists.sort(key=lambda x: x[0])
    (
        d_min,
        norm_d_min,
        t_min,
        _,
        Tc_f,
        muBc_f,
        eH_f,
        nBH_f,
        eQ_f,
        nBQ_f,
        pH_f,
        sH_f,
        pQ_f,
        sQ_f,
    ) = dists[0]

    if (seg_norm_dist_max is not None) and (norm_d_min > seg_norm_dist_max):
        return False, False, False, False, False, False, None

    chi = float(np.clip(t_min, 0.0, 1))

    P = S = None
    if chi is not None and (pH_f is not None) and (pQ_f is not None):
        P = (1 - chi) * pH_f + chi * pQ_f
    if chi is not None and (sH_f is not None) and (sQ_f is not None):
        S = (1 - chi) * sH_f + chi * sQ_f

    return Tc_f, muBc_f, P, S, None, None, chi


def create_interpolator(points, values):
    return LinearNDInterpolator(points, values, fill_value=np.nan)


def run_merger(
    EoS_table,
    TrLine,
    datS,
    out_path,
    Ne,
    Nn,
    show_progress,
    seg_norm_dist_max=0.15,
):
    """
    Sloučí jednotnou EoS tabulku s daty Region S.

    EoS_table formát:
      e(GeV/fm3) nB(1/fm3) T(GeV) muB(GeV) P(GeV/fm3) s(1/fm3) cs2 chi2(GeV^2) hyper_index
    """
    # Načtení vstupních tabulek
    EoS_data = readTable(EoS_table)
    datTRLine = readTable(TrLine)
    datS = readTable(datS)  # T muB eH nBH eQ nBQ pH sH pQ sQ

    # Sloupce z jednotné EoS tabulky
    e_eos = EoS_data[:, 0]
    nB_eos = EoS_data[:, 1]
    T_eos = EoS_data[:, 2]
    muB_eos = EoS_data[:, 3]
    P_eos = EoS_data[:, 4]
    s_eos = EoS_data[:, 5]
    cs2_eos = EoS_data[:, 6]
    chi2_eos = EoS_data[:, 7]

    muB_axis = datTRLine[:, 0] / 1000  # muB [GeV]
    T_axis = datTRLine[:, 1] / 1000  # T [GeV]
    TrLine = interp1d(
        muB_axis,
        T_axis,
        kind="linear",
        bounds_error=False,
        fill_value=(T_axis[0], T_axis[-1]),
    )

    points_en = np.column_stack((e_eos, nB_eos))

    T_interp = create_interpolator(points_en, T_eos)
    muB_interp = create_interpolator(points_en, muB_eos)
    P_interp = create_interpolator(points_en, P_eos) if P_eos is not None else None
    S_interp = create_interpolator(points_en, s_eos) if s_eos is not None else None
    cs2_interp = (
        create_interpolator(points_en, cs2_eos) if cs2_eos is not None else None
    )
    chi2_interp = (
        create_interpolator(points_en, chi2_eos) if chi2_eos is not None else None
    )

    # Rozsahy v tilde rovině
    TTil_eos, muBBtil_eos = Get2DTilde(e_eos, nB_eos)
    TTil_min = float(TTil_eos.min())
    TTil_max = float(TTil_eos.max())
    muBBtil_min = float(muBBtil_eos.min())
    muBBtil_max = float(muBBtil_eos.max())

    TTILDEArr = np.linspace(TTil_min, TTil_max, Ne)
    mubTILDEArr = np.linspace(muBBtil_min, muBBtil_max, Nn)

    print(f"Writing merged EoS to {out_path} ...")
    print(f"Grid size: {Ne} x {Nn} = {Ne * Nn} points")

    with open(out_path, "w") as f:
        f.write(
            "# Ttilde(GeV) muBtilde(GeV) e(GeV^4) nB(GeV^3) T(GeV) muB(GeV) "
            "P(GeV^4) s(GeV^3) cs2 chi2(GeV^2) chi\n"
        )

        for i_T, Ttilde in enumerate(TTILDEArr):
            if show_progress and (i_T % max(1, Ne // 10) == 0):
                print(f"Progress: {i_T}/{Ne} rows ({100.0 * i_T / max(1, Ne):.1f}%)")

            for i_mub, muBtilde in enumerate(mubTILDEArr):
                e, nB = ToEN(Ttilde, muBtilde)

                T = muB = P = S = np.nan
                cs2_val = chi2_val = np.nan
                chi_val = np.nan

                res = find_closest_transition(
                    e,
                    nB,
                    datS,
                    seg_norm_dist_max=seg_norm_dist_max,
                )

                if res[0] is not False:
                    T_tr, muB_tr, P, S, cs2_val, chi2_val, chi = res

                    if cs2_interp is not None:
                        try:
                            cs2_val = float(cs2_interp(e, nB))
                        except Exception:
                            cs2_val = np.nan
                    if chi2_interp is not None:
                        try:
                            chi2_val = float(chi2_interp(e, nB))
                        except Exception:
                            chi2_val = np.nan

                    try:
                        T_tmp = float(T_interp(e, nB))
                        muB_tmp = float(muB_interp(e, nB))
                    except Exception:
                        T_tmp, muB_tmp = np.nan, np.nan

                    if np.isfinite(T_tmp) and np.isfinite(muB_tmp):
                        T, muB = T_tmp, muB_tmp
                    else:
                        T, muB = float(T_tr), float(muB_tr)

                    if chi is not None:
                        chi_val = float(chi)
                else:
                    (
                        T_raw,
                        muB_raw,
                        P_raw,
                        S_raw,
                        cs2_raw,
                        chi2_raw,
                        region,
                    ) = GetGoodTmuB(
                        T_interp,
                        muB_interp,
                        P_interp,
                        S_interp,
                        cs2_interp,
                        chi2_interp,
                        TrLine,
                        e,
                        nB,
                    )

                    if T_raw is not False and muB_raw is not False:
                        T = float(T_raw)
                        muB = float(muB_raw)
                        P = (
                            float(P_raw)
                            if P_raw is not False and P_raw is not None
                            else np.nan
                        )
                        S = (
                            float(S_raw)
                            if S_raw is not False and S_raw is not None
                            else np.nan
                        )
                        cs2_val = (
                            float(cs2_raw)
                            if cs2_raw is not False and cs2_raw is not None
                            else np.nan
                        )
                        chi2_val = (
                            float(chi2_raw)
                            if chi2_raw is not False and chi2_raw is not None
                            else np.nan
                        )

                        chi_val = float(region) if region is not None else np.nan
                P_val = P if P is not None else np.nan
                S_val = S if S is not None else np.nan
                cs2_val = cs2_val if cs2_val is not None else np.nan
                chi2_val = chi2_val if chi2_val is not None else np.nan

                f.write(
                    f"{Ttilde:.6e} {muBtilde:.6e} {e:.6e} {nB:.6e} {T:.6e} {muB:.6e} "
                    f"{P_val:.6e} {S_val:.6e} {cs2_val:.6e} {chi2_val:.6e} {chi_val:.6e}\n"
                )

    print("... done.")


def main(argv):
    if not argv:
        print("Usage: python merger_single_table.py <parameter_file>")
        sys.exit(1)

    param_path = argv[0]

    if not os.path.exists(param_path):
        print(f"Parameter file not found: {param_path}")
        sys.exit(1)

    Param = read_parameters(param_path)

    EoS_table = os.path.join(Param["OutputFolder"], "EoS_all.dat")

    if not os.path.exists(EoS_table):
        print(f"EoS table not found: {EoS_table}")
        sys.exit(1)

    TrLine = os.path.join(Param["TransitionLine"])
    datS = os.path.join(Param["RegionS"])
    out_path = Param["OutputMergedEoS"]
    Ne = Param["NTildeT"]
    Nn = Param["NTildemuB"]

    print("Running EoS merger with single unified table...")
    print(f"Input EoS: {EoS_table}")
    print(f"Grid size: {Ne} x {Nn} = {Ne * Nn}")

    run_merger(
        EoS_table,
        TrLine,
        datS,
        out_path,
        Ne,
        Nn,
        show_progress=True,
        seg_norm_dist_max=0.15,
    )


if __name__ == "__main__":
    try:
        status = main(sys.argv[1:])
    except Exception as e:
        print(f"Error: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)

    sys.exit(0)
