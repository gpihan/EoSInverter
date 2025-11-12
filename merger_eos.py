# Copyright Tomáš Poledníček, Gregoire Pihan @ 2025
import os
import sys

import numpy as np
from scipy.interpolate import LinearNDInterpolator, interp1d
from scipy.spatial import Delaunay
from shapely.geometry import MultiPoint

from utils import read_parameters, readTable


def ToEN(Tt, mbt):
    e = 19 * np.pi**2 / 12 * Tt**4
    nb = 1 / 3 * (mbt * Tt**2)
    return e, nb


def Get2DTilde(e, nb):
    Ttilde = (12 / (19 * np.pi**2) * e) ** 0.25
    muBtilde = 5 * nb / Ttilde**2
    return Ttilde, muBtilde


def CheckB(T, mub, TrLine):
    try:
        return (
            (not np.isnan(T))
            and (not np.isnan(mub))
            and (mub >= min_muB_B)
            and (T < TrLine(mub))
        )
    except Exception:
        return False


def CheckA(T, mub, TrLine):
    try:
        return (
            (not np.isnan(T))
            and (not np.isnan(mub))
            and (mub >= min_muB_A)
            and (T > TrLine(mub))
        )
    except Exception:
        return False


def CheckC(T, muB):
    return (
        (not np.isnan(T))
        and (not np.isnan(muB))
        and (min_T_C < T < max_T_C)
        and (min_muB_C < muB < max_muB_C)
    )


def GetGoodTmuB(
    TC,
    TA,
    TB,
    muBC,
    muBA,
    muBB,
    TrLine,
    e,
    nB,
    triA,
    triB,
    triC,
    PC,
    PB,
    PA,
    SC,
    SB,
    SA,
):
    TTilde, muBTilde = Get2DTilde(e, nB)

    p = np.array([TTilde, muBTilde])
    try:
        if triC is not None and (triC.find_simplex(p) >= 0):
            Tc = TC(e, nB)
            muBc = muBC(e, nB)
            Pc = PC(e, nB) if PC is not None else None
            Sc = SC(e, nB) if SC is not None else None
            if CheckC(Tc, muBc):
                return Tc, muBc, Pc, Sc, -1
    except Exception:
        pass

    try:
        if triB is not None and (triB.find_simplex(p) >= 0):
            Tb = TB(e, nB)
            muBb = muBB(e, nB)
            Pb = PB(e, nB) if PB is not None else None
            Sb = SB(e, nB) if SB is not None else None
            if CheckB(Tb, muBb, TrLine):
                return Tb, muBb, Pb, Sb, 0
    except Exception:
        pass

    try:
        if triA is not None and (triA.find_simplex(p) >= 0):
            Ta = TA(e, nB)
            muBa = muBA(e, nB)
            Pa = PA(e, nB) if PA is not None else None
            Sa = SA(e, nB) if SA is not None else None
            if CheckA(Ta, muBa, TrLine):
                return Ta, muBa, Pa, Sa, 1
    except Exception:
        pass

    return False, False, False, False, False


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
            # datS can contain: T muB eH nBH eQ nBQ pH sH pQ sQ
            Tc, muBc, eH, nBH, eQ, nBQ = row[:6]
            # Optional pressure/entropy at hadron/QGP endpoints
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
        return False, False, False, False, False, None, None

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
    # Apply distance thresholds (absolute and normalized)
    if (seg_norm_dist_max is not None) and (norm_d_min > seg_norm_dist_max):
        return False, False, False, False, False, None, None
    chi = float(np.clip(t_min, 0.0, 1.0))  # Use actual projection parameter
    chi_e = None
    chi_nB = None
    # Interpolate pressure and entropy along the segment if available
    P = S = None
    if chi is not None and (pH_f is not None) and (pQ_f is not None):
        P = (1 - chi) * pH_f + chi * pQ_f
    if chi is not None and (sH_f is not None) and (sQ_f is not None):
        S = (1 - chi) * sH_f + chi * sQ_f
    return Tc_f, muBc_f, P, S, chi, chi_e, chi_nB


def create_polygons_tilde(eos_above, eos_below, eos_cross):
    """Create convex hull polygons in (T~, muB~) plane using shapely only (no GeoPandas).

    Returns shapely geometries (poly_above, poly_below, poly_cross).
    """
    TildeT_above, muBtilde_above = Get2DTilde(eos_above[:, 0], eos_above[:, 1])
    TildeT_below, muBtilde_below = Get2DTilde(eos_below[:, 0], eos_below[:, 1])
    TildeT_cross, muBtilde_cross = Get2DTilde(eos_cross[:, 0], eos_cross[:, 1])

    pts_above = MultiPoint(list(zip(TildeT_above, muBtilde_above)))
    pts_below = MultiPoint(list(zip(TildeT_below, muBtilde_below)))
    pts_cross = MultiPoint(list(zip(TildeT_cross, muBtilde_cross)))

    poly_above = pts_above.convex_hull
    poly_below = pts_below.convex_hull
    poly_cross = pts_cross.convex_hull

    return poly_above, poly_below, poly_cross


def create_interpolator(points, values):
    return LinearNDInterpolator(points, values, fill_value=np.nan)


def _make_triangulation(pts):
    try:
        if pts.shape[0] >= 3:
            return Delaunay(pts)
    except Exception:
        return None
    return None


def find_min_max(A, B, C):
    global \
        min_muB_A, \
        max_muB_A, \
        min_muB_B, \
        max_muB_B, \
        min_muB_C, \
        max_muB_C, \
        min_T_C, \
        max_T_C
    min_muB_A, max_muB_A = np.nanmin(A[:, 3]), np.nanmax(A[:, 3])
    min_muB_B, max_muB_B = np.nanmin(B[:, 3]), np.nanmax(B[:, 3])
    min_muB_C, max_muB_C = np.nanmin(C[:, 3]), np.nanmax(C[:, 3])

    min_T_C, max_T_C = np.nanmin(C[:, 2]), np.nanmax(C[:, 2])


def run_merger(
    EoS_above,
    EoS_below,
    EoS_cross,
    TrLine,
    datS,
    out_path,
    Ne,
    Nn,
    show_progress,
    seg_norm_dist_max=0.15,
):
    A = readTable(EoS_above)
    B = readTable(EoS_below)
    C = readTable(EoS_cross)
    datTRLine = readTable(TrLine)
    datS = readTable(
        datS
    )  # T [GeV] muB [GeV] eH [GeV/fm3] nBH [1/fm3] eQ [GeV/fm3] nBQ [1/fm3] pH [GeV/fm3] sH [1/fm3] pQ [GeV/fm3] sQ [1/fm3]

    # define interpolators
    # Transition line (do not extrapolate: return NaN outside domain)
    TrLine = interp1d(
        datTRLine[:, 0] / 1000,
        datTRLine[:, 1] / 1000,
        kind="linear",
        bounds_error=False,
        fill_value=np.nan,
    )  # convert MeV -> GeV, NaN outside

    TC = create_interpolator(C[:, :2], C[:, 2])
    muBC = create_interpolator(C[:, :2], C[:, 3])
    PC = create_interpolator(C[:, :2], C[:, 4]) if C.shape[1] >= 6 else None
    SC = create_interpolator(C[:, :2], C[:, 5]) if C.shape[1] >= 6 else None

    TB = create_interpolator(B[:, :2], B[:, 2])
    muBB = create_interpolator(B[:, :2], B[:, 3])
    PB = create_interpolator(B[:, :2], B[:, 4]) if B.shape[1] >= 6 else None
    SB = create_interpolator(B[:, :2], B[:, 5]) if B.shape[1] >= 6 else None

    TA = create_interpolator(A[:, :2], A[:, 2])
    muBA = create_interpolator(A[:, :2], A[:, 3])
    PA = create_interpolator(A[:, :2], A[:, 4]) if A.shape[1] >= 6 else None
    SA = create_interpolator(A[:, :2], A[:, 5]) if A.shape[1] >= 6 else None

    find_min_max(A, B, C)

    TTilA, muBBtilA = Get2DTilde(A[:, 0], A[:, 1])
    TTilB, muBBtilB = Get2DTilde(B[:, 0], B[:, 1])
    TTilC, muBBtilC = Get2DTilde(C[:, 0], C[:, 1])

    # Precompute points in the tilde plane and Delaunay/KDTree structures for fast membership/nearest queries
    ptsA_tilde = np.column_stack((TTilA, muBBtilA))
    ptsB_tilde = np.column_stack((TTilB, muBBtilB))
    ptsC_tilde = np.column_stack((TTilC, muBBtilC))

    triA = _make_triangulation(ptsA_tilde)
    triB = _make_triangulation(ptsB_tilde)
    triC = _make_triangulation(ptsC_tilde)

    TTil_min = float(min(TTilA.min(), TTilB.min(), TTilC.min()))
    TTil_max = float(max(TTilA.max(), TTilB.max(), TTilC.max()))
    muBBtil_min = float(min(muBBtilA.min(), muBBtilB.min(), muBBtilC.min()))
    muBBtil_max = float(max(muBBtilA.max(), muBBtilB.max(), muBBtilC.max()))

    TTILDEArr = np.linspace(TTil_min, TTil_max, Ne)
    mubTILDEArr = np.linspace(muBBtil_min, muBBtil_max, Nn)

    # Create polygons in tilde plane (for a few geometric checks)
    # poly_above, poly_below, poly_cross = create_polygons_tilde(A, B, C)

    print(f"Writing merged EoS to {out_path} ...")

    with open(out_path, "w") as f:
        f.write("# e nb Ttilde muBtilde T muB P S chi chi_e chi_n\n")
        for i_T, Ttilde in enumerate(TTILDEArr):
            if show_progress and (i_T % max(1, Ne // 10) == 0):
                print(f"Progress: {i_T}/{Ne} rows ({100.0 * i_T / max(1, Ne):.1f}%)")
            for i_mub, muBtilde in enumerate(mubTILDEArr):
                e, nB = ToEN(Ttilde, muBtilde)
                res = find_closest_transition(
                    e,
                    nB,
                    datS,
                    seg_norm_dist_max=seg_norm_dist_max,
                )
                if res[0] is not False:
                    T, muB, P, S, chi, chi_e, chi_n = res
                else:
                    # fallback: use interpolators to decide region

                    T, muB, P, S, chi = GetGoodTmuB(
                        TC,
                        TA,
                        TB,
                        muBC,
                        muBA,
                        muBB,
                        TrLine,
                        e,
                        nB,
                        triA,
                        triB,
                        triC,
                        PC,
                        PB,
                        PA,
                        SC,
                        SB,
                        SA,
                    )
                    chi_e = None
                    chi_n = None
                if T is False:
                    continue

                if muB < 0.02 and T < 0.13:
                    continue
                elif muB < 0.4 and T < 0.075:
                    continue

                chi_e_val = chi_e if chi_e is not None else np.nan
                chi_n_val = chi_n if chi_n is not None else np.nan
                P_val = P if P is not None else np.nan
                S_val = S if S is not None else np.nan

                f.write(
                    f"{e:.6e} {nB:.6e} {Ttilde:.6e} {muBtilde:.6e} {T:.6e} {muB:.6e} {P_val:.6e} {S_val:.6e} {chi:.6e} {chi_e_val:.6e} {chi_n_val:.6e} \n"
                )

    print("... done.")


def main(argv):
    try:
        param_path = sys.argv[1]
    except FileNotFoundError:
        print("Parameter file is not found or specified.")
        sys.exit()

    Param = read_parameters(param_path)

    premerger = Param["premerger_eos"]

    if premerger:
        print("Running premerger EoS cleaning...")
        EoS_above = os.path.join(Param["premerger_output"], "EoSAbove_clean.dat")
        EoS_below = os.path.join(Param["premerger_output"], "EoSBelow_clean.dat")
        EoS_cross = os.path.join(Param["premerger_output"], "EoSCross_clean.dat")
    else:
        EoS_above = os.path.join(Param["EoS_above"])
        EoS_below = os.path.join(Param["EoS_below"])
        EoS_cross = os.path.join(Param["EoS_cross"])
    TrLine = os.path.join(Param["TransitionLine"])
    datS = os.path.join(Param["RegionS"])
    out_path = Param["OutputMergedEoS"]
    Ne = Param["Ne"]
    Nn = Param["Nn"]

    print("Running EoS merger...")
    print("Grid size: ", Ne * Nn)

    run_merger(
        EoS_above,
        EoS_below,
        EoS_cross,
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
        print("Invalid arguments.")

    sys.exit()
