# Copyright Tomáš Poledníček, Gregoire Pihan @ 2025
import argparse
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator, interp1d
from scipy.spatial import Delaunay
from shapely.geometry import MultiPoint, Point

"""
Merge three EoS tables (above/below/crossover) into a unified map (e, nB) -> thermodynamic fields.

Example:
    python3 merger_eos.py \
        --above output/EoSAbove_clean.dat \
        --below output/EoSBelow_clean.dat \
        --cross output/EoSCross_clean.dat \
        --trline input/TransitionLine.dat \
        --regions input/RegionS.dat \
        --out output/MergedEoS.dat \
        --Ne 100 --Nn 100

Inputs (text files):
    - EoS_above.dat : columns [e, nB, T, muB] with optional [P, S] in columns 5-6
    - EoS_below.dat : columns [e, nB, T, muB] with optional [P, S]
    - EoS_cross.dat : columns [e, nB, T, muB] with optional [P, S]
    - TransitionLine.dat : columns [muB(MeV), T(MeV)] - definition of 1st-order line (loaded in MeV, converted to GeV)
    - RegionS.dat : columns [Tc(GeV), muBc(GeV), eH, nBH, eQ, nBQ] with optional [pH, sH, pQ, sQ]

Output:
    - MergedEoS.dat : columns [e, nB, Ttilde, muBtilde, T, muB, P, S, chi, chi_e, chi_n]
        chi = 0.0  -> hadron (below)
        chi = 1.0  -> QGP (above)
        chi = -1.0 -> crossover
        0..1       -> smooth mixture fraction along H-Q transition segment (projection parameter)

Notes:
    - Temperatures and chemical potentials are in GeV.
    - TransitionLine.dat is loaded in MeV and converted to GeV (division by 1000).
"""


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Merge EoS tables into MergedEoS.dat")
    p.add_argument("--above", default="EoS_above.dat", help="file with 'above' table")
    p.add_argument("--below", default="EoS_below.dat", help="file with 'below' table")
    p.add_argument("--cross", default="EoS_cross.dat", help="file with 'cross' table")
    p.add_argument(
        "--trline", default="TransitionLine.dat", help="file with 1st order line"
    )
    p.add_argument(
        "--regions", default="RegionS.dat", help="file with H–Q segments (FindVals.py)"
    )
    p.add_argument("--out", default="MergedEoS.dat", help="output file")
    p.add_argument("--Ne", type=int, default=400, help="number of steps in e axis")
    p.add_argument("--Nn", type=int, default=200, help="number of steps in nB axis")

    # Misc
    p.add_argument(
        "--progress",
        action="store_true",
        help="Print simple progress updates during merging",
    )

    p.add_argument(
        "--seg-norm-dist-max",
        type=float,
        default=0.15,
        help="Maximum normalized distance d/L (distance to segment divided by its length) to apply RegionS (default: 0.15)",
    )
    return p.parse_args(argv)


def readTable(EoS_path):
    df = pd.read_csv(EoS_path, sep="\s+", comment="#", header=None)
    return df.to_numpy()


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
            and (mub >= 0.4)
            and (T < TrLine(mub))
        )
    except Exception:
        return False


def CheckA(T, mub, TrLine):
    try:
        return (
            (not np.isnan(T))
            and (not np.isnan(mub))
            and (mub >= 0.4)
            and (T > TrLine(mub))
        )
    except Exception:
        return False


def CheckC(T, muB, Tmin=0.0, Tmax=0.6, muBmin=0.0, muBmax=0.399):
    return (
        (not np.isnan(T))
        and (not np.isnan(muB))
        and (Tmin < T < Tmax)
        and (muBmin < muB < muBmax)
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


def run_merger(
    EoS_above,
    EoS_below,
    EoS_cross,
    TrLine,
    datS,
    out_path,
    Ne,
    Nn,
    show_progress=False,
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
    # Optional P,S columns: only build interpolators if data has >= 6 columns
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
    poly_above, poly_below, poly_cross = create_polygons_tilde(A, B, C)

    print(f"Writing merged EoS to {out_path} ...")

    with open(out_path, "w") as f:
        f.write("# e nb Ttilde muBtilde T muB P S chi chi_e chi_n\n")
        for i_T, Ttilde in enumerate(TTILDEArr):
            if show_progress and (i_T % max(1, Ne // 10) == 0):
                print(f"Progress: {i_T}/{Ne} rows ({100.0 * i_T / max(1, Ne):.1f}%)")
            for muBtilde in mubTILDEArr:
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

                point = Point(Get2DTilde(e, nB))
                # If classified as purely hadron or crossover while lying inside the below convex hull,
                # skip as likely inconsistent
                if poly_below.contains(point) and (chi == 1 or chi == -1):
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
    args = parse_args(argv)
    path_above = args.above
    path_below = args.below
    path_cross = args.cross
    path_trline = args.trline
    path_regions = args.regions
    out_path = args.out
    Ne = args.Ne
    Nn = args.Nn

    run_merger(
        path_above,
        path_below,
        path_cross,
        path_trline,
        path_regions,
        out_path,
        Ne,
        Nn,
        show_progress=args.progress,
        seg_norm_dist_max=args.seg_norm_dist_max,
    )


if __name__ == "__main__":
    status = main(sys.argv[1:])
    sys.exit(status)
