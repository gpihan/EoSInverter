# Copyright Tomáš Poledníček, Gregoire Pihan @ 2025
import argparse
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator, interp1d
from scipy.spatial import Delaunay
from shapely.geometry import Point

"""
Merge three EoS tables (above/below/crossover) into a unified map (e, nB) -> (T, muB, chi).
python3 merger_eos.py \
  --above output/EoSAbove_clean.dat \
  --below output/EoSBelow_clean.dat \
  --cross output/EoSCross_clean.dat \
  --trline input/TransitionLine.dat \
  --regions input/RegionS.dat \
  --out output/MergedEoS.dat \
  --Ne 100 --Nn 100
Inputs (text files):
  - EoS_above.dat   : columns [e, nB, T, muB]
  - EoS_below.dat   : columns [e, nB, T, muB]
  - EoS_cross.dat   : columns [e, nB, T, muB]
  - TransitionLine.dat : columns [muB(MeV), T(MeV)]  — definition of 1st order line
  - RegionS.dat     : columns [Tc(GeV), muBc(GeV), eH, nBH, eQ, nBQ] — from FindVals.py

Output:
    - MergedEoS.dat : columns [e, nB, T, muB, P, S, chi]
      chi = 0.0  -> hadron (below)
      chi = 1.0  -> QGP (above)
      chi = -1.0 -> crossover
      0..1       -> smooth mixture fraction along H–Q transition segment

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

    return p.parse_args(argv)


def readTable(EoS_path):
    df = pd.read_csv(EoS_path, sep="\s+", comment="#", header=None)
    return df.to_numpy()


def ToEN(Tt, mbt):
    e = 19 * np.pi**2 / 12 * Tt**4
    nb = 1 / 5 * (mbt * Tt**2)
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
            and (mub > 0.4)
            and (T < TrLine(mub))
        )
    except Exception:
        return False


def CheckA(T, mub, TrLine):
    try:
        return (
            (not np.isnan(T))
            and (not np.isnan(mub))
            and (mub > 0.4)
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
            Tc, muBc, Pc, Sc = TC(e, nB), muBC(e, nB), PC(e, nB), SC(e, nB)
            if CheckC(Tc, muBc):
                return Tc, muBc, Pc, Sc, -1
    except Exception:
        pass

    try:
        if triB is not None and (triB.find_simplex(p) >= 0):
            Tb, muBb, Pb, Sb = TB(e, nB), muBB(e, nB), PB(e, nB), SB(e, nB)
            if CheckB(Tb, muBb, TrLine):
                return Tb, muBb, Pb, Sb, 0
    except Exception:
        pass

    try:
        if triA is not None and (triA.find_simplex(p) >= 0):
            Ta, muBa, Pa, Sa = TA(e, nB), muBA(e, nB), PA(e, nB), SA(e, nB)
            if CheckA(Ta, muBa, TrLine):
                return Ta, muBa, Pa, Sa, 1
    except Exception:
        pass

    return False, False, False, False, False


def compute_chi(e, nB, eH, nH, eQ, nQ):
    chi_e = None
    chi_n = None
    if eQ != eH:
        chi_e = (e - eH) / (eQ - eH)
    if nQ != nH:
        chi_n = (nB - nH) / (nQ - nH)
    # average if both are valid
    if chi_e is not None and chi_n is not None:
        chi = 0.3  # * (chi_e + chi_n)
    else:
        chi = chi_e if chi_e is not None else chi_n
    return chi, chi_e, chi_n


def point_to_segment_distance(px, py, x1, y1, x2, y2):
    seg = np.array([x2 - x1, y2 - y1])
    pt = np.array([px - x1, py - y1])
    seg_len2 = np.dot(seg, seg)
    if seg_len2 == 0:  # degenerate segment
        return np.linalg.norm(pt)
    t = max(0, min(1, np.dot(pt, seg) / seg_len2))
    proj = np.array([x1, y1]) + t * seg
    return np.linalg.norm([px - proj[0], py - proj[1]])


def find_closest_transition(e, nB, table):
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
        if not (
            (min(eH_f, eQ_f) <= e <= max(eH_f, eQ_f))
            and (min(nBH_f, nBQ_f) <= nB <= max(nBH_f, nBQ_f))
        ):
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
        # T, muB, P, S, chi, chi_e, chi_n
        return False, False, False, False, False, None, None

    # compute distances
    dists = []
    for i, Tc_f, muBc_f, eH_f, nBH_f, eQ_f, nBQ_f, pH_f, sH_f, pQ_f, sQ_f in rows:
        d = point_to_segment_distance(e, nB, eH_f, nBH_f, eQ_f, nBQ_f)
        dists.append(
            (d, i, Tc_f, muBc_f, eH_f, nBH_f, eQ_f, nBQ_f, pH_f, sH_f, pQ_f, sQ_f)
        )

    dists.sort(key=lambda x: x[0])
    (
        _,
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
    chi, chi_e, chi_nB = compute_chi(e, nB, eH_f, nBH_f, eQ_f, nBQ_f)
    # Interpolate pressure and entropy along the segment if available
    P = S = None
    if chi is not None and (pH_f is not None) and (pQ_f is not None):
        P = (1 - chi) * pH_f + chi * pQ_f
    if chi is not None and (sH_f is not None) and (sQ_f is not None):
        S = (1 - chi) * sH_f + chi * sQ_f
    return Tc_f, muBc_f, P, S, chi, chi_e, chi_nB


def create_geodataframes_and_polygons(eos_above, eos_below, eos_cross):
    TildeT_above, muBtilde_above = Get2DTilde(eos_above[:, 0], eos_above[:, 1])
    TildeT_below, muBtilde_below = Get2DTilde(eos_below[:, 0], eos_below[:, 1])
    TildeT_cross, muBtilde_cross = Get2DTilde(eos_cross[:, 0], eos_cross[:, 1])

    geo_df_above = gpd.GeoDataFrame(
        eos_above.copy(), geometry=gpd.points_from_xy(TildeT_above, muBtilde_above)
    )
    geo_df_below = gpd.GeoDataFrame(
        eos_below.copy(), geometry=gpd.points_from_xy(TildeT_below, muBtilde_below)
    )
    geo_df_cross = gpd.GeoDataFrame(
        eos_cross.copy(), geometry=gpd.points_from_xy(TildeT_cross, muBtilde_cross)
    )

    poly_above = geo_df_above.union_all().convex_hull
    poly_below = geo_df_below.union_all().convex_hull
    poly_cross = geo_df_cross.union_all().convex_hull

    poly_below_clean = poly_below

    poly_above_clean = poly_above
    poly_cross_clean = poly_cross

    poly_above_gdf = gpd.GeoDataFrame(geometry=[poly_above_clean])
    poly_below_gdf = gpd.GeoDataFrame(geometry=[poly_below_clean])
    poly_cross_gdf = gpd.GeoDataFrame(geometry=[poly_cross_clean])

    return poly_above_gdf, poly_below_gdf, poly_cross_gdf


def create_interpolator(points, values):
    return LinearNDInterpolator(points, values, fill_value=np.nan)


def run_merger(EoS_above, EoS_below, EoS_cross, TrLine, datS, out_path, Ne, Nn):
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
    PC = create_interpolator(C[:, :2], C[:, 4])
    SC = create_interpolator(C[:, :2], C[:, 5])

    TB = create_interpolator(B[:, :2], B[:, 2])
    muBB = create_interpolator(B[:, :2], B[:, 3])
    PB = create_interpolator(B[:, :2], B[:, 4])
    SB = create_interpolator(B[:, :2], B[:, 5])

    TA = create_interpolator(A[:, :2], A[:, 2])
    muBA = create_interpolator(A[:, :2], A[:, 3])
    PA = create_interpolator(A[:, :2], A[:, 4])
    SA = create_interpolator(A[:, :2], A[:, 5])

    TTilA, muBBtilA = Get2DTilde(A[:, 0], A[:, 1])
    TTilB, muBBtilB = Get2DTilde(B[:, 0], B[:, 1])
    TTilC, muBBtilC = Get2DTilde(C[:, 0], C[:, 1])

    # Precompute points in the tilde plane and Delaunay/KDTree structures for fast membership/nearest queries
    ptsA_tilde = np.column_stack((TTilA, muBBtilA))
    ptsB_tilde = np.column_stack((TTilB, muBBtilB))
    ptsC_tilde = np.column_stack((TTilC, muBBtilC))

    def _make_triangulation(pts):
        try:
            if pts.shape[0] >= 3:
                return Delaunay(pts)
        except Exception:
            return None
        return None

    triA = _make_triangulation(ptsA_tilde)
    triB = _make_triangulation(ptsB_tilde)
    triC = _make_triangulation(ptsC_tilde)

    TTil_min = float(min(TTilA.min(), TTilB.min(), TTilC.min()))
    TTil_max = float(max(TTilA.max(), TTilB.max(), TTilC.max()))
    muBBtil_min = float(min(muBBtilA.min(), muBBtilB.min(), muBBtilC.min()))
    muBBtil_max = float(max(muBBtilA.max(), muBBtilB.max(), muBBtilC.max()))

    TTILDEArr = np.linspace(TTil_min, TTil_max, Ne)
    mubTILDEArr = np.linspace(muBBtil_min, muBBtil_max, Nn)

    # Create GeoDataFrames and polygons
    poly_above_gdf, poly_below_gdf, poly_cross_gdf = create_geodataframes_and_polygons(
        A, B, C
    )

    # plot_T_mub_plane(A, B, C)
    # Main loop - create MergedEoS.dat
    print(f"Writing merged EoS to {out_path} ...")

    with open(out_path, "w") as f:
        f.write("# e nb Ttilde muBtilde T muB P S chi chi_e chi_n\n")
        for Ttilde in TTILDEArr:
            for muBtilde in mubTILDEArr:
                e, nB = ToEN(Ttilde, muBtilde)
                res = find_closest_transition(e, nB, datS)
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
                    point = Point(Get2DTilde(e, nB))
                    if poly_below_gdf.geometry.iloc[0].contains(point) and (
                        chi == 1 or chi == -1
                    ):
                        continue

                    chi_e = None
                    chi_n = None
                    if T is False:
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
    )


if __name__ == "__main__":
    status = main(sys.argv[1:])
    sys.exit(status)
