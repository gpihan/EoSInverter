# Copyright Tomáš Poledníček, Gregoire Pihan @ 2025
import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator, interp1d

"""
Merge three EoS tables (above/below/crossover) into a unified map (e, nB) -> (T, muB, chi).
python3 merger_eos.py \
  --above input/EoSAbove_clean.dat \
  --below input/EoSBelow_clean.dat \
  --cross input/EoS_cross.dat \
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
  - MergedEoS.dat : columns [e, nB, T, muB, chi]
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
    nb = 1 / 3 * mbt * Tt**2
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


def CheckC(T, muB, Tmin=0.029, Tmax=0.4, muBmin=0.0, muBmax=0.4):
    return (
        (not np.isnan(T))
        and (not np.isnan(muB))
        and (Tmin < T < Tmax)
        and (muBmin < muB < muBmax)
    )


def GetGoodTmuB(TC, TA, TB, muBC, muBA, muBB, TrLine, e, nB, A, B, C):
    e_A, n_A = A[:, 0], A[:, 1]
    e_B, n_B = B[:, 0], B[:, 1]
    e_C, n_C = C[:, 0], C[:, 1]

    # check if inside convex hulls
    if (e > e_A.min() and e < e_A.max()) and (nB > n_A.min() and nB < n_A.max()):
        Ta, muBa = TA(e, nB), muBA(e, nB)
    if (e > e_B.min() and e < e_B.max()) and (nB > n_B.min() and nB < n_B.max()):
        Tb, muBb = TB(e, nB), muBB(e, nB)
    if (e > e_C.min() and e < e_C.max()) and (nB > n_C.min() and nB < n_C.max()):
        Tc, muBc = TC(e, nB), muBC(e, nB)

    if (
        (e > e_A.min() and e < e_A.max())
        and (nB > n_A.min() and nB < n_A.max())
        and CheckA(Ta, muBa, TrLine)
    ):
        return Ta, muBa, 1
    if (
        (e > e_B.min() and e < e_B.max())
        and (nB > n_B.min() and nB < n_B.max())
        and CheckB(Tb, muBb, TrLine)
    ):
        return Tb, muBb, 0
    if (
        (e > e_C.min() and e < e_C.max())
        and (nB > n_C.min() and nB < n_C.max())
        and CheckC(Tc, muBc)
    ):
        return Tc, muBc, -1

    return False, False, False


def compute_chi(e, nB, eH, nH, eQ, nQ):
    chi_e = None
    chi_n = None
    if eQ != eH:
        chi_e = (e - eH) / (eQ - eH)
    if nQ != nH:
        chi_n = (nB - nH) / (nQ - nH)
    # average if both are valid
    if chi_e is not None and chi_n is not None:
        chi = 0.5 * (chi_e + chi_n)
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
            Tc, muBc, eH, nBH, eQ, nBQ = row
            eH_f, nBH_f, eQ_f, nBQ_f = float(eH), float(nBH), float(eQ), float(nBQ)
        except Exception:
            continue
        if not (
            (min(eH_f, eQ_f) <= e <= max(eH_f, eQ_f))
            and (min(nBH_f, nBQ_f) <= nB <= max(nBH_f, nBQ_f))
        ):
            continue
        rows.append((i, float(Tc), float(muBc), eH_f, nBH_f, eQ_f, nBQ_f))

    if len(rows) == 0:
        return False, False, False, None, None

    # compute distances
    dists = []
    for i, Tc_f, muBc_f, eH_f, nBH_f, eQ_f, nBQ_f in rows:
        d = point_to_segment_distance(e, nB, eH_f, nBH_f, eQ_f, nBQ_f)
        dists.append((d, i, Tc_f, muBc_f, eH_f, nBH_f, eQ_f, nBQ_f))

    dists.sort(key=lambda x: x[0])
    _, _, Tc_f, muBc_f, eH_f, nBH_f, eQ_f, nBQ_f = dists[0]
    chi, chi_e, chi_nB = compute_chi(e, nB, eH_f, nBH_f, eQ_f, nBQ_f)
    return Tc_f, muBc_f, chi, chi_e, chi_nB


def plot_T_mub_plane(A, B, C):
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(A[:, 3], A[:, 2], s=1, color="red", label="Above")
    ax.scatter(B[:, 3], B[:, 2], s=1, color="blue", label="Below")
    ax.scatter(C[:, 3], C[:, 2], s=1, color="green", label="Cross")
    ax.set_xlabel("muB [GeV]")
    ax.set_ylabel("T [GeV]")
    ax.legend()
    plt.show()


def plot_tilde_plane(
    TTilA, muBBtilA, TTilB, muBBtilB, TTilC, muBBtilC, TTILDEArr, mubTILDEArr
):
    TT_mesh, mub_mesh = np.meshgrid(TTILDEArr, mubTILDEArr, indexing="xy")
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(
        TT_mesh.ravel(),
        mub_mesh.ravel(),
        s=0.3,
        color="black",
        alpha=0.1,
        label="Grid points",
    )
    ax.scatter(TTilA, muBBtilA, s=1, color="red", label="Above")
    ax.scatter(TTilB, muBBtilB, s=1, color="blue", label="Below")
    ax.scatter(TTilC, muBBtilC, s=1, color="green", label="Cross")
    ax.set_xlabel("Ttilde [GeV]")
    ax.set_ylabel("muBtilde [GeV]")
    ax.legend()
    plt.show()


def run_merger(
    EoS_above, EoS_below, EoS_cross, TrLine, datS, out_path, Ne, Nn, no_plot=False
):
    A = readTable(EoS_above)
    B = readTable(EoS_below)
    C = readTable(EoS_cross)
    datTRLine = readTable(TrLine)
    datS = readTable(datS)

    # define interpolators
    # Transition line
    TrLine = interp1d(
        datTRLine[:, 0] / 1000, datTRLine[:, 1] / 1000
    )  # convert MeV -> GeV

    # Cross
    pointsC = C[:, :2]
    valuesC = C[:, 2]
    TC = LinearNDInterpolator(pointsC, valuesC, fill_value=np.nan)
    valuesMC = C[:, 3]
    muBC = LinearNDInterpolator(pointsC, valuesMC, fill_value=np.nan)

    # Below
    pointsB = B[:, :2]
    valuesB = B[:, 2]
    TB = LinearNDInterpolator(pointsB, valuesB, fill_value=np.nan)
    valuesMB = B[:, 3]
    muBB = LinearNDInterpolator(pointsB, valuesMB, fill_value=np.nan)

    # Above
    pointsA = A[:, :2]
    valuesA = A[:, 2]
    TA = LinearNDInterpolator(pointsA, valuesA, fill_value=np.nan)
    valuesMA = A[:, 3]
    muBA = LinearNDInterpolator(pointsA, valuesMA, fill_value=np.nan)

    TTilA, muBBtilA = Get2DTilde(A[:, 0], A[:, 1])
    TTilB, muBBtilB = Get2DTilde(B[:, 0], B[:, 1])
    TTilC, muBBtilC = Get2DTilde(C[:, 0], C[:, 1])

    TTil_min = float(min(TTilA.min(), TTilB.min(), TTilC.min()))
    TTil_max = float(max(TTilA.max(), TTilB.max(), TTilC.max()))
    muBBtil_min = float(min(muBBtilA.min(), muBBtilB.min(), muBBtilC.min()))
    muBBtil_max = float(max(muBBtilA.max(), muBBtilB.max(), muBBtilC.max()))

    TTILDEArr = np.linspace(TTil_min, TTil_max, Ne)
    mubTILDEArr = np.linspace(muBBtil_min, muBBtil_max, Nn)

    # plot_T_mub_plane(A, B, C)
    # Main loop - create MergedEoS.dat
    print(f"Writing merged EoS to {out_path} ...")

    with open(out_path, "w") as f:
        f.write("# e nb Ttilde muBtilde T muB chi chi_e chi_n\n")
        for Ttilde in TTILDEArr:
            for muBtilde in mubTILDEArr:
                e, nB = ToEN(Ttilde, muBtilde)
                res = find_closest_transition(e, nB, datS)
                if res[0] is not False:
                    T, muB, chi, chi_e, chi_n = res
                else:
                    # fallback: use interpolators to decide region
                    T, muB, chi = GetGoodTmuB(
                        TC, TA, TB, muBC, muBA, muBB, TrLine, e, nB, A, B, C
                    )
                    chi_e = None
                    chi_n = None
                    if T is False:
                        continue
                chi_e_val = chi_e if chi_e is not None else np.nan
                chi_n_val = chi_n if chi_n is not None else np.nan

                f.write(
                    f"{e:.6e} {nB:.6e} {Ttilde:.6e} {muBtilde:.6e} {T:.6e} {muB:.6e} {chi:.6e} {chi_e_val:.6e} {chi_n_val:.6e}\n"
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
