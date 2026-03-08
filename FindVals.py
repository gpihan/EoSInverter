import os
import sys

import numpy as np
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.ndimage import gaussian_filter1d
from utils import read_parameters


def ToTilde(e, nB):
    Tt = (12 / 19 / np.pi**2 * e) ** (1 / 4)
    mbt = 5 * nB / Tt**2
    return Tt, mbt


def ToEN(Tt, mbt):
    e = 19 / 12 * np.pi**2 * Tt**4
    nb = 1 / 5 * mbt * Tt**2
    return e, nb


def smooth_gaussian_uniform(x, y, sigma_x):
    x = np.asarray(x)
    y = np.asarray(y)

    dx = np.mean(np.diff(x))
    sigma_idx = sigma_x / dx

    y_smooth = gaussian_filter1d(y, sigma=sigma_idx, mode="nearest")
    return y_smooth


def FindVals(EoS, TrLine, RegionS):
    dat = np.loadtxt(EoS)

    NT = len(np.unique(dat[:, 0]))
    NmuB = len(np.unique(dat[:, 1]))
    dat = dat.reshape(NT, NmuB, -1)

    TrLine = np.loadtxt(TrLine)
    # Build interpolation with safe boundary handling to avoid tiny FP overflow
    _x_tr = TrLine[:, 0] / 1000
    _y_tr = TrLine[:, 1] / 1000
    TFUNC = interp1d(
        _x_tr,
        _y_tr,
        bounds_error=False,
        fill_value=(_y_tr[0], _y_tr[-1]),
    )

    shiftT = 0.002

    OUT = []

    for i in range(NmuB):
        arr = dat[:, i]
        T = arr[:, 0]
        muB = arr[:, 1]
        e = arr[:, 2] * T**4
        nb = arr[:, 3] * T**3

        mask = (muB >= 0.4) & (muB <= 0.7)
        Ttrunc = T[mask]
        etrunc = e[mask]
        nbtrunc = nb[mask]

        mask2 = Ttrunc <= TFUNC([muB[0] for i in range(len(Ttrunc))]) - shiftT
        Ttrunc2 = Ttrunc[mask2]
        etrunc2 = etrunc[mask2]
        nbtrunc2 = nbtrunc[mask2]

        if len(Ttrunc2) == 0 or len(etrunc2) == 0:
            continue

        currSpline = UnivariateSpline(Ttrunc2, etrunc2 / Ttrunc2**4, s=0.0)
        currSplinenb = UnivariateSpline(
            Ttrunc2, nbtrunc2 / Ttrunc2**3, s=0.0
        )  # s = 0.0 = linear extrapolation
        # for the hadron gas, the linear extrapolation yields the best result for parallelism and reproduction of eH, nBH

        Tc = TFUNC(muB[0])
        eH = currSpline(Tc) * Tc**4
        nbH = currSplinenb(Tc) * Tc**3

        etrunc_sm = smooth_gaussian_uniform(Ttrunc, etrunc, 0.0008)

        mask4 = Ttrunc >= TFUNC([muB[0] for i in range(len(Ttrunc))]) + 0.002
        Ttrunc4 = Ttrunc[mask4]
        etrunc4 = etrunc_sm[mask4]
        nbtrunc4 = nbtrunc[mask4]

        currSplineA = UnivariateSpline(Ttrunc4, etrunc4 / Ttrunc4**4, s=0.01)
        currSplinenbA = UnivariateSpline(Ttrunc4, nbtrunc4 / Ttrunc4**3, s=0.0)

        eQ = currSplineA(Tc) * Tc**4
        nbQ = currSplinenbA(Tc) * Tc**3

        TtH, mbtH = ToTilde(eH, nbH)
        TtQ, mbtQ = ToTilde(eQ, nbQ)

        OUT.append([TtH, mbtH, TtQ, mbtQ, eH, nbH, eQ, nbQ, Tc, muB[0]])

    OUT = np.array(OUT)

    SortIndexes = np.argsort(OUT[:, 9])
    X, Y = OUT[SortIndexes, 9], OUT[SortIndexes, 2]
    splMBTt = UnivariateSpline(X, Y, s=0.00001)

    SortIndexes = np.argsort(OUT[:, 8])
    X, Y = OUT[SortIndexes, 8], OUT[SortIndexes, 3]
    splTmbt = UnivariateSpline(X, Y, s=0.00008)

    MT = np.array([[splMBTt(mub), splTmbt(T)] for mub, T in zip(OUT[:, 9], OUT[:, 8])])

    OUT[:, 2] = MT[:, 0]
    OUT[:, 3] = MT[:, 1]
    E, N = ToEN(MT[:, 0], MT[:, 1])
    OUT[:, 6] = E
    OUT[:, 7] = N

    with open(RegionS, "w") as f:
        f.write(
            "# TtildeHG muBTildeHG TtildeQGP muBtildeQGP eHG nBHG eQGP nBQGP T muB\n"
        )
        for i in range(len(OUT[:, 0])):
            T, muB, eH, nBH, eQ, nBQ, TtH, mbtH, TtQ, mbtQ = (
                OUT[i, 8],
                OUT[i, 9],
                OUT[i, 4],
                OUT[i, 5],
                OUT[i, 6],
                OUT[i, 7],
                OUT[i, 0],
                OUT[i, 1],
                OUT[i, 2],
                OUT[i, 3],
            )
            f.write(
                f"{TtH:5} {mbtH:5} {TtQ:5} {mbtQ:5} {eH:5} {nBH:5} {eQ:5} {nBQ:5} {T:5} {muB:5}\n"
            )

    print(f"Transition region data written to {RegionS}")


def main(argv):
    try:
        param_path = sys.argv[1]
    except FileNotFoundError:
        print("Parameter file is not found or specified.")
        sys.exit()

    Param = read_parameters(param_path)

    EoS = Param["EoS_table"]

    TrLine = os.path.join(Param["TransitionLine"])
    RegionS = os.path.join(Param["RegionS"])

    print(f"Using EoS table: {EoS}")
    print(f"Using Transition Line data: {TrLine}")
    print(f"Outputting transition region data to: {RegionS}")
    FindVals(EoS, TrLine, RegionS)


if __name__ == "__main__":
    try:
        status = main(sys.argv[1:])
    except Exception as e:
        print(f"Error: {e}")
        print("Invalid arguments.")

    sys.exit()
