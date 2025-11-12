import os
import sys

import numpy as np
from scipy.interpolate import interp1d

from utils import read_parameters, readTable


def first_index_zero(lst, tol: float = 1e-12):
    arr = np.asarray(lst)
    idx = np.where(np.abs(arr) <= tol)[0]
    return int(idx[0]) if idx.size else None


def first_index_non_zero(lst, tol: float = 1e-12):
    arr = np.asarray(lst)
    idx = np.where(np.abs(arr) > tol)[0]
    return int(idx[0]) if idx.size else None


def deriv(vec):
    return np.diff(vec)


def GetFirst(lst):
    arr = np.asarray(lst)
    diffs = np.abs(arr - arr[0])
    return first_index_non_zero(diffs, tol=1e-12)


def infer_grid_shape(arr, t_col=0, mub_col=1, tol=1e-10):
    nrows = arr.shape[0]
    T = arr[:, t_col]
    muB = arr[:, mub_col]

    diffs = np.abs(np.diff(T))
    brk = np.where(diffs > tol)[0]
    NmuB_run = brk[0] + 1 if brk.size else nrows

    if nrows % NmuB_run == 0:
        NT = nrows // NmuB_run
        return NT, NmuB_run

    uniq_muB = np.unique(np.round(muB, 12))
    NmuB_uniq = uniq_muB.size
    if nrows % NmuB_uniq == 0:
        NT = nrows // NmuB_uniq
        return NT, NmuB_uniq

    raise ValueError(
        f"Cannot infer grid: nrows={nrows}, NmuB_run={NmuB_run} not dividing, unique muB={NmuB_uniq} not dividing."
    )


def main(argv):
    try:
        param_path = sys.argv[1]
    except FileNotFoundError:
        print("Parameter file is not found or specified.")
        sys.exit()

    Param = read_parameters(param_path)
    script_dir = os.path.dirname(os.path.abspath(__file__))

    EoS_A = os.path.join(script_dir, Param["EoS_above_input"])
    EoS_B = os.path.join(script_dir, Param["EoS_below_input"])
    TrLine = os.path.join(script_dir, Param["TransitionLine"])

    print("Running FindVals.py to generate RegionS.dat...")
    print(f"EoS above: {EoS_A}")
    print(f"EoS below: {EoS_B}")
    print(f"Transition line: {TrLine}")

    Trline = readTable(TrLine)
    TL = interp1d(Trline[:, 0] / 1000, Trline[:, 1] / 1000)

    EoSB = readTable(EoS_B)

    EoSA = readTable(EoS_A)
    try:
        NTB, NmuB = infer_grid_shape(EoSB, t_col=0, mub_col=1)

        EoSB = EoSB.reshape(NTB, NmuB, -1)
    except ValueError:
        print("Error inferring grid for EoSBelowGeV.dat:")
        raise

    TListH = []
    muBListH = []
    eHList = []
    nBHList = []
    # new: pressure and entropy at H side (below)
    pHList = []
    sHList = []

    for imuB in range(NmuB):
        T = EoSB[:, imuB][:, 0]
        eps = EoSB[:, imuB][:, 2]

        muB_val = EoSB[:, imuB][:, 1][0]

        ind = len(eps) - GetFirst(eps[::-1]) - 1
        eH = eps[:ind][-1] * T[:ind][-1] ** 4
        eHList.append(eH)

        nB = EoSB[:, imuB][:, 3]
        ind = len(nB) - GetFirst(nB[::-1])
        nBH = nB[:ind][-1] * T[:ind][-1] ** 3
        nBHList.append(nBH)

        # pressure p/T^4 -> p [GeV^4]
        p = EoSB[:, imuB][:, 5]
        ind_p = len(p) - GetFirst(p[::-1]) - 1
        pH = p[:ind_p][-1] * T[:ind_p][-1] ** 4
        pHList.append(pH)

        # entropy s/T^3 -> s [GeV^3]
        s = EoSB[:, imuB][:, 4]
        ind_s = len(s) - GetFirst(s[::-1]) - 1
        sH = s[:ind_s][-1] * T[:ind_s][-1] ** 3
        sHList.append(sH)

        muBListH.append(muB_val)
        TListH.append(TL(muB_val))

    try:
        NTB_A, NmuB_A = infer_grid_shape(EoSA, t_col=0, mub_col=1)

        EoSA = EoSA.reshape(NTB_A, NmuB_A, -1)
    except ValueError:
        print("Error inferring grid for EoSAboveGeV.dat:")
        raise

    eQList = []
    nBQList = []
    # new: pressure and entropy at Q side (above)
    pQList = []
    sQList = []

    for imuB in range(NmuB_A):
        T = EoSA[:, imuB][:, 0]
        eps = EoSA[:, imuB][:, 2]

        ind = GetFirst(eps) - 1
        eQ = eps[ind:][0] * T[ind:][0] ** 4
        eQList.append(eQ)

        nB = EoSA[:, imuB][:, 3]
        ind = GetFirst(nB) - 1
        nBQ = nB[ind:][0] * T[ind:][0] ** 3
        nBQList.append(nBQ)

        # pressure p/T^4 -> p [GeV^4]
        p = EoSA[:, imuB][:, 5]
        ind_p = GetFirst(p) - 1
        pQ = p[ind_p:][0] * T[ind_p:][0] ** 4
        pQList.append(pQ)

        # entropy s/T^3 -> s [GeV^3]
        s = EoSA[:, imuB][:, 4]
        ind_s = GetFirst(s) - 1
        sQ = s[ind_s:][0] * T[ind_s:][0] ** 3
        sQList.append(sQ)

    # arrays

    eHList = np.array(eHList)
    nBHList = np.array(nBHList)
    eQList = np.array(eQList)
    nBQList = np.array(nBQList)
    # arrays for new quantities
    pHList = np.array(pHList)
    sHList = np.array(sHList)
    pQList = np.array(pQList)
    sQList = np.array(sQList)

    with open("RegionS.dat", "w") as f:
        print(
            "#T [GeV] muB [GeV] eH [GeV/fm3] nBH [1/fm3] eQ [GeV/fm3] nBQ [1/fm3] pH [GeV/fm3] sH [1/fm3] pQ [GeV/fm3] sQ [1/fm3]",
            file=f,
        )
        for T, muB, eH, nBH, eQ, nBQ, pH, sH, pQ, sQ in zip(
            TListH,
            muBListH,
            eHList,
            nBHList,
            eQList,
            nBQList,
            pHList,
            sHList,
            pQList,
            sQList,
        ):
            f.write(
                f"{T:5} {muB:5} {eH:5} {nBH:5} {eQ:5} {nBQ:5} {pH:5} {sH:5} {pQ:5} {sQ:5}\n"
            )

    print("Data saved to RegionS.dat")


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except Exception as e:
        print(f"Error: {e}")
        print("Invalid arguments.")

    sys.exit()
