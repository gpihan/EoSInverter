import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from utils import read_parameters


def set_panel_border(ax):
    for spine in ax.spines.values():
        spine.set_linewidth(2.5)


def create_interpolator(points, values):
    return LinearNDInterpolator(points, values, fill_value=np.nan)


def readTable(path):
    df = pd.read_csv(path, sep=r"\s+", comment="#", header=None)
    return df.to_numpy()


try:
    param_path = sys.argv[1]
except FileNotFoundError:
    print("Parameter file is not found or specified.")
    sys.exit()

# ============================================================
# Load inverse EoS
# columns (HydroIsingEos/EoS.dat):
# Ttilde(GeV) muBtilde(GeV) e(GeV^4) nB(GeV^3) T(GeV) muB(GeV) P(GeV^4) s(GeV^3) cs2 chi2(GeV^2) chi
# ============================================================

Param = read_parameters(param_path)
inv = os.path.join(os.path.dirname(param_path), Param["OutputMergedEoS"])

inv = np.loadtxt(inv, skiprows=2)


# ============================================================
# Load forward EoS
# columns (EoS2DIsing*.dat or similar):
# T(GeV) muB(GeV) e/T^4 nB/T^3 P/T^4 S/T^3 Cs2 chi2/T^2 chi3/T chi4
# Only the first four columns are used here.
# ============================================================

inp_raw = os.path.join(os.path.dirname(param_path), Param["EoS_table"])
inp_raw = np.loadtxt(inp_raw)

T_GeV = inp_raw[:, 0]
muB_GeV = inp_raw[:, 1]
e_T4 = inp_raw[:, 2]
nB_T3 = inp_raw[:, 3]

inp = np.column_stack([T_GeV, muB_GeV, e_T4, nB_T3])


def build_interpolators(inv_table, inp_table):
    # inverse: (e, nB) -> (T, muB)
    T_from_e_nb = create_interpolator(inv_table[:, 2:4], inv_table[:, 4])
    muB_from_e_nb = create_interpolator(inv_table[:, 2:4], inv_table[:, 5])

    # forward: (T, muB) -> (e/T^4, nB/T^3)
    eT4_from_T_muB = create_interpolator(inp_table[:, 0:2], inp_table[:, 2])
    nBT3_from_T_muB = create_interpolator(inp_table[:, 0:2], inp_table[:, 3])

    return T_from_e_nb, muB_from_e_nb, eT4_from_T_muB, nBT3_from_T_muB


def roundtrip_for_muB(muB0, T_grid, eT4_fwd, nBT3_fwd, T_inv, muB_inv):
    pts = np.column_stack([T_grid, np.full_like(T_grid, muB0)])

    eT4 = eT4_fwd(pts)
    nBT3 = nBT3_fwd(pts)

    mask = np.isfinite(eT4) & np.isfinite(nBT3)
    if not np.any(mask):
        return np.array([]), np.array([]), np.array([]), np.array([])

    T_in = T_grid[mask]
    muB_in = np.full_like(T_in, muB0)

    e = eT4[mask] * T_in**4
    nB = nBT3[mask] * T_in**3

    pts_enb = np.column_stack([e, nB])

    T_rt = T_inv(pts_enb)
    muB_rt = muB_inv(pts_enb)

    mask_rt = np.isfinite(T_rt) & np.isfinite(muB_rt)

    return T_in[mask_rt], muB_in[mask_rt], T_rt[mask_rt], muB_rt[mask_rt]


def plot_roundtrip():
    T_inv, muB_inv, eT4_fwd, nBT3_fwd = build_interpolators(inv, inp)
    TrLine = os.path.join(os.path.dirname(param_path), Param["TransitionLine"])
    tr = np.loadtxt(TrLine)

    Tmin, Tmax = np.nanmin(inp[:, 0]), np.nanmax(inp[:, 0])
    T_grid = np.linspace(Tmin, Tmax, 300)

    muB_list = [
        0.0,
        0.05,
        0.1,
        0.15,
        0.2,
        0.25,
        0.3,
        0.35,
        0.4,
        0.45,
        0.5,
        0.55,
        0.6,
        0.65,
    ]

    plt.figure(figsize=(8, 6))

    plt.plot(
        tr[:, 0] / 1000, tr[:, 1] / 1000, color="black", lw=2, label="Transition line"
    )

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    lines = []
    labels = []
    for i, muB0 in enumerate(muB_list):
        T_in, muB_in, T_rt, muB_rt = roundtrip_for_muB(
            muB0, T_grid, eT4_fwd, nBT3_fwd, T_inv, muB_inv
        )
        if T_in.size == 0:
            continue

        c = colors[i % len(colors)]
        (l1,) = plt.plot(muB_in, T_in, "-", color=c, lw=1)
        (l2,) = plt.plot(muB_rt, T_rt, "--", color=c, lw=1.8)
        lines.append(l1)
        labels.append(f"$\\mu_B = {muB0:.2f}$ GeV")

    plt.xlabel("μB [GeV]", fontsize=22)
    plt.ylabel("T [GeV]", fontsize=22)
    plt.tick_params(labelsize=18)
    plt.grid(alpha=0.3)
    set_panel_border(plt.gca())
    plt.legend(
        lines,
        labels,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=10,
    )
    plt.tight_layout(rect=[0, 0, 0.82, 1])
    plt.show()


def roundtrip_for_T(T0, muB_grid, eT4_fwd, nBT3_fwd, T_inv, muB_inv):
    pts = np.column_stack([np.full_like(muB_grid, T0), muB_grid])

    eT4 = eT4_fwd(pts)
    nBT3 = nBT3_fwd(pts)

    mask = np.isfinite(eT4) & np.isfinite(nBT3)
    if not np.any(mask):
        return np.array([]), np.array([]), np.array([]), np.array([])

    muB_in = muB_grid[mask]
    T_in = np.full_like(muB_in, T0)

    e = eT4[mask] * T_in**4
    nB = nBT3[mask] * T_in**3

    pts_enb = np.column_stack([e, nB])

    T_rt = T_inv(pts_enb)
    muB_rt = muB_inv(pts_enb)

    mask_rt = np.isfinite(T_rt) & np.isfinite(muB_rt)

    return muB_in[mask_rt], T_in[mask_rt], muB_rt[mask_rt], T_rt[mask_rt]


def plot_roundtrip_fixed_T():
    T_inv, muB_inv, eT4_fwd, nBT3_fwd = build_interpolators(inv, inp)
    TrLine = os.path.join(os.path.dirname(param_path), Param["TransitionLine"])
    tr = np.loadtxt(TrLine)

    muBmin, muBmax = np.nanmin(inp[:, 1]), np.nanmax(inp[:, 1])
    muB_grid = np.linspace(muBmin, muBmax, 300)

    T_list = [0.05, 0.075, 0.10, 0.12, 0.13, 0.14, 0.15, 0.20, 0.25, 0.30, 0.35]

    plt.figure(figsize=(8, 6))

    plt.plot(tr[:, 0] / 1000, tr[:, 1] / 1000, color="black", lw=2)

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    lines = []
    labels = []
    for i, T0 in enumerate(T_list):
        muB_in, T_in, muB_rt, T_rt = roundtrip_for_T(
            T0, muB_grid, eT4_fwd, nBT3_fwd, T_inv, muB_inv
        )
        if muB_in.size == 0:
            continue

        c = colors[i % len(colors)]
        (l1,) = plt.plot(muB_in, T_in, "-", color=c, lw=1)
        (l2,) = plt.plot(muB_rt, T_rt, "--", color=c, lw=1.8)
        lines.append(l1)
        labels.append(f"$T = {T0:.3f}$ GeV")

    plt.xlabel("μB [GeV]", fontsize=22)
    plt.ylabel("T [GeV]", fontsize=22)
    plt.tick_params(labelsize=18)
    plt.grid(alpha=0.3)
    set_panel_border(plt.gca())
    plt.legend(
        lines,
        labels,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=10,
    )
    plt.tight_layout(rect=[0, 0, 0.82, 1])
    plt.show()


if __name__ == "__main__":
    plot_roundtrip()
    plot_roundtrip_fixed_T()
