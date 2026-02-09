import os
import sys

import numpy as np
from matplotlib import pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from utils import read_parameters

OUTDIR = "plots_fin"
os.makedirs(OUTDIR, exist_ok=True)


def set_panel_border(ax):
    """Set thick border around each panel."""
    for spine in ax.spines.values():
        spine.set_linewidth(2.5)


try:
    param_path = sys.argv[1]
except FileNotFoundError:
    print("Parameter file is not found or specified.")
    sys.exit()

Param = read_parameters(param_path)
data = os.path.join(os.path.dirname(param_path), Param["OutputMergedEoS"])


# Columns in invMergedEoS_test.dat:
# Ttilde(GeV) muBtilde(GeV) e(GeV^4) nB(GeV^3) T(GeV) muB(GeV)
# P(GeV^4) s(GeV^3) cs2 chi2(GeV^2) chi

EoS = np.loadtxt(data, skiprows=1)

# Ensure that the number of rows matches the expected NTildeT * NTildemuB
NTildeT = Param["NTildeT"]
NTildemuB = Param["NTildemuB"]

nrows, ncols = EoS.shape
expected_rows = NTildeT * NTildemuB

if nrows != expected_rows:
    print(
        f"Warning: EoS rows ({nrows}) != expected grid ({expected_rows} = "
        f"{NTildeT}x{NTildemuB}). Adjusting data to match grid."
    )
    if nrows > expected_rows:
        # Too many rows: truncate extras
        EoS = EoS[:expected_rows]
    else:
        # Too few rows: pad with NaNs so that reshape works
        pad_rows = expected_rows - nrows
        pad = np.full((pad_rows, ncols), np.nan)
        EoS = np.vstack((EoS, pad))
    nrows, ncols = EoS.shape

# Extract columns (flat arrays)
Ttilde = EoS[:, 0]
mubtilde = EoS[:, 1]
e = EoS[:, 2]
nb = EoS[:, 3]
T = EoS[:, 4]
muB = EoS[:, 5]
P = EoS[:, 6]
S = EoS[:, 7]
chi2_flat = EoS[:, 9]
chi = EoS[:, 10]

plt.figure(figsize=(12, 5))
sc = plt.scatter(muB, T, c=chi, cmap="plasma", s=2)
cbar = plt.colorbar(sc, label=r"$\chi$ (fraction)")
cbar.set_label(r"$\chi$ (fraction)", fontsize=21)
plt.xlabel(r"$\mu_B$ [GeV]", fontsize=21)
plt.ylabel(r"T[GeV]", fontsize=21)
cbar.ax.tick_params(axis="both", which="major", labelsize=22)
plt.tick_params(axis="both", which="major", labelsize=18)
plt.grid(True)
plt.tight_layout()
set_panel_border(plt.gca())
plt.savefig(os.path.join(OUTDIR, "muB_vs_T_colored_by_chi.png"), dpi=300)
plt.show()


plt.figure(figsize=(12, 5))
sc = plt.scatter(nb, e, c=chi, cmap="plasma", s=2)
cbar = plt.colorbar(sc, label=r"$\chi$ (fraction)")
cbar.set_label(r"$\chi$ (fraction)", fontsize=21)
plt.xlabel(r"$n_B$ [GeV$^3$]", fontsize=21)
plt.ylabel(r"$e$ [GeV$^4$]", fontsize=21)
cbar.ax.tick_params(axis="both", which="major", labelsize=22)
plt.tick_params(axis="both", which="major", labelsize=18)
plt.grid(True)
plt.tight_layout()
set_panel_border(plt.gca())
plt.savefig(os.path.join(OUTDIR, "nB_vs_e_colored_by_chi.png"), dpi=300)
plt.show()


plt.figure(figsize=(10, 6))
sc = plt.scatter(mubtilde, Ttilde, c=chi, cmap="plasma", s=2)
cbar = plt.colorbar(sc, label=r"$\chi$ (fraction)")
cbar.set_label(r"$\chi$ (fraction)", fontsize=21)
plt.xlabel(r"$\tilde{\mu}_B$ [GeV]", fontsize=21)
plt.ylabel(r"$\tilde{T}$ [GeV]", fontsize=21)
plt.tick_params(axis="both", which="major", labelsize=18)
cbar.ax.tick_params(axis="both", which="major", labelsize=22)
plt.grid(True)
plt.tight_layout()
set_panel_border(plt.gca())
plt.savefig(os.path.join(OUTDIR, "muBtilde_vs_Ttilde_colored_by_chi.png"), dpi=300)
plt.show()


# Plots analogous to HydroIsingEos/plotEoS.py

# Reshape directly using NTildeT and NTildemuB from the parameter file.
nrows, ncols = EoS.shape
dat = EoS.reshape(NTildeT, NTildemuB, ncols)


fig, ax = plt.subplots(1, 1, figsize=(10, 8))
for i in range(NTildemuB):
    arr = dat[:, i]
    T_grid, e_grid = arr[:, 4], arr[:, 2]
    ax.plot(e_grid / T_grid**4, T_grid, lw=1)
ax.set_ylabel(r"$T$ (GeV)", fontsize=21)
ax.set_xlabel(r"$e/T^4$", fontsize=21)
ax.grid()
set_panel_border(ax)
plt.tight_layout()
plt.tick_params(axis="both", which="major", labelsize=18)
plt.savefig(os.path.join(OUTDIR, "Teps.png"), format="png", dpi=300)
plt.show()


fig, ax = plt.subplots(1, 1, figsize=(10, 8))
for i in range(NTildeT):  # To see the 1st order line under the folding.
    arr = dat[i, :]
    mub_grid, nb_grid, T_grid = arr[:, 5], arr[:, 3], arr[:, 4]
    ax.plot(nb_grid / T_grid**3, mub_grid, lw=1)

ax.set_xlabel(r"$n_B/T^3$", fontsize=21)
ax.set_ylabel(r"$\mu_B$ [GeV]", fontsize=21)
ax.set_xlim(-0.05, 1.64)
ax.grid()
set_panel_border(ax)
plt.tight_layout()
plt.tick_params(axis="both", which="major", labelsize=18)
plt.savefig(os.path.join(OUTDIR, "mubnb.png"), format="png", dpi=300)
plt.show()


fig, ax = plt.subplots(1, 1, figsize=(10, 8))
for i in range(NTildeT):
    arr = dat[i, :]
    chi2_grid, nb_grid, T_grid = arr[:, 9], arr[:, 3], arr[:, 4]
    ax.plot(nb_grid / T_grid**3, chi2_grid / T_grid**2, lw=1)

ax.set_xlabel(r"$n_B/T^3$", fontsize=21)
ax.set_ylabel(r"$\chi_2/T^2$", fontsize=21)
ax.set_xlim(-0.05, 1.37)
ax.grid()
ax.set_ylim(-0.01, 0.9)
set_panel_border(ax)
plt.tight_layout()
plt.tick_params(axis="both", which="major", labelsize=18)
plt.savefig(os.path.join(OUTDIR, "chinb.png"), format="png", dpi=300)
plt.show()
