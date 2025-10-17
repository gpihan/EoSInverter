# Copyright Tomáš Poledníček, Gregoire Pihan @ 2025
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from utils import readTable


def filter_eos(eos_df, transition_table, keep_above):
    eos_df = pd.DataFrame(eos_df)

    mub_arr = transition_table[:, 0] / 1000.0
    T_arr = transition_table[:, 1] / 1000.0
    sorter = np.argsort(mub_arr)

    mub_sorted = mub_arr[sorter]
    T_sorted = T_arr[sorter]

    T_transition_interp = interp1d(
        mub_sorted, T_sorted, bounds_error=False, fill_value="extrapolate"
    )

    muB_vals = eos_df.iloc[:, 3].values
    T_vals = eos_df.iloc[:, 2].values
    T_transition_vals = T_transition_interp(muB_vals)

    if keep_above:
        mask = T_vals > T_transition_vals
    else:
        mask = T_vals <= T_transition_vals

    eos_clean = eos_df[mask].reset_index(drop=True)
    return eos_clean.to_numpy()


def plot(
    eoscross,
    eos_below_raw,
    eos_above_raw,
    eos_below_clean,
    eos_above_clean,
    transition_line,
):
    fig, ax = plt.subplots(1, 2, figsize=(10, 8))
    ax[0].scatter(
        eoscross[:, 3],
        eoscross[:, 2],
        color="green",
        alpha=0.5,
        label="EoS cross",
    )
    ax[0].scatter(
        eos_below_raw[:, 3],
        eos_below_raw[:, 2],
        color="blue",
        alpha=0.4,
        label="EoS below raw",
    )
    ax[0].scatter(
        eos_above_raw[:, 3],
        eos_above_raw[:, 2],
        color="red",
        alpha=0.6,
        label="EoS above raw",
    )
    ax[0].plot(
        transition_line[:, 0] / 1000,
        transition_line[:, 1] / 1000,
        color="black",
        linewidth=2,
        label="Transition Line",
    )
    ax[0].set_xlabel(r"$\mu_B$ [GeV]")
    ax[0].set_ylabel(r"$T$ [GeV]")
    ax[0].grid(True)

    ax[1].scatter(
        eoscross[:, 3],
        eoscross[:, 2],
        color="green",
        alpha=0.5,
        label="EoS cross",
    )

    ax[1].scatter(
        eos_below_clean[:, 3],
        eos_below_clean[:, 2],
        color="blue",
        alpha=0.4,
        label="EoS below clean",
    )
    ax[1].scatter(
        eos_above_clean[:, 3],
        eos_above_clean[:, 2],
        color="red",
        alpha=0.6,
        label="EoS above clean",
    )
    ax[1].plot(
        transition_line[:, 0] / 1000,
        transition_line[:, 1] / 1000,
        color="black",
        linewidth=2,
        label="Transition Line",
    )
    ax[1].set_xlabel(r"$\mu_B$ [GeV]")
    ax[1].set_ylabel(r"$T$ [GeV]")
    ax[1].legend()
    ax[1].grid(True)

    plt.tight_layout()

    plt.show()


def main():
    base_dir = os.getcwd()
    output_dir = os.path.join(base_dir, "output/")
    os.makedirs(output_dir, exist_ok=True)

    transition_line = readTable(base_dir + "/input/TransitionLine.dat")

    eos_above_raw = readTable(base_dir + "/EoSAboveGeV/EoS_all.dat")
    eos_below_raw = readTable(base_dir + "/EoSBelowGeV/EoS_all.dat")
    eos_cross = readTable(base_dir + "/EoSCrossGeV/EoS_all.dat")

    eos_above_clean = filter_eos(eos_above_raw, transition_line, keep_above=True)
    eos_below_clean = filter_eos(eos_below_raw, transition_line, keep_above=False)

    pd.DataFrame(eos_above_clean).to_csv(
        os.path.join(output_dir, "EoSAbove_clean.dat"),
        index=False,
        header=False,
        sep="\t",
    )

    pd.DataFrame(eos_below_clean).to_csv(
        os.path.join(output_dir, "EoSBelow_clean.dat"),
        index=False,
        header=False,
        sep="\t",
    )
    pd.DataFrame(eos_cross).to_csv(
        os.path.join(output_dir, "EoSCross_clean.dat"),
        index=False,
        header=False,
        sep="\t",
    )

    plot(
        eos_cross,
        eos_below_raw,
        eos_above_raw,
        eos_below_clean,
        eos_above_clean,
        transition_line,
    )

    print(f"Number of points in EoSAbove (clean): {len(eos_above_clean)}")
    print(f"Number of points in EoSBelow (clean): {len(eos_below_clean)}")
    print(f"Files saved to {output_dir}")


if __name__ == "__main__":
    main()
