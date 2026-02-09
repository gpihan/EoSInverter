import argparse
import os

import numpy as np

from utils import readTable


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Preprocess EoS tables for inversion.")
    p.add_argument(
        "--table-input",
        default="Ising-TExS_output_thermodynamics.csv",
        help="input table",
    )
    p.add_argument(
        "--table-output",
        default="Ising-TExS_output_thermodynamics.dat",
        help="output table",
    )
    return p.parse_args(argv)


def main():
    args = parse_args()

    base_dir = os.path.dirname(os.path.abspath(__file__))
    in_path = os.path.join(base_dir, args.table_input)
    out_path = os.path.join(base_dir, args.table_output)

    # Use appropriate reader depending on file format
    if in_path.endswith(".csv"):
        data = np.loadtxt(in_path, delimiter=",", comments="#")
    else:
        data = readTable(in_path)

    # Expected columns mapping
    T_mev = data[:, 0]
    muB_mev = data[:, 1]
    B_dens = data[:, 2]
    E_dens = data[:, 3]
    P = data[:, 4]
    c_s = data[:, 5]  # Not used
    chi2B = data[:, 6]  # Not used
    s_dens = data[:, 7]

    T_gev = T_mev / 1000.0
    muB_gev = muB_mev / 1000.0

    T3 = T_mev**3
    T4 = T_mev**4

    e_T4 = E_dens / T4
    nB_T3 = B_dens / T3
    P_T4 = P / T4
    s_T3 = s_dens / T3

    sort_idx = np.lexsort((muB_gev, T_gev))
    output = np.column_stack([T_gev, muB_gev, e_T4, nB_T3, P_T4, s_T3, c_s, chi2B])[
        sort_idx
    ]
    header = "#T[GeV] muB[GeV] e/T^4 nB/T^3 P/T^4 s/T^3 c_s chi2B"

    np.savetxt(out_path, output, fmt="%.8e", header=header, comments="")
    print(f"Saved to '{args.table_output}'")


if __name__ == "__main__":
    main()
