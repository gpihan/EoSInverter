import os
import sys

import numpy as np
from scipy.interpolate import interp1d
from utils import read_parameters, readTable


def project_eos_columns(arr, e_i, nb_i, T_i, muB_i):
    ncols = arr.shape[1] if arr.ndim == 2 else 0
    idxs = [e_i, nb_i, T_i, muB_i]
    if any(i < 0 or i >= ncols for i in idxs):
        raise ValueError(
            f"Column index out of range for array with {ncols} columns: {idxs}"
        )
    return arr[:, idxs]


def Get2DTilde(e, nb):
    Ttilde = (12 / (19 * np.pi**2) * e) ** 0.25
    muBtilde = 5 * nb / Ttilde**2
    return Ttilde, muBtilde


def above_line_mask(
    T,
    muB,
    TrLine_interpolator,
    tol=0.0,
):
    """Vectorized mask: point is above the transition line and muB in domain."""
    Tc = TrLine_interpolator(muB)
    tc_ok = np.isfinite(Tc)
    return tc_ok & (T > (Tc + tol))


def below_line_mask(
    T,
    muB,
    TrLine_interpolator,
    tol=0.0,
):
    """Vectorized mask: point is below the transition line and muB in domain."""
    Tc = TrLine_interpolator(muB)
    tc_ok = np.isfinite(Tc)
    return tc_ok & (T < (Tc - tol))


def check_boundary_overlap(
    eos_above_raw,
    eos_cross_raw,
    *,
    e_i=0,
    nb_i=1,
    T_i=2,
    muB_i=3,
    nb_tolerance=1e-5,
    T_tolerance=1e-5,
    muB_atol=1e-9,
):
    """Find boundary overlap between above and cross for minimum/maximum muB.

    Returns a set of keys (e, nb, T) of points from ABOVE that should be removed
    due to overlap with CROSS at boundary slices.
    """
    if eos_above_raw.size == 0 or eos_cross_raw.size == 0:
        return set()

    mub_above_min = float(np.nanmin(eos_above_raw[:, muB_i]))
    mub_cross_max = float(np.nanmax(eos_cross_raw[:, muB_i]))

    # Tolerance-aware slicing
    above_mask = np.isclose(eos_above_raw[:, muB_i], mub_above_min, atol=muB_atol)
    cross_mask = np.isclose(eos_cross_raw[:, muB_i], mub_cross_max, atol=muB_atol)
    above_boundary = eos_above_raw[above_mask]
    cross_boundary = eos_cross_raw[cross_mask]

    if len(above_boundary) == 0 or len(cross_boundary) == 0:
        return set()

    # near zero along nb or T
    above_near_zero_nb = above_boundary[np.abs(above_boundary[:, nb_i]) < nb_tolerance]
    cross_near_zero_nb = cross_boundary[np.abs(cross_boundary[:, nb_i]) < nb_tolerance]
    above_near_zero_T = above_boundary[np.abs(above_boundary[:, T_i]) < T_tolerance]
    cross_near_zero_T = cross_boundary[np.abs(cross_boundary[:, T_i]) < T_tolerance]

    remove_points = set()

    # nb-based overlap
    if len(above_near_zero_nb) > 0 and len(cross_near_zero_nb) > 0:
        max_cross_nb = float(np.max(np.abs(cross_near_zero_nb[:, nb_i])))
        min_above_nb = float(np.min(np.abs(above_near_zero_nb[:, nb_i])))
        boundary_nb = 0.5 * (max_cross_nb + min_above_nb)
        # Mark above points with |nb| <= boundary_nb
        sel = np.abs(above_boundary[:, nb_i]) <= boundary_nb
        for row in above_boundary[sel]:
            e, nb, T = float(row[e_i]), float(row[nb_i]), float(row[T_i])
            remove_points.add((e, nb, T))

    # T-based overlap
    if len(above_near_zero_T) > 0 and len(cross_near_zero_T) > 0:
        max_cross_T = float(np.max(np.abs(cross_near_zero_T[:, T_i])))
        min_above_T = float(np.min(np.abs(above_near_zero_T[:, T_i])))
        boundary_T = 0.5 * (max_cross_T + min_above_T)
        sel = np.abs(above_boundary[:, T_i]) <= boundary_T
        for row in above_boundary[sel]:
            e, nb, T = float(row[e_i]), float(row[nb_i]), float(row[T_i])
            remove_points.add((e, nb, T))

    return remove_points


def remove_rows_by_keys(arr, remove_keys, e_i, nb_i, T_i):
    """Remove rows from numpy array based on keys (e, nb, T)."""
    if not remove_keys or len(arr) == 0:
        return arr

    mask = np.ones(len(arr), dtype=bool)
    for i, row in enumerate(arr):
        key = (float(row[e_i]), float(row[nb_i]), float(row[T_i]))
        if key in remove_keys:
            mask[i] = False

    return arr[mask]


def filter_array_by_mask(arr, mask):
    """Filter numpy array by mask."""
    return arr[mask] if len(arr) > 0 else arr


def pick_extreme_nb_by_e(arr, e_i, nb_i, pick_max=True):
    """For each e value, select the row with max/min nb."""
    if len(arr) == 0:
        return arr

    # Group by e values
    unique_e = np.unique(arr[:, e_i])
    selected_rows = []

    for e_val in unique_e:
        e_mask = arr[:, e_i] == e_val
        e_group = arr[e_mask]

        if len(e_group) > 0:
            if pick_max:
                idx = np.argmax(e_group[:, nb_i])
            else:
                idx = np.argmin(e_group[:, nb_i])
            selected_rows.append(e_group[idx])

    return np.array(selected_rows) if selected_rows else np.empty((0, arr.shape[1]))


def save_array_to_file(arr, file_path):
    """Save numpy array to file as TSV."""
    np.savetxt(file_path, arr, delimiter="\t", fmt="%.12g")


def main(argv):
    """Clean and merge EoS data with robust filters and parametrized paths."""
    script_dir = os.path.dirname(os.path.abspath(__file__))

    try:
        param_path = sys.argv[1]
    except FileNotFoundError:
        print("Parameter file is not found or specified.")
        sys.exit()

    Param = read_parameters(param_path)

    EoS_above = os.path.join(Param["EoS_above"])
    EoS_below = Param["EoS_below"]
    EoS_cross = Param["EoS_cross"]
    TrLine = Param["TransitionLine"]
    output = Param["premerger_output"]
    transition_path = os.path.join(script_dir, TrLine)
    above_path = os.path.join(script_dir, EoS_above)
    below_path = os.path.join(script_dir, EoS_below)
    cross_path = os.path.join(script_dir, EoS_cross)
    output_dir = os.path.join(script_dir, output)

    col_e = 0
    col_nb = 1
    col_T = 2
    col_muB = 3

    muB_atol = 1e-9
    nb_tolerance = 1e-4
    T_tolerance = 1e-4
    line_tol = 10e-5
    T_max = None

    os.makedirs(output_dir, exist_ok=True)

    # Transition line T(muB); input given in MeV -> convert to GeV
    datTRLine = readTable(transition_path)
    mu_vals = datTRLine[:, 0] / 1000.0
    T_vals = datTRLine[:, 1] / 1000.0

    # Create interpolator for transition line (do not extrapolate: return NaN outside domain)
    TrLine = interp1d(
        mu_vals,
        T_vals,
        kind="linear",
        bounds_error=False,
        fill_value=np.nan,
    )

    eos_above_raw_full = readTable(above_path)
    eos_below_raw_full = readTable(below_path)
    eos_cross_raw_full = readTable(cross_path)

    # Preserve ALL columns, but remember indices of e, nb, T, muB
    ncols_above = eos_above_raw_full.shape[1] if eos_above_raw_full.ndim == 2 else 0
    ncols_below = eos_below_raw_full.shape[1] if eos_below_raw_full.ndim == 2 else 0
    ncols_cross = eos_cross_raw_full.shape[1] if eos_cross_raw_full.ndim == 2 else 0

    # Overlap removal between above and cross at boundary
    overlapping_points = check_boundary_overlap(
        eos_above_raw_full,
        eos_cross_raw_full,
        nb_tolerance=nb_tolerance,
        T_tolerance=T_tolerance,
        muB_atol=muB_atol,
        e_i=col_e,
        nb_i=col_nb,
        T_i=col_T,
        muB_i=col_muB,
    )

    # Apply masks for above/below relative to line
    if len(eos_above_raw_full) > 0:
        mask_above = above_line_mask(
            eos_above_raw_full[:, col_T],
            eos_above_raw_full[:, col_muB],
            TrLine,
            tol=line_tol,
        )
        eos_above_filtered = filter_array_by_mask(eos_above_raw_full, mask_above)
    else:
        eos_above_filtered = eos_above_raw_full

    if len(eos_below_raw_full) > 0:
        mask_below = below_line_mask(
            eos_below_raw_full[:, col_T],
            eos_below_raw_full[:, col_muB],
            TrLine,
            tol=line_tol,
        )
        eos_below_filtered = filter_array_by_mask(eos_below_raw_full, mask_below)
    else:
        eos_below_filtered = eos_below_raw_full

    # Remove overlap points from ABOVE by matching (e, nb, T)
    eos_above_filtered = remove_rows_by_keys(
        eos_above_filtered, overlapping_points, col_e, col_nb, col_T
    )

    # Apply optional T_max cutoff before boundary slicing
    if T_max is not None:
        if len(eos_above_filtered) > 0:
            T_mask = eos_above_filtered[:, col_T] <= T_max
            eos_above_filtered = filter_array_by_mask(eos_above_filtered, T_mask)

        if len(eos_cross_raw_full) > 0:
            T_mask = eos_cross_raw_full[:, col_T] <= T_max
            eos_cross_filtered = filter_array_by_mask(eos_cross_raw_full, T_mask)
        else:
            eos_cross_filtered = eos_cross_raw_full
    else:
        eos_cross_filtered = eos_cross_raw_full

    # Removed specific filtering for BELOW - keeping only general line-based filtering

    # muB boundary values
    if len(eos_above_filtered) > 0:
        mub_above_min = float(np.nanmin(eos_above_filtered[:, col_muB]))
        mub_above_max = float(np.nanmax(eos_above_filtered[:, col_muB]))
    else:
        mub_above_min = np.nan
        mub_above_max = np.nan

    if len(eos_below_filtered) > 0:
        mub_below_min = float(np.nanmin(eos_below_filtered[:, col_muB]))
        mub_below_max = float(np.nanmax(eos_below_filtered[:, col_muB]))
    else:
        mub_below_min = np.nan
        mub_below_max = np.nan

    if len(eos_cross_filtered) > 0:
        mub_cross_min = float(np.nanmin(eos_cross_filtered[:, col_muB]))
        mub_cross_max = float(np.nanmax(eos_cross_filtered[:, col_muB]))

    # ABOVE: Split boundary slices
    if len(eos_above_filtered) > 0:
        above_mu_min_mask = np.isclose(
            eos_above_filtered[:, col_muB], mub_above_min, atol=muB_atol
        )
        above_mu_max_mask = np.isclose(
            eos_above_filtered[:, col_muB], mub_above_max, atol=muB_atol
        )

        above_mu_min_slice = filter_array_by_mask(eos_above_filtered, above_mu_min_mask)
        above_mu_max_slice = filter_array_by_mask(eos_above_filtered, above_mu_max_mask)
        above_core = filter_array_by_mask(
            eos_above_filtered, ~(above_mu_min_mask | above_mu_max_mask)
        )

        above_mu_min_picked = pick_extreme_nb_by_e(
            above_mu_min_slice, col_e, col_nb, pick_max=True
        )
        above_mu_max_picked = pick_extreme_nb_by_e(
            above_mu_max_slice, col_e, col_nb, pick_max=False
        )
    else:
        above_mu_min_picked = np.empty((0, ncols_above))
        above_mu_max_picked = np.empty((0, ncols_above))
        above_core = np.empty((0, ncols_above))

    # BELOW: Split boundary slices
    if len(eos_below_filtered) > 0:
        below_mu_min_mask = np.isclose(
            eos_below_filtered[:, col_muB], mub_below_min, atol=muB_atol
        )
        below_mu_max_mask = np.isclose(
            eos_below_filtered[:, col_muB], mub_below_max, atol=muB_atol
        )

        below_mu_min_slice = filter_array_by_mask(eos_below_filtered, below_mu_min_mask)
        below_mu_max_slice = filter_array_by_mask(eos_below_filtered, below_mu_max_mask)
        below_core = filter_array_by_mask(
            eos_below_filtered, ~(below_mu_min_mask | below_mu_max_mask)
        )

        below_mu_min_picked = pick_extreme_nb_by_e(
            below_mu_min_slice, col_e, col_nb, pick_max=True
        )
        below_mu_max_picked = pick_extreme_nb_by_e(
            below_mu_max_slice, col_e, col_nb, pick_max=False
        )
    else:
        below_mu_min_picked = np.empty((0, ncols_below))
        below_mu_max_picked = np.empty((0, ncols_below))
        below_core = np.empty((0, ncols_below))

    # CROSS: Only remove muB min slice and do deduplication for max slice
    if len(eos_cross_filtered) > 0:
        cross_mu_min_mask = np.isclose(
            eos_cross_filtered[:, col_muB], mub_cross_min, atol=muB_atol
        )
        cross_mu_max_mask = np.isclose(
            eos_cross_filtered[:, col_muB], mub_cross_max, atol=muB_atol
        )

        cross_mu_max_slice = filter_array_by_mask(eos_cross_filtered, cross_mu_max_mask)
        cross_core = filter_array_by_mask(
            eos_cross_filtered, ~(cross_mu_min_mask | cross_mu_max_mask)
        )

        # Deduplication: pick extreme nb for each e value at max muB boundary
        cross_mu_max_picked = pick_extreme_nb_by_e(
            cross_mu_max_slice, col_e, col_nb, pick_max=False
        )
    else:
        cross_mu_max_picked = np.empty((0, ncols_cross))
        cross_core = np.empty((0, ncols_cross))

    # Concat in the original order logic: below, above, cross
    eos_below_final = (
        np.vstack([below_mu_min_picked, below_core, below_mu_max_picked])
        if len(below_mu_min_picked) + len(below_core) + len(below_mu_max_picked) > 0
        else np.empty((0, ncols_below))
    )
    eos_above_final = (
        np.vstack([above_mu_min_picked, above_core, above_mu_max_picked])
        if len(above_mu_min_picked) + len(above_core) + len(above_mu_max_picked) > 0
        else np.empty((0, ncols_above))
    )
    eos_cross_final = (
        np.vstack([cross_core, cross_mu_max_picked])
        if len(cross_core) + len(cross_mu_max_picked) > 0
        else np.empty((0, ncols_cross))
    )

    # Sort by e (col 0)
    if len(eos_above_final) > 0:
        sort_indices = np.argsort(eos_above_final[:, col_e], kind="mergesort")
        eos_above_final = eos_above_final[sort_indices]

    if len(eos_below_final) > 0:
        sort_indices = np.argsort(eos_below_final[:, col_e], kind="mergesort")
        eos_below_final = eos_below_final[sort_indices]

    if len(eos_cross_final) > 0:
        sort_indices = np.argsort(eos_cross_final[:, col_e], kind="mergesort")
        eos_cross_final = eos_cross_final[sort_indices]

    # Write outputs
    save_array_to_file(eos_above_final, os.path.join(output_dir, "EoSAbove_clean.dat"))
    save_array_to_file(eos_below_final, os.path.join(output_dir, "EoSBelow_clean.dat"))
    save_array_to_file(eos_cross_final, os.path.join(output_dir, "EoSCross_clean.dat"))


if __name__ == "__main__":
    try:
        status = main(sys.argv[1:])
    except Exception as e:
        print(f"Error: {e}")
        print("Invalid arguments.")

    sys.exit()
