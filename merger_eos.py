import os
import sys

import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator, interp1d
from utils import read_parameters, readTable


def ToEN(Tt, mbt):
    e = 19 * np.pi**2 / 12 * Tt**4
    nb = 1 / 3 * (mbt * Tt**2)
    return e, nb


def Get2DTilde(e, nb):
    Ttilde = (12 / (19 * np.pi**2) * e) ** 0.25
    muBtilde = 3 * (nb / Ttilde**2)
    return Ttilde, muBtilde


def create_interpolator(points, values):
    lin = LinearNDInterpolator(points, values)
    near = NearestNDInterpolator(points, values)

    def interp(x, y):
        val = lin(x, y)

        if val is None or not np.isfinite(val):
            val = near(x, y)

        return val

    return interp


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


def find_closest_transition(e, nB, table, seg_norm_dist_max=10, blend_width=0.05):
    """Find closest transition and interpolate thermodynamic quantities.

    Parameters
    ----------
    e : float
        Energy density
    nB : float
        Baryon density
    table : array-like
        Transition line data (T, muB, eH, nBH, eQ, nBQ, [pH, sH, pQ, sQ])
    seg_norm_dist_max : float
        Maximum normalized distance to transition line
    blend_width : float
        Width of blending region for smooth transition (fraction of seg_norm_dist_max)

    Returns
    -------
    tuple : (T, muB, P, S, chi, chi_e, chi_nB, blend_factor)
        blend_factor: 1.0 = full RegionS, 0.0 = use interpolators
    """
    # First, try quadrilateral containment to directly assign muB

    # Fallback: previous closest-segment logic with smooth blending
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
        return False, False, False, False, False, None, None, 0.0

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

    # Calculate blend factor for smooth transition
    # blend_factor = 1.0 at norm_d=0, smoothly goes to 0 at norm_d=seg_norm_dist_max
    if seg_norm_dist_max is not None and seg_norm_dist_max > 0:
        # Smooth blending using cosine interpolation
        if norm_d_min <= seg_norm_dist_max * (1 - blend_width):
            blend_factor = 1.0
        elif norm_d_min >= seg_norm_dist_max:
            blend_factor = 0.0
        else:
            # Smooth transition in the blend region
            x = (norm_d_min - seg_norm_dist_max * (1 - blend_width)) / (
                seg_norm_dist_max * blend_width
            )
            blend_factor = 0.5 * (1.0 + np.cos(np.pi * x))
    else:
        blend_factor = 1.0 if norm_d_min < float("inf") else 0.0

    # Apply distance thresholds (absolute and normalized)
    if (seg_norm_dist_max is not None) and (norm_d_min > seg_norm_dist_max):
        return False, False, False, False, False, None, None, 0.0

    chi = float(np.clip(t_min, 0.0, 1.0))  # Use actual projection parameter
    chi_e = None
    chi_nB = None

    # Interpolate T and muB along the segment
    T = Tc_f
    muB = muBc_f

    # Interpolate pressure and entropy along the segment if available
    P = S = None
    if chi is not None and (pH_f is not None) and (pQ_f is not None):
        P = (1 - chi) * pH_f + chi * pQ_f
    if chi is not None and (sH_f is not None) and (sQ_f is not None):
        S = (1 - chi) * sH_f + chi * sQ_f

    return T, muB, P, S, chi, chi_e, chi_nB, blend_factor


def run_merger(
    EoS,
    TrLine,
    datS,
    out_path,
    Ne,
    Nn,
    show_progress,
    seg_norm_dist_max=0.15,
    blend_width=0.05,
):
    """Merge EoS regions with smooth blending.

    Parameters
    ----------
    EoS : str
        Path to EoS table file
    TrLine : str
        Path to transition line file
    datS : str
        Path to RegionS data file
    out_path : str
        Output file path
    Ne : int
        Number of energy density points
    Nn : int
        Number of baryon density points
    show_progress : bool
        Whether to show progress
    seg_norm_dist_max : float
        Maximum normalized distance to RegionS
    blend_width : float
        Width of smooth blending region (fraction of seg_norm_dist_max)
    """
    EoS = readTable(EoS)
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

    # EoS columns (assumed):
    # 0: e, 1: nB, 2: T, 3: muB, 4: P, 5: s, 6: cs2, 7: chi2
    TEoS = create_interpolator(EoS[:, :2], EoS[:, 2])
    muBEoS = create_interpolator(EoS[:, :2], EoS[:, 3])
    PEoS = create_interpolator(EoS[:, :2], EoS[:, 4]) if EoS.shape[1] >= 6 else None
    SEoS = create_interpolator(EoS[:, :2], EoS[:, 5]) if EoS.shape[1] >= 6 else None
    cs2EoS = create_interpolator(EoS[:, :2], EoS[:, 6]) if EoS.shape[1] >= 7 else None
    chi2EoS = create_interpolator(EoS[:, :2], EoS[:, 7]) if EoS.shape[1] >= 8 else None

    T_min = (EoS[:, 2]).min()
    T_max = (EoS[:, 2]).max()

    mub_min = (EoS[:, 3]).min()
    mub_max = (EoS[:, 3]).max()
    print(T_min, T_max, mub_min, mub_max)

    TTilEoS, muBBtilEoS = Get2DTilde(EoS[:, 0], EoS[:, 1])

    TTil_min = float(TTilEoS.min())
    TTil_max = float(TTilEoS.max())
    muBBtil_min = float(muBBtilEoS.min())
    muBBtil_max = float(muBBtilEoS.max())

    TTILDEArr = np.linspace(TTil_min, TTil_max, Ne)
    mubTILDEArr = np.linspace(muBBtil_min, muBBtil_max, Nn)

    print(f"Writing merged EoS to {out_path} ...")
    print(
        f"Blend width: {blend_width * seg_norm_dist_max:.3f} (relative: {blend_width})"
    )

    with open(out_path, "w") as f:
        # nový header
        f.write(
            "# Ttilde(GeV) muBtilde(GeV) e(GeV^4) nB(GeV^3) "
            "T(GeV) muB(GeV) P(GeV^4) s(GeV^3) cs2 chi2(GeV^2) chi\n"
        )

        for i_T, Ttilde in enumerate(TTILDEArr):
            if show_progress and (i_T % max(1, Ne // 10) == 0):
                print(f"Progress: {i_T}/{Ne} rows ({100.0 * i_T / max(1, Ne):.1f}%)")
            for i_mub, muBtilde in enumerate(mubTILDEArr):
                e, nB = ToEN(Ttilde, muBtilde)
                res = find_closest_transition(
                    e,
                    nB,
                    datS,
                    seg_norm_dist_max=seg_norm_dist_max,
                    blend_width=blend_width,
                )

                if res[0] is not False:
                    (
                        T_regionS,
                        muB_regionS,
                        P_regionS,
                        S_regionS,
                        chi,
                        chi_e,
                        chi_n,
                        blend_factor,
                    ) = res

                    # Get values from interpolators
                    T_interp = TEoS(e, nB)
                    muB_interp = muBEoS(e, nB)
                    P_interp = PEoS(e, nB) if PEoS is not None else None
                    S_interp = SEoS(e, nB) if SEoS is not None else None
                    cs2_interp = cs2EoS(e, nB) if cs2EoS is not None else None
                    chi2_interp = chi2EoS(e, nB) if chi2EoS is not None else None

                    # Blend between RegionS and interpolators
                    T = blend_factor * T_regionS + (1 - blend_factor) * T_interp
                    muB = blend_factor * muB_regionS + (1 - blend_factor) * muB_interp

                    # Blend P and S if available
                    P = None
                    S = None
                    if P_regionS is not None and P_interp is not None:
                        P = blend_factor * P_regionS + (1 - blend_factor) * P_interp
                    elif P_interp is not None:
                        P = P_interp
                    elif P_regionS is not None:
                        P = P_regionS

                    if S_regionS is not None and S_interp is not None:
                        S = blend_factor * S_regionS + (1 - blend_factor) * S_interp
                    elif S_interp is not None:
                        S = S_interp
                    elif S_regionS is not None:
                        S = S_regionS

                    # cs2, chi2 bereme z EoS interpolátorů
                    cs2 = cs2_interp
                    chi2 = chi2_interp

                else:
                    # fallback: use interpolators only
                    T = TEoS(e, nB)
                    muB = muBEoS(e, nB)
                    P = PEoS(e, nB) if PEoS is not None else None
                    S = SEoS(e, nB) if SEoS is not None else None
                    cs2 = cs2EoS(e, nB) if cs2EoS is not None else None
                    chi2 = chi2EoS(e, nB) if chi2EoS is not None else None

                    if muB < 0.400:
                        chi = -1.0
                    elif muB >= 0.400 and T < TrLine(muB):
                        chi = 0.0
                    elif muB >= 0.400 and T > TrLine(muB):
                        chi = 1.0

                f.write(
                    f" {Ttilde:.8e} {muBtilde:.8e} {e:.8e} {nB:.8e} "
                    f"{T:.8e} {muB:.8e} {P:.8e} {S:.8e} {cs2:.8e} {chi2:.8e} {chi:.8e}\n"
                )

    print("... done.")
    print(f"Output written to: {out_path}")


def main(argv):
    try:
        param_path = sys.argv[1]
    except FileNotFoundError:
        print("Parameter file is not found or specified.")
        sys.exit()

    Param = read_parameters(param_path)

    EoS = Param["OutputFolder"] + "/EoS_all.dat"

    TrLine = os.path.join(Param["TransitionLine"])
    datS = os.path.join(Param["RegionS"])
    out_path = Param["OutputMergedEoS"]
    Ne = Param["NTildeT"]
    Nn = Param["NTildemuB"]

    print("Running EoS merger...")
    print("Grid size: ", Ne * Nn)

    run_merger(
        EoS,
        TrLine,
        datS,
        out_path,
        Ne,
        Nn,
        show_progress=True,
        seg_norm_dist_max=0.15,
        blend_width=0.05,  # 5% blending region for smooth transition
    )


if __name__ == "__main__":
    try:
        status = main(sys.argv[1:])
    except Exception as e:
        print(f"Error: {e}")
        print("Invalid arguments.")

    sys.exit(0)
