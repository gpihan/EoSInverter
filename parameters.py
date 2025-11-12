general_parameters = {
    "RunMode": 0,  # Run Mode: 0 -> Run locally
    # Run Mode: 1 -> Run on slurm cluster
    "Dimension": 2,  # Dimensions = 1 gives T(e)
    # Dimensions = 2 gives T(e, nb), mub(T, nb)
    # Dimensions = 3 gives T(e, nb, nq), mub(e, nb, nq), muq(e, nb, nq)
    # Dimensions = 4 gives T(e, nb, nq, ns), mub(e, nb, nq, ns), muq(e, nb, nq, ns), mus(e, nb, nq, ns)
    "AutoSetBoundaries": True,
    "NT": 200,
    "NB": 200,
    "NQ": 2,
    "NS": 2,
    "Ttilde": [0.009777504201390019, 0.31403115649771435, 10],
    "muBtilde": [-2.2732322158914076, 2.2732322158914076, 10],
    "muQtilde": [-0.6292193138512704, 0.6292193138512704, 10],
    "muStilde": [-0.6880440931493973, 0.6880440931493973, 10],
    "Accuracy": 1e-6,
    "MAXITER": 100,
    "Number_of_cores": 12,
    "hydro_model": "vHLLE",  # MUSIC
    "EoS_table": "EoSCrossGeV.dat",
    "OutputFolder": "EoSCrossGeV_test",
    # Paths to EoS files for merger
    "EoS_above_input": "EoSAboveGeV.dat",  # input table for inversion above transition line
    "EoS_below_input": "EoSBelowGeV.dat",  # input table for inversion below transition line
    "EoS_cross_input": "EoSCrossGeV.dat",  # input table for inversion crossing transition line
    "premerger_eos": True,  # whether to run premerger EoS inversion (it cleans unphysical regions first)
    "premerger_output": "output",  # output folder for premerger EoS inversion
    "EoS_above": "EoS_inv_above.dat",  # output table from inversion above transition line
    "EoS_below": "EoS_inv_below.dat",  # output table from inversion below transition line
    "EoS_cross": "EoS_inv_cross.dat",  # output table from inversion crossing transition line
    "TransitionLine": "TransitionLine.dat",  # path to transition line data
    "RegionS": "RegionS.dat",  # path to transition region S data
    "OutputMergedEoS": "invMergedEoS_new.dat",  # output merged EoS table
    "Ne": 400,  # number of points in energy density direction for merged EoS
    "Nn": 400,  # number of points in baryon density direction for merged EoS
}
