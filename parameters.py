general_parameters = {
    "RunMode": 2,  # Run Mode: 0 -> Run locally
    # Run Mode: 1 -> Run on slurm cluster
    "Dimension": 2,  # Dimensions = 1 gives T(e)
    # Dimensions = 2 gives T(e, nb), mub(T, nb)
    # Dimensions = 3 gives T(e, nb, nq), mub(e, nb, nq), muq(e, nb, nq)
    # Dimensions = 4 gives T(e, nb, nq, ns), mub(e, nb, nq, ns), muq(e, nb, nq, ns), mus(e, nb, nq, ns)
    "AutoSetBoundaries": True,
    "NT": 500,
    "NB": 500,
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
    "EoS_table": "EoS2DISINGnew.dat",
    "OutputFolder": "EoS2DISINGnew",
    "TransitionLine": "TransitionlineFilePar.dat",  # Path to transition line data (muB, T format)
    "RegionS": "RegionS.dat",  # path to transition region S data
    "OutputMergedEoS": "invIsing-2DTExS.dat",  # output merged EoS table
    "NTildeT": 500,  # number of points in energy density direction for merged EoS
    "NTildemuB": 500,  # number of points in baryon density direction for merged ES
}
