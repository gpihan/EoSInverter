general_parameters = {
        "RunMode":1,   # Run Mode: 0 -> Run locally
                       # Run Mode: 1 -> Run on slurm cluster

        "Dimension":4, # Dimensions = 1 gives T(e)
                       # Dimensions = 2 gives T(e, nb), mub(T, nb)
                       # Dimensions = 3 gives T(e, nb, nq), mub(e, nb, nq), muq(e, nb, nq)
                       # Dimensions = 4 gives T(e, nb, nq, ns), mub(e, nb, nq, ns), muq(e, nb, nq, ns), mus(e, nb, nq, ns)

        "AutoSetBoundaries":True,
        "NT":2,
        "NB":2,
        "NQ":2,
        "NS":2,

        "Ttilde": [0.009777504201390019, 0.31403115649771435, 10], 
        "muBtilde":[-2.2732322158914076, 2.2732322158914076, 10],
        "muQtilde":[-0.6292193138512704, 0.6292193138512704, 10],
        "muStilde":[-0.6880440931493973, 0.6880440931493973, 10],
        "Accuracy": 1e-6,
        "MAXITER": 100,
        "Number_of_cores": 12,
        "EoS_table": "HRG4D.dat",
        "OutputFolder":"OutTest"
}
