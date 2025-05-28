.. _plsnet-label:

PLSNET
======

PLSNET uses a partial least squares feature selection algorithm to predict
interacting genes. It is published in [Guo2016]_ .

Running PLSNET
^^^^^^^^^^^^^^^

PLSNET needs a minimum of two input files:

* ``-i, --infile``: An expression matrix (genes are columns, samples are rows) without headers.
* ``-g, --genes``: A file containing gene names that correspond to columns in the expression matrix.

Here is an example matrix containing expression data for five genes in ten samples::

    0.4254475 0.0178292 0.9079888 0.4482474 0.1723238
    0.4424002 0.0505248 0.8693676 0.4458513 0.1733112
    1.0568470 0.2084539 0.4674478 0.5050774 0.2448833
    1.1172264 0.0030010 0.3176543 0.3872039 0.2537921
    0.9710677 0.0010565 0.3546514 0.4745322 0.2077183
    1.1393856 0.1220468 0.4024654 0.3484362 0.1686139
    1.0648694 0.1405077 0.4817628 0.4748571 0.1826433
    0.8761173 0.0738140 1.0582917 0.7303661 0.0536562
    1.2059661 0.1534070 0.7608608 0.6558457 0.1577311
    1.0006755 0.0789863 0.8036309 0.8389751 0.0883061

In the genes files, we provide the column headers for the expression matrix *in order*::

    G1
    G2
    G3
    G4
    G5

With that, we can run PCor::

    plsnet -i expr_mat.tsv -g genes.txt

The output is a square matrix of scores::

    0       10661.7 9103.01 3781.48 8553.33
    1672.03 0       3808.75 4130.81 24318.7
    2783.44 6850.63 0       2885.31 23882.1
    1683.27 11640.3 4560.34 0       63590
    1218.51 20635.5 9127.68 14218.4 0

Optional arguments for PLSNET
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``-S, --no-scale``: This turns off `feature scaling <https://en.wikipedia.org/wiki/Feature_scaling#Standardization>`_ of the expression matrix before the correlation calculation. By default scaling is enabled.
* ``-e, --ensemble``: Perform this many resampling iterations for each gene.
* ``-c, --components``: The number of PLS components to be considered.
* ``-p,  --predictor-size``: The number of predictors (genes) to be sampled at each iteration.

Running PLSNET for a subset of genes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Often we have only a small number of genes of interest. We can instruct 
PLSNET to only calculate interactions involving those genes by 
providing a ``-t, --targets`` file containing these gene names::

    G3
    G4

And running it with the ``-t, --targets`` options::

    plsnet -i expr_mat.tsv -g genes.txt -t targets.txt

In this case we will receive an edge list as output::

    G3  G1  1560.18
    G3  G2  892.019
    G3  G4  1471.69
    G3  G5  203666
    G4  G1  943.506
    G4  G2  1515.68
    G4  G3  5611.66
    G4  G5  542294

Running PLSNET in MPI mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PLSNET can use parallel processing. For general info
on how to run parallel algorithms in ``seidr``, please see :ref:`mpirun-label`
