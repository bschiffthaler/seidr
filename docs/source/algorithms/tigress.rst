.. _tigress-label:

TIGRESS
=======

TIGRESS uses an ensemble approach (here called stability selection) to reduce prediction
variance in a LASSO model. It works somewhat similar to the other :ref:`ensemble-label` .
TIGRESS is published in [Haury2012]_ . The LASSO uses the GLMNET Fortran backend in 
[Friedman2010]_ .

Running TIGRESS
^^^^^^^^^^^^^^^

TIGRESS needs a minimum of two input files:

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

    tigress -i expr_mat.tsv -g genes.txt

The output is a square matrix of scores::

    0       0.5197  0.6558  0.2139  0.0838
    0.6909  0       0.169   0.2216  0.3628
    0.4819  0.1774  0       0.3084  0.7075
    0.137   0.1418  0.3349  0       0.675
    0.0182  0.1736  0.7209  0.6138  0

Optional arguments for TIGRESS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``-S, --no-scale``: This turns off `feature scaling <https://en.wikipedia.org/wiki/Feature_scaling#Standardization>`_ of the expression matrix before the correlation calculation. By default scaling is enabled.
* ``-B, --nbootstrap``: Perform this many resampling iterations for each gene.
* ``-n, --nlambda``: Consider this many shrinkage lambdas.
* ``-l, --min-lambda``: The minimum lambda value considered is this fraction of the maximum.

Running TIGRESS for a subset of genes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Often we have only a small number of genes of interest. We can instruct 
TIGRESS to only calculate interactions involving those genes by 
providing a ``-t, --targets`` file containing these gene names::

    G3
    G4

And running it with the ``-t, --targets`` options::

    tigress -i expr_mat.tsv -g genes.txt -t targets.txt

In this case we will receive an edge list as output::

    G3  G1  0.4819
    G3  G2  0.1774
    G3  G4  0.3084
    G3  G5  0.7075
    G4  G1  0.137
    G4  G2  0.1418
    G4  G3  0.3349
    G4  G5  0.675

Running TIGRESS in MPI mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TIGRESS can use parallel processing. For general info
on how to run parallel algorithms in ``seidr``, please see :ref:`mpirun-label`