.. _anoverence-label:

ANOVERENCE
==========

**Please not that it is currently not recommended to run ANOVERENCE due to inconsitencies with the original implementation that we were not able to clarify with the original author**

ANOVERENCE ([Kueffner2012]_) employs the :math:`\eta^2` metric, a nonlinear correlation coefficient derived from an analysis of variance (ANOVA) ([Cohen1973]_). It is one
of the few methods that make direct use of experiment metadata, like perturbations,
knockouts and overexpressions.

Running ANOVERENCE
^^^^^^^^^^^^^^^^^^

ANOVERENCE needs a minimum of three input files:

* ``-i, --infile``: An expression matrix (genes are columns, samples are rows) without headers.
* ``-g, --genes``: A file containing gene names that correspond to columns in the expression matrix.
* ``-e, --features``: A file that contains experiment metadata.

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

The metadata file contains eight columns plus one row for each sample. If a 
column is not applicable, provide ``NA`` as input. Note that this file has
headers::

    Experiment Perturbations PerturbationLevels  Treatment DeletedGenes  OverexpressedGenes  Time  Repeat
    1 NA  NA  NA  NA  NA  NA  1
    1 NA  NA  NA  NA  NA  NA  2
    2 NA  NA  NA  NA  NA  NA  1
    3 NA  NA  NA  NA  NA  NA  1
    3 NA  NA  NA  NA  NA  NA  2
    4 NA  NA  NA  NA  NA  NA  1
    4 NA  NA  NA  NA  NA  NA  2
    5 NA  NA  NA  NA  NA  NA  1
    5 NA  NA  NA  NA  NA  NA  2
    5 NA  NA  NA  NA  NA  NA  3

Further we need to provide a ``-w, --weight``, typically an integral value between
10 and 1000 that controls how much more weight we give to perturbation experiments that involve the genes that are tested. Once we have those four parameters, we are
ready to run ``anoverence``::

    anoverence -i exr_mat.tsv -g genes.txt -e meta.tsv -w 500

As output we receive a lower triangular matrix of interaction scores::

    0.288087
    0.388856        0.405731
    0.459865        0.276648        0.336653
    0.432748        0.374432        0.397973        0.403535

Running ANOVERENCE for a subset of genes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Often we have only a small number of genes of interest. We can instruct 
``ANOVERENCE`` to only calculate interactions involving those genes by 
providing a ``-t, --targets`` file containing these gene names::

    G3
    G4

And running it with the ``-t, --targets`` options::

    anoverence -i expr_mat.tsv -g genes.txt -e meta.tsv -w 500 -t targets.txt

In this case we will receive an edge list as output::

    G3  G1  0.388856
    G4  G1  0.459865
    G3  G2  0.405731
    G4  G2  0.276648
    G4  G3  0.336653
    G3  G5  0.397973 
    G4  G5  0.403535