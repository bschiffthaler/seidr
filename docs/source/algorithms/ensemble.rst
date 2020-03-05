.. _ensemble-label:

Ensemble methods
================

The ensemble methods are based on [Ruyssinck2014]_ . The three main executables
work the same way and have the same options. They all work by resampling the expression data along samples and genes, which often reduces variance in their predictions:

* ``el-ensemble``: Uses an ensemble of Elastic Net regression predictors.
* ``svm-ensemble``: Uses an ensemble of Support Vector Machine predictors.
* ``llr-ensemble``: Uses an ensemble of Support Vector Machine predictors.

The Elastic Net code uses the GLMNET Fortran backend from [Friedman2010]_ .

Running Ensembles
^^^^^^^^^^^^^^^^^^

Each ensemble needs a minimum of two input files:

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

With that, we can run the ensembles::

    el-ensemble -i expr_mat.tsv -g genes.txt
    svm-ensemble -i expr_mat.tsv -g genes.txt
    llr-ensemble -i expr_mat.tsv -g genes.txt

The output is a square matrix of scores::

    0       0   0.876   0.124   0
    0.894   0   0.106   0       0
    0.894   0   0       0.106   0
    0.894   0   0.106   0       0
    0.894   0   0.106   0       0


Optional arguments for the Ensemble methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``-s, --scale``: This triggers `feature scaling <https://en.wikipedia.org/wiki/Feature_scaling#Standardization>`_ of the expression matrix before the regression calculation. Generally this should be *on*.
* ``-X, --max-experiment-size``: In each resampling iteration, choose maximally this many samples along rows (experiments) of the dataset.
* ``-x, --min-experiment-size``: In each resampling iteration, choose minimally this many samples along rows (experiments) of the dataset.
* ``-P, --max-predictor-size``: In each resampling iteration, choose maximally this many genes along columns (predictors) of the dataset.
* ``-p, --min-predictor-size``: In each resampling iteration, choose minimally this many genes along columns (predictors) of the dataset.
* ``-e, --ensemble``: Perform this many resampling iterations for each gene.

The sampling boundaries ``-X``, ``-x``, ``-P`` and ``-p`` default to 4/5th of 
samples/predictors for the upper bound and 1/5th for the lower. In runs with small
experiment sizes (<50) one should set this manually higher to avoid undersampling.
In these cases, I suggest 90% for the upper boundary and 50% for the lower (in experiments).

Running ensembles for a subset of genes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Often we have only a small number of genes of interest. We can instruct 
the ensembles to only calculate interactions involving those genes by 
providing a ``-t, --targets`` file containing these gene names::

    G3
    G4

And running it with the ``-t, --targets`` options::

    llr-ensemble -i expr_mat.tsv -g genes.txt -t targets.txt
    svm-ensemble -i expr_mat.tsv -g genes.txt -t targets.txt
    el-ensemble -i expr_mat.tsv -g genes.txt -t targets.txt

In this case we will receive an edge list as output::

    G3  G1  0.894
    G3  G2  0
    G3  G4  0.106
    G3  G5  0
    G4  G1  0.894
    G4  G2  0
    G4  G3  0.106
    G4  G5  0

Running Ensembles in MPI mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each ensemble can use parallel processing. For general info
on how to run parallel algorithms in ``seidr``, please see :ref:`mpirun-label`

The difference between SVM and LLR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

LLR and SVM are based on different implementations of SVMs in C. One is based on
`LibLinear <https://www.csie.ntu.edu.tw/~cjlin/liblinear/>`_ , the other on 
`LibSVM <https://www.csie.ntu.edu.tw/~cjlin/libsvm>`_ using a linear kernel. While 
they should in general agree most of the time, coefficients are handled differently.
SVM is closer to the reference implementation by [Ruyssinck2014]_ , but LLR is 
much faster.

