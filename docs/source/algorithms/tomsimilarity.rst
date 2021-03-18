.. _tom-label:

TOM Similarity
==============

The toplogical overlap matrix (TOM) is the similarity measure implemented by
WGCNA [Langfelder2008]_. It calculates a correlation matrix from the expression 
data, calculates a soft threshold and assigns two genes a high topological
overlap if they share common neighbourhoods.

Running TOM similarity
^^^^^^^^^^^^^^^^^^^^^^

``tomsimilarity`` needs a minimum of two input files:

* ``-i, --infile``: An expression matrix (genes are columns, samples are rows) without headers.
* ``-g, --genes``: A file containing gene names that correspond to columns in the expression matrix.

Here is an example matrix containing expression data for five genes in ten samples::

    6.107967  7.188796  7.139945  9.417835   6.195927
    8.602925  9.134458  8.630118  10.695973  6.930023
    6.699199  8.307864  8.174942  10.874148  7.143233
    7.661777  8.891523  8.348661  10.439793  6.868748
    7.031853  9.019152  8.539557  10.726523  7.461354
    8.931517  9.246769  8.944240  10.774747  6.729316
    6.815357  9.209684  8.607074  9.574451   7.400409
    7.424712  9.603071  8.347164  10.609222  7.168921
    8.465108  8.788967  8.875855  10.537852  6.628380
    8.559188  8.992996  8.279209  10.640245  6.744078

In the genes files, we provide the column headers for the expression matrix *in order*::

    G1
    G2
    G3
    G4
    G5

With that, we can run PCor::

    tomsimilarity -i expr_mat.tsv -g genes.txt -b 4

The output is a lower triangular matrix of scores::

    0.44357
    0.486974  0.504881
    0.370446  0.408224  0.42039
    0.225011  0.465292  0.396999  0.252425

Optional arguments for ``tomsimilarity``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``-S, --no-scale``: This turns off `feature scaling <https://en.wikipedia.org/wiki/Feature_scaling#Standardization>`_ of the expression matrix before the correlation calculation. By default scaling is enabled.
* ``-m, --method``: Choose between "pearson" or "bicor" (`biweight midcorrelation <https://en.wikipedia.org/wiki/Biweight_midcorrelation>`_. The latter is typically a good choice unless you have a lot of outliers.)
* ``-b, --sft``: The soft threshold power. This is the exponent for soft thresholding the correlation matrix. Unless you know why, leave it default.
* ``-M, --max-power``: When auto-detecting the soft threshold power, this is the maximum value that will be tested. It's usually not a good idea to go above 30. If you cannot get a good fit, decrease the cutoff instead.
* ``-S, --sft-cutoff``: When the network reaches this scale free fit R^2 value, stop testing powers. Sometimes, you cannot get a good fit (>0.8) on larger datasets. In this case, decrease this value.
* ``-T, --tom-type``: "unsigned", "signed", or "signed-hybrid". This defines how to score the TOM. "unsigned" is :math:`\vert a_{ij} \vert`, "signed" is :math:`\frac{a_{ij} + 1}{2}` and "signed-hybrid" is :math:`\vert a_{ij} \vert` for positive correlation, 0 otherwise.

Running ``tomsimilarity`` for a subset of genes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Often we have only a small number of genes of interest. We can instruct 
``tomsimilarity`` to only calculate interactions involving those genes by 
providing a ``-t, --targets`` file containing these gene names::

    G3
    G4

And running it with the ``-t, --targets`` options::

    tomsimilarity -i expr_mat.tsv -g genes.txt -t targets.txt -b 4

In this case we will receive an edge list as output::

    G3  G1  0.486974
    G3  G2  0.504881
    G3  G4  0.42039
    G3  G5  0.396999
    G4  G1  0.370446
    G4  G2  0.408224
    G4  G3  0.42039
    G4  G5  0.252425