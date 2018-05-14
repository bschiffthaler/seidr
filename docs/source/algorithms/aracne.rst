.. _aracne-label:

ARACNE
==========

ARACNE ([Margolin2006]_) is an inference algorithm based on mutual information
and applies data processing inequality to delete most indirect edges.

Our implementation differs to the original in that it estimates the initial
mutual information using a B-spline approach as described in [Daub2004]_ .

Running ARACNE
^^^^^^^^^^^^^^^^^^

ARACNE needs a minimum of three input files:

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

With that, we can run ARACNE::

    mi -m ARACNE -i expr_mat.tsv -g genes.txt

The output is a lower triangular matrix of scores::

    0
    0.798215    0.874873
    0   0.889133    0
    0   0   0.860645    0.95965


Tuning the number of bins and spline degree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Estimating mutual infofrmation from discrete data is well defined, but normalized
expression data is usually continuous. To estimate the MI from continuous data, each
data point is usually assigned to one bin. This can lead to a loss of information.

The B-Spline estimator for MI therefore performs fuzzy assignment of the data to 
bins. The ``-s, --spline`` parameter controls the spline degree (therefore 
the shape) of the indicator function. For ``s=1`` the indicator function is the
same as for simple binning. Improvements in the MI beyond a degree of ``s=3``
are rarely seen, therefore it is a good choice as a default.

The number of bins used in the assignment can be controlled with the ``-b, --bins``
option. By default it is automatically inferred from the data, but this can lead
to high memory requirements on large datasets. Generally, the number of bins is
assumed not to influence the MI much as long as it's within a reasonable range. A
value between 5 and 10 is a good starting point for typically sized datasets from RNA-Seq.

Running ARACNE for a subset of genes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Often we have only a small number of genes of interest. We can instruct 
ARACNE to only calculate interactions involving those genes by 
providing a ``-t, --targets`` file containing these gene names::

    G3
    G4

And running it with the ``-t, --targets`` options::

    mi -m ARACNE -i expr_mat.tsv -g genes.txt -t targets.txt

In this case we will receive an edge list as output::

    G3  G1  0.798215
    G3  G2  0.874873
    G3  G4  0
    G3  G5  0.860645
    G4  G1  0
    G4  G2  0.889133
    G4  G3  0
    G4  G5  0.95965

Running ARACNE in MPI mode
^^^^^^^^^^^^^^^^^^^^^^^^^^

ARACNE can use parallel processing in the MI estimation step. For general info
on how to run parallel algorithms in ``seidr``, please see this page: 

References
^^^^^^^^^^

.. [Margolin2006] Margolin, A. A., Nemenman, I., Basso, K., Wiggins, C., Stolovitzky, G., Dalla Favera, R., & Califano, A. (2006). ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. BMC Bioinformatics, 7 Suppl 1, S7. https://doi.org/10.1186/1471-2105-7-S1-S7
.. [Daub2004] Daub, C. O., Steuer, R., Selbig, J., & Kloska, S. (2004). Estimating mutual information using B-spline functions--an improved similarity measure for analysing gene expression data. BMC Bioinformatics, 5, 118. https://doi.org/10.1186/1471-2105-5-118
