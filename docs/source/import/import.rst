.. _import-label:

Importing text based formats into ``Seidr``
===========================================

Introduction
^^^^^^^^^^^^

Seidr works with its own binary file format, ``SeidrFile`` (see :ref:`seidrfile-label`).
In order to convert text based formats to a ``SeidrFile`` we use the ``seidr import``
command. Format conversion is not the only thing ``seidr import`` does, it also
ranks edge weights in the input files according to user parameters.

Import formats
^^^^^^^^^^^^^^

``seidr import`` currently supports three text based input formats.

**Lower triangular matrix** (``--format "lm"``)

A lower triangular matrix represents the lower half of a symmetric matrix. This
is particularly useful for non directional inference algorithms. Take Pearson
correlation as an example: If we correlate two vectors :math:`x, y`, it does
not matter if we check :math:`x \sim y` or :math:`y \sim x`. It would therefore
be a waste of space to store and a waste of computational resources to compute
the second comparison. A lower trinagular matrix will have the result of each
comparison exactly once::

      G1  G2  G3  G4  G5
  G1
  G2  0
  G3  1   2
  G5  3   4   5
  G5  6   7   8   9

All cells with an index in the above matrix exist in the lower triangular. Note that
the input is expected to be without headers, e.g.::

  0
  1   2
  3   4   5
  6   7   8   9

**Matrix** (``--format "m"``)

Opposed to the lower triangular, a matrix input is a square of all nodes vs
all nodes (including self-self, which is ignored). This is the output of
several machine learning algorithms which are non-symmetric (e.g. GENIE3, ELNET)::

      G1  G2  G3  G4  G5
  G1  0   1   2   3   4
  G2  5   6   7   8   9
  G3  10  11  12  13  14
  G4  15  16  17  18  19
  G5  20  21  22  23  24

Same as before, the input is expected to be without headers::

  0   1   2   3   4
  5   6   7   8   9
  10  11  12  13  14
  15  16  17  18  19
  20  21  22  23  24

**Edge lists** (``--format "el"``)

Edge lists are simple TAB separated files, which describe one edge one line. They
are relatively inefficient and store a lot of repetitive information, but they 
are very convenient for sparse networks and for humans to read::

  G1  G2  0
  G1  G3  1
  G1  G4  2
  G1  G5  3
  G2  G3  4
  G2  G4  5
  G2  G5  6
  G3  G4  7
  G3  G5  8
  G4  G5  9

**ARACNE2** (``--format "ara"``)

While ``seidr`` can read output in the format of the original ARACNE2 output,
the feature is not very well tested and should be considered experimental.

Importing your data
^^^^^^^^^^^^^^^^^^^

As a minimum, three arguments are required:

* ``-i, --infile``: The input text file
* ``-g, --genes`` : A file containing all genes (nodes) in the input file
* ``-F, --format``: The input format (``lm``, ``m``, ``el``)

Let's assume we have the lower triangular from before as output from our 
algorithm **lm.txt**::

  0
  1   2
  3   4   5
  6   7   8   9

We also have a file containg the names of nodes in the same order as the matrix.
Note that for the lower triangular this file is assumed to be sorted as if it
was column headers for the full (square) matrix. Generally, this is the same
as column headers for the input data matrix of the algorithms **nodes.txt**::

  G1
  G2
  G3
  G4
  G5

We can then run::

  seidr import -i lm.txt -g nodes.txt -F lm

Once it finishes, we can view the output with::

  $ seidr view elranks.sf
  G2  G1  Directed  0;1
  G3  G1  Directed  1;2
  G3  G2  Directed  2;3
  G4  G1  Directed  3;4
  G4  G2  Directed  4;5
  G4  G3  Directed  5;6
  G5  G1  Directed  6;7
  G5  G2  Directed  7;8
  G5  G3  Directed  8;9
  G5  G4  Directed  9;10

Adjusting import behaviour
^^^^^^^^^^^^^^^^^^^^^^^^^^

Depending on the algorithm, the default behaviour might need to be adjusted. In
the last example, we imported a lower triangular matrix, which by default creates
all directed edges. In many cases, this might not be true as the lower triangular
is likely to stem from a symmetric inference algorithm. The ``-u, --undirected``
option would do just that. Here are all modifiers:

* ``-u, --undirected``: Forces all edges to be interpreted as undirected. Use when source data is from a symmetric method
* ``-z, --drop-zero`` : Regards edges with a score of 0 as missing. Use for sparse methods.
* ``-r, --reverse``   : Considers higher edge weights better. Use when a higher score means a more confident prediction. Most methods implemented in ``seidr`` work that way, e.g. an edge weight of 0.6 is better than one of 0.2. If you import data from an algorithm that computes e.g. P-values, you need to omit this flag, as lower P-values are better.
* ``-A, --absolute``  : Computes the ranking using absolute values. A good example for this is Pearson correlation. Both 1 and -1 are perfect correlations, but they tell different stories. We want to keep the sign intact, but give both edges the highest rank for aggregation, therefore we use this flag.


Naming imports
^^^^^^^^^^^^^^

The last flag (``-n, --name``) lets you provide an internal name to the
``SeidrFile`` you are creating. Later, when you aggregate several ``SeidrFiles``
this will let you recognize the source of each score/rank column in the aggregated
network.

A note on parallelism
^^^^^^^^^^^^^^^^^^^^^

If ``seidr`` was compiled with a compiler that supports OpenMP, ``seidr import``
will carry out some steps in parallel. You can control how many CPUs it should 
use with the ``OMP_NUM_THREAD`` environment variable. If you would like to turn
multithreading for for example::

  OMP_NUM_THREAD=1 seidr import -f lm.txt -g nodes.txt -F lm ...