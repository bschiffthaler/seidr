.. _aggregate-label:

Aggregating networks into a crowd network
===========================================

Introduction
^^^^^^^^^^^^

Given a number of networks in ``SeidrFile`` format, seidr can aggregate those
into a crowd network. The basic syntax is::

  seidr aggregate <SeidrFile> <SeidrFile> ...

There are currently four methods of aggregation implemented:

* ``-m borda``: This will ouput a mean of ranks.
* ``-m top1``: This will ouput the edge with the highest score (==lowest rank) of all methods
* ``-m top2``: This will ouput the middle of the two highest scores (==lowest ranks) of all methods
* ``-m irp``: This will calculate the inverse rank product.

From a real example::

  seidr aggregate -m irp ../elnet/elnet_scores.sf ../narromi/narromi_scores.sf ../pearson/pearson_scores.sf ../spearman/spearman_scores.sf ../plsnet/plsnet_scores.sf ../aracne/aracne_scores.sf ../tigress/tigress_scores.sf ../clr/clr_scores.sf ../genenet/genenet_scores.sf ../svm/svm_scores.sf ../llr/llr_scores.sf ../genie3/genie3_scores.sf ../anova/anova_scores.sf

Without specifying an output file, this will create a file ``aggregated.sf`` in the
current working directory. Each column after the third (excluding the supplementary)
column stores the score and rank for each edge (if present) in all aggregated methods.
Converted to text (with ``seidr view``) the file looks like this::

  Source  Target  Type  ELNET_score;ELNET_rank  Narromi_score;Narromi_rank  Pearson_score;Pearson_rank  Spearman_score;Spearman_rank  PLSNET_score;PLSNET_rank  ARACNE_score;ARACNE_rank  TIGRESS_score;TIGRESS_rank  CLR_score;CLR_rank  PCor_score;PCor_rank  SVM_score;SVM_rank  LLR_score;LLR_rank  GENIE3_score;GENIE3_rank  ANOVA_score;ANOVA_rank  irp_score;irp_rank
  G2  G1  Undirected  0.004;334084  0.0128741;202752  -0.159435;202751  -0.00225177;1.32058e+06 1.07712e-05;360264nan;nan nan;nan 1.87357;106802  -0.018736;243746  0.152;26168 0.244;37455.5 0.0904447;42007 0.288087;1.30856e+06  0.176275;129253
  G3  G1  Undirected  0.334;22729.5 0.0381324;38394 -0.270978;44973 -0.214385;48864 3.2165e-05;61265  nan;nan 0.0028;78346.5  2.27349;70552.5 -0.021059;184389  0.077;91342.5 0.203;48670.5 0.215094;12249  0.388856;608154 0.299126;27713

We note that the final column stores the score of the aggregated network (IRP method).
For all future purposes, this is the representative score unless otherwise specified.

Directionality
^^^^^^^^^^^^^^

When there are directional algorithms in the aggregated set, seidr will attempt to create a directed aggregated edge as well. The exact procedure used depends on the aggregation algorithm:

* ``-m irp`` and ``-m borda``: When an edge is aggregated with the IRP or Borda schemes, the resulting edge takes the direction of the two-thirds majority of the input set *for methods that have called a direction*. Example: We aggregate 10 methods. 5 of these methods predict A->B, 2 methods predict A<-B, 3 methods predict A-B. 7 methods call a direction, therefore a 2/3 majority is reached with :math:`2 \cdot \lfloor\frac{7}{3}\rfloor = 4`. The 5 votes of A->B are therefore taken as the final direction.

* ``-m top1`` and ``-m top2``: For these methods the directionality of the selected edges is taken over into the final edge. For ``top2`` both methods have to agree, otherwise the edge is set to undirected.


Flags
^^^^^

When running seidr aggregate with the ``-k, --keep`` flag, supplementary data with flags for each algorithm will be kept. This suplementary data is a bitwise flag describing the edge. Currently 3 bits are used::

  000
    ^____ Edge exists
   ^_____ Edge is in direction A->B
  ^______ Edge is in direction A<-B

Resolving these bits to unsigned integers therefore results in:

* ``0``: Edge does not exist
* ``1``: Edge exists and is undirected
* ``3``: Edge exists and is A->B
* ``5``: Edge exists and is A<-B