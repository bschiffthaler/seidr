.. _getting-started-label:

Getting Started
===============

Tutorial Data
^^^^^^^^^^^^^

The tutorial data is a set of 500 `salmon <https://combine-lab.github.io/salmon/>`_ pseudo-alignments targeting *Saccharomyces cerevisiae*. They originate from public `SRA <https://www.ncbi.nlm.nih.gov/sra/>`_ data. In fact the files names are simply their SRA accessions. In this tutorial, we will go through a simple pre-processing step and then create individual and aggregate networks. First, let's get the data:

.. code-block:: Bash
  
  mkdir seidr_tutorial
  cd seidr_tutorial
  wget https://bschiffthaler.s3-eu-west-1.amazonaws.com/SeidrPublic/seidr_tutorial.tar.gz
  tar -xvf seidr_tutorial.tar.gz

Pre-processing
^^^^^^^^^^^^^^

Data pre-processing is done in R. We'll use `tximport <https://bioconductor.org/packages/release/bioc/html/tximport.html>`_ to load the data into R, and then use `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ to variance-stabilize.

.. code-block:: R

  library(tximport)
  library(DESeq2)
  library(readr)

  # This load the mapping of genes to transcripts into R so that tximport can
  # summarize counts
  tx2g <- read_tsv('seidr_tutorial/tx2gene.tsv',
                   col_names=c('Transcript', 'Gene'))

  # Now let's find all our count files and load them into R using tximport
  input_files <- dir('seidr_tutorial', pattern='*_quant.sf$', full.names=TRUE)
  txi <- tximport(input_files, 'salmon', tx2gene=tx2g)

  # In order to make a DESeq2 data set, we need some metadata. For now we'll
  # just use some dummy data
  dummy_meta <- data.frame(N = seq_along(input_files))
  dds <- DESeqDataSetFromTximport(txi, dummy_meta, ~1)

  # Now we run variance stabilization and get the stabilized data as a matrix. If
  # you have good metadata, you can use the experimental design in the DESeqDataSet
  # and set blind=FALSE here
  vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
  vst <- t(assay(vsd))

  # Genes that do not vary at all create problems down the line, so it's better
  # to drop them
  vars <- apply(vst, 2, var)
  filt_id <- which(is.finite(vars))
  vst <- vst[, filt_id]

  # Let's also center samples around their median, which has been shown to
  # improve reconstruction accuracy
  medians <- apply(vst, 1, median)
  vst <- sweep(vst, MARGIN=1, FUN='-', STATS=medians)

  # MASS's write.matrix function is a bit faster and better suited for our task
  # compared to write.table. Don't forgot to unname(), otherwise you will have
  # column headers in the output
  MASS::write.matrix(x=unname(vst), sep='\t', file='expression.tsv')

  # Finally, let's write the column headers (== gene names) as a text file
  write(colnames(vst), file="genes.txt")

All versus all networks
^^^^^^^^^^^^^^^^^^^^^^^

Sub-setting the data
""""""""""""""""""""

In this section, we will create an all vs. all comparison, meaning we will estimate connectivity of all genes to all other genes. This approach is the most resource demanding, so we'll create a smaller subset of the tutorial data first. We'll take 100 samples and the first 1000 genes. Our output network will therefore be :math:`\frac{1000 \cdot 999}{2}`.

.. code-block:: Bash

  tail -n 100 expression.tsv | cut -f 1-1000 > expression_sub.tsv
  head -n 1000 genes.txt > genes_sub.txt

Network inference
"""""""""""""""""

Now that we have a data subset, we can get started with the inference. In this step, we'll create 13 different all vs. all networks using algorithms that seidr ships with. If you have any inference algorithm you would like to include that is not yet implemented in seidr, you can run that as well, but make sure its output is in a format seidr can import (:ref:`import-label`). Even though this is a sub-set, you'll probably need to set aside an hour for the inference. If you want a quicker run, you can leave out ``el-ensemble``.

.. code-block:: Bash

  # fast
  correlation -m pearson -i expression_sub.tsv -g genes_sub.txt
  correlation -m spearman -i expression_sub.tsv -g genes_sub.txt
  pcor -i expression_sub.tsv -g genes_sub.txt
  
  # medium
  mi -m RAW -i expression_sub.tsv -g genes_sub.txt -o mi_scores.tsv
  mi -m CLR -i expression_sub.tsv -g genes_sub.txt -M mi_scores.tsv -o clr_scores.tsv
  mi -m ARACNE -i expression_sub.tsv -g genes_sub.txt -M mi_scores.tsv -o aracne_scores.tsv

  # slow
  narromi -m interior-point -i expression_sub.tsv -g genes_sub.txt -o narromi_scores.tsv
  plsnet -i expression_sub.tsv -g genes_sub.txt -o plsnet_scores.tsv
  llr-ensemble -i expression_sub.tsv -g genes_sub.txt -o llr_scores.tsv
  svm-ensemble -k POLY -i expression_sub.tsv -g genes_sub.txt -o svm_scores.tsv
  genie3 -i expression_sub.tsv -g genes_sub.txt -o genie3_scores.tsv
  tigress -i expression_sub.tsv -g genes_sub.txt -o tigress_scores.tsv

  # very slow
  el-ensemble -i expression_sub.tsv -g genes_sub.txt -o elnet_scores.tsv

Network ranking
"""""""""""""""

Different inference algorithms output networks with different metrics for edge weights. A correlation network, will assign scores anywhere in :math:`[-1, ..., 1]`, whereas mututal information is in :math:`[0, ..., N]`, and many of the regression algorithms in :math:`[0, ..., 1]`. We can therefore not just sum the weights to get a final, community network. In order to do that, we want to convert the scores to `ranks <https://en.wikipedia.org/wiki/Ranking#Ranking_in_statistics>`_. The command ``seidr import`` takes care of that. Some useful options:

* ``-A``: This option computes the rank on the **absolute** value of the score, so -1 and +1 would get the same rank.
* ``-r``: This option indicates that **higher scores are better**. A score of 1 would get a lower (== better) rank than a score of 0.5.
* ``-u``: This option creates an undirected network. We use this in algorithms where we know the output is symmetric (A->B and B->A are the same), but only have the lower triangular matrix. Examples are all correlation and mutual information based methods.
* ``-z``: This option drops edges with a score of 0. By default we keep all edges, but this will create sparser networks for methods that output 0-valued edges. 

.. code-block:: Bash

  seidr import -A -r -u -n PEARSON -o person_scores.sf -F lm -i pearson_scores.tsv -g genes_sub.txt
  seidr import -A -r -u -n SPEARMAN -o spearman_scores.sf -F lm -i spearman_scores.tsv -g genes_sub.txt
  seidr import -A -r -u -n PCOR -o pcor_scores.sf -F lm -i pcor_scores.tsv -g genes_sub.txt

  seidr import -r -u -n MI -o mi_scores.sf -F lm -i mi_scores.tsv -g genes_sub.txt
  seidr import -r -u -z -n CLR -o clr_scores.sf -F lm -i clr_scores.tsv -g genes_sub.txt
  seidr import -r -u -z -n ARACNE -o aracne_scores.sf -F lm -i aracne_scores.tsv -g genes_sub.txt

  seidr import -r -z -n NARROMI -o narromi_scores.sf -F m -i narromi_scores.tsv -g genes_sub.txt
  seidr import -r -z -n PLSNET -o plsnet_scores.sf -F m -i plsnet_scores.tsv -g genes_sub.txt
  seidr import -r -z -n LLR -o llr_scores.sf -F m -i llr_scores.tsv -g genes_sub.txt
  seidr import -r -z -n SVM -o svm_scores.sf -F m -i svm_scores.tsv -g genes_sub.txt
  seidr import -r -z -n GENIE3 -o genie3_scores.sf -F m -i genie3_scores.tsv -g genes_sub.txt
  seidr import -r -z -n TIGRESS -o tigress_scores.sf -F m -i tigress_scores.tsv -g genes_sub.txt
  seidr import -r -z -n ELNET -o elnet_scores.sf -F m -i elnet_scores.tsv -g genes_sub.txt

Aggregating
"""""""""""

Aggregating refers to the construction of a community network from the individual networks created before. Note that there are several aggregation methods available. We will use the "Inverse Rank Product" method described in [Zhong2014]_.

.. code-block:: Bash
  
  seidr aggregate -m irp aracne_scores.sf clr_scores.sf elnet_scores.sf genie3_scores.sf llr_scores.sf mi_scores.sf narromi_scores.sf pcor_scores.sf person_scores.sf plsnet_scores.sf spearman_scores.sf svm_scores.sf tigress_scores.sf

This creates a community network of all the 1000 genes in our sample data. If you don't want to learn how you can create a network for a group of genes (e.g. only transcription factors), jump right to :ref:`post-processing-label`.

Taking a look at the final network
""""""""""""""""""""""""""""""""""

We can have a look at the top three edges in the network:

.. code-block:: Bash

  seidr top -n 3 aggregated.sf | column -t

.. code-block:: none

  YDL039C    YDL037C    Undirected  0.777779;83  14.1252;1   0.511;433  3.61084;7   0.138;61223.5  0.777779;204  0.709272;39   0.0787986;93     0.837201;940  1.5595;25      0.780375;1512  nan;nan        1;1.5       0.949186;3
  YDL025W-A  YBL006W-A  Undirected  1.05879;2    8.3004;29   0.543;4    3.30464;21  0.496;2775     1.05879;2     0.372718;529  0.0172536;31848  0.909562;219  0.654752;1758  0.847192;364   0.463;4415     0.8677;144  0.967855;2
  YAL037C-B  YCR013C    Directed    1.00539;7    9.08174;14  0.519;168  2.82549;92  0.517;407      1.00539;7     0.777048;28   0.032911;3601    0.928263;137  0.991901;297   0.846061;373   0.204;12028.5  0.99325;16  1;1

Most of these are from dubious ORFs (which should have maybe been filtered beforehand). The one that is not, is definitely a good result, YDL039C and YDL037C as both these genes form the `IMI1 protein <https://www.yeastgenome.org/locus/S000149345>`_. 

Creating targeted networks
^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes, we are not interested in the interactions of all genes, we just want to know what our genes of interest look like in the network. We can then run seidr in targeted mode, which will compute only what's necessary to understand that particular group of genes. The slowest of the bunch will probably be the mutual information based algorithms CLR and ARACNe, since they are context dependent and the full mutual information matrix needs to be computed first.

Making targets
""""""""""""""

Since we need some targets to look at, I select a single transcription factor FZF1, and store the gene identifier in a file.

.. code-block:: Bash

  echo "YGL254W" > FZF1.txt

Inferring sub-networks
""""""""""""""""""""""

The network inference step is nearly the same, but now we use the full expression set (all ~6500 genes and 500 samples) as well as the ``FZF1.txt`` targets file.

.. code-block:: Bash

  # fast
  correlation -t FZF1.txt -m pearson -i expression.tsv -g genes.txt --scale -o pearson_fzf1_scores.tsv
  correlation -t FZF1.txt -m spearman -i expression.tsv -g genes.txt -o spearman_fzf1_scores.tsv
  pcor -t FZF1.txt -i expression.tsv -g genes.txt --scale -o pcor_fzf1_scores.tsv
  
  # medium
  mi -t FZF1.txt -m RAW -i expression.tsv -g genes.txt -M mi_full_scores.tsv -o mi_fzf1_scores.tsv
  mi -t FZF1.txt -m CLR -i expression.tsv -g genes.txt -M mi_full_scores.tsv -o clr_fzf1_scores.tsv
  mi -t FZF1.txt -m ARACNE -i expression.tsv -g genes.txt -M mi_full_scores.tsv -o aracne_fzf1_scores.tsv

  # slow
  narromi -t FZF1.txt -m interior-point -i expression.tsv -g genes.txt -o narromi_fzf1_scores.tsv
  plsnet -t FZF1.txt -i expression.tsv -g genes.txt -o plsnet_fzf1_scores.tsv --scale
  llr-ensemble -t FZF1.txt -i expression.tsv -g genes.txt -o llr_fzf1_scores.tsv --scale
  svm-ensemble -t FZF1.txt -k POLY -i expression.tsv -g genes.txt -o svm_fzf1_scores.tsv --scale
  genie3 -t FZF1.txt -i expression.tsv -g genes.txt -o genie3_fzf1_scores.tsv --scale
  tigress -t FZF1.txt -i expression.tsv -g genes.txt -o tigress_fzf1_scores.tsv --scale

  el-ensemble -t FZF1.txt -i expression.tsv -g genes.txt -o elnet_fzf1_scores.tsv --scale


Importing
"""""""""

Targeted mode outputs results in edge list format, so all out imports now contain ``-F el`` instead of ``-F lm`` or ``-F m``.

.. code-block:: Bash

  seidr import -A -r -u -n PEARSON -o person_fzf1_scores.sf -F el -i pearson_fzf1_scores.tsv -g genes.txt
  seidr import -A -r -u -n SPEARMAN -o spearman_fzf1_scores.sf -F el -i spearman_fzf1_scores.tsv -g genes.txt
  seidr import -A -r -u -n PCOR -o pcor_fzf1_scores.sf -F el -i pcor_fzf1_scores.tsv -g genes.txt

  seidr import -r -u -n MI -o mi_fzf1_scores.sf -F el -i mi_fzf1_scores.tsv -g genes.txt
  seidr import -r -u -z -n CLR -o clr_fzf1_scores.sf -F el -i clr_fzf1_scores.tsv -g genes.txt
  seidr import -r -u -z -n ARACNE -o aracne_fzf1_scores.sf -F el -i aracne_fzf1_scores.tsv -g genes.txt

  seidr import -r -z -n NARROMI -o narromi_fzf1_scores.sf -F el -i narromi_fzf1_scores.tsv -g genes.txt
  seidr import -r -z -n PLSNET -o plsnet_fzf1_scores.sf -F el -i plsnet_fzf1_scores.tsv -g genes.txt
  seidr import -r -z -n LLR -o llr_fzf1_scores.sf -F el -i llr_fzf1_scores.tsv -g genes.txt
  seidr import -r -z -n SVM -o svm_fzf1_scores.sf -F el -i svm_fzf1_scores.tsv -g genes.txt
  seidr import -r -z -n GENIE3 -o genie3_fzf1_scores.sf -F el -i genie3_fzf1_scores.tsv -g genes.txt
  seidr import -r -z -n TIGRESS -o tigress_fzf1_scores.sf -F el -i tigress_fzf1_scores.tsv -g genes.txt
  seidr import -r -z -n ELNET -o elnet_fzf1_scores.sf -F el -i elnet_fzf1_scores.tsv -g genes.txt

Aggregating
"""""""""""

This is exactly the same as for the full network.

.. code-block:: Bash

  seidr aggregate -m irp -o aggregated_fzf1.sf aracne_fzf1_scores.sf clr_fzf1_scores.sf elnet_fzf1_scores.sf genie3_fzf1_scores.sf llr_fzf1_scores.sf mi_fzf1_scores.sf narromi_fzf1_scores.sf pcor_fzf1_scores.sf person_fzf1_scores.sf plsnet_fzf1_scores.sf spearman_fzf1_scores.sf svm_fzf1_scores.sf tigress_fzf1_scores.sf

Taking a look at the final network
""""""""""""""""""""""""""""""""""

Just as before, let's look at the top three connections of our TF.

.. code-block:: Bash

  seidr top -n 3 aggregated_fzf1.sf

.. code-block:: none

  YGL254W  YEL051W  Directed  nan;nan  3.11877;43     0.458;2    1.81032;4     0.004;3226  0.387233;48    nan;nan       -0.00566415;1138  -0.57406;1    0.035306;28   -0.462468;6  0.001;3050.5  0.74525;2  0.874875;3
  YGL254W  YGL128C  Directed  nan;nan  0.886213;1151  0.443;4.5  0.266012;279  0.248;86    0.265754;1222  0.0648351;43  0.010073;254      0.453637;131  0.0382971;21  0.404141;49  0.207;210.5   0.19945;7  0.888152;2
  YGL254W  YFL044C  Directed  nan;nan  3.24149;31     0.451;3    0.581017;113  0.352;17    0.383645;52    0.1711;4      0.00935282;329    0.477231;78   0.0320321;45  0.444953;11  0.041;1019    0.25265;4  1;1

The top connections are OTU1 (YFL044C), CWC23 (YGL128C), and VMA8 (YEL051W).

.. _post-processing-label:

Post processing
^^^^^^^^^^^^^^^

Pruning noisy edges
"""""""""""""""""""

In most cases, the community network will be fully dense, meaning every gene is connected to every other gene with a certain score. Many of these edges are just noise and we would like to prune them. [Coscia2017]_ have developed a smart approach to pruning noisy edges called "Network backboning". We can apply this to our community network as:

.. code-block:: Bash

  seidr backbone -F 1.28 aggregated.sf

Viewing edges in the network
""""""""""""""""""""""""""""

The ``seidr view`` command offers an interface to query the seidr output. Let's look at a few edges.

.. code-block:: Bash

  seidr view --column-headers aggregated.bb.sf | head -n 3 | column -t

.. code-block:: none
  
  Source  Target  Type        ARACNE_score;ARACNE_rank  CLR_score;CLR_rank  ELNET_score;ELNET_rank  GENIE3_score;GENIE3_rank  LLR_score;LLR_rank  MI_score;MI_rank  NARROMI_score;NARROMI_rank  PCOR_score;PCOR_rank  PEARSON_score;PEARSON_rank  PLSNET_score;PLSNET_rank  SPEARMAN_score;SPEARMAN_rank  SVM_score;SVM_rank  TIGRESS_score;TIGRESS_rank  irp_score;irp_rank  NC_Score;NC_SDev;SEC;EBC
  Q0017   Q0010   Undirected  nan;nan                   3.58054;6569        0.317;14818             0.587286;20930            0.221;43021         0.255911;100939   0.364063;574                0.0246169;10270       0.644044;14824              1.01683;268               0.509318;33915                0.154;13896         0.17055;3855.5              0.440763;1637       0.602226;0.375774;0.459647;130
  Q0032   Q0010   Undirected  nan;nan                   2.38815;27858       0.269;18028             0.905035;9646             0.138;61223.5       0.138116;379828   nan;nan                     0.0481595;843         0.679379;10339              1.19476;122               0.386693;81327                0.287;9567          0.0787;8143.5               0.37659;3358        0.621852;0.388965;0.364875;218

Querying specific nodes or edges
""""""""""""""""""""""""""""""""

If the seidr output is indexed with the ``seidr index`` command, we can query specific nodes and edges.

.. code-block:: Bash

  seidr index aggregated.bb.sf
  # Node
  seidr view -n YBR142W aggregated.bb.sf
  # Edge
  seidr view -n YBR142W:YDL063C aggregated.bb.sf

Graph and centrality statistics
"""""""""""""""""""""""""""""""

Seidr can compute statistics on the entire graph and some node centrality measures. Before we do that, it's best to make sure we have no disconnected nodes in the graph, which we drop with::

  seidr reheader aggregated.bb.sf

Then, we can use ``seidr graphstats`` to compute graph summary stats.

.. code-block:: Bash
  
  seidr graphstats aggregated.bb.sf

.. code-block:: none

  Number of Nodes:        974
  Number of Edges:        4150
  Number of Connected Components: 2
  Global clustering coefficient:  0.338051
  Scale free fit: 0.0555276
  Average degree: 8.52156
  Average weighted degree:        3.71159
  Network diameter:       4.00379
  Average path length:    1.66305

Finally, we can compute node centrality statistics with ``seidr stats``

.. code-block:: Bash

  seidr stats --exact aggregated.bb.sf
  seidr view --centrality aggregated.bb.sf | sort -k2g | tail -n 5 | column -t

.. code-block:: none

  YBL006W-A  0.002748    1344   14.8438  0.140839     0.0339691
  YBL039C    0.00279961  1322   14.4697  0.000898011  0.0339036
  YBR142W    0.00282644  6460   14.3319  0.00113027   0.03388
  YBL012C    0.00296141  10198  16.4288  0.156134     0.0342437
  YAL045C    0.00306967  4364   17.8655  0.234332     0.0344952
