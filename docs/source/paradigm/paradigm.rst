.. _intro-label:

The Crowd Network Idea
======================

``Seidr`` is a product of an idea presented in the DREAM 5 network challenge 
[Marbach2012]_. In it, the authors show that gene regulatory network 
inference algorithms tend to suffer from biases towards specific interaction patters. 
They suggested a way to get around this by creating an aggregate of all the 
methods used in the study: a crowd network.

While the paper is widely cited, there is little software that attempts to
integrate the findings. ``Seidr`` is an attempt to create a toolbox that 
simplifies the laborious effort of creating crowd networks.

The basic pipeline
^^^^^^^^^^^^^^^^^^

A typical run of ``seidr`` has three steps:

* Infer: In the inference step, independent gene-gene networks are created by a multitude of algorithms.
* Import: In order to merge these networks, they are first sorted and ranked. To achieve this ``seidr`` uses its own file format: SeidrFiles (see :ref:`seidrfile-label`).
* Aggregate: Once all methods are ready, ``seidr`` can aggregate them to a crowd network.

.. image:: pipeline.png

In principle, any network can be input into ``seidr``, as long as it was constructed
under similar assumptions as all other networks. For example, it would be a bad idea
to take a subset of genes and create a network, which is then aggregated with another
subset using *different* genes. ``Seidr`` provides a number of algorithms as native
applications written in C++:

+----------------------+--------------------------+--------------------------+------------------+---------------------------+----------------+----------------+
| Name                 | Published                | Type                     | Orig. Lang.      | Seidr Lang.               | Orig. Parallel | Seidr Parallel |
+======================+==========================+==========================+==================+===========================+================+================+
| ANOVERENCE           | [Kueffner2012]_          | ANOVA                    | C++              | C++                       | No             | No             |
+----------------------+--------------------------+--------------------------+------------------+---------------------------+----------------+----------------+
| ARACNE               | [Margolin2006]_          | MI + DPI                 | C++              | C++                       | Yes            | Yes            |
+----------------------+--------------------------+--------------------------+------------------+---------------------------+----------------+----------------+
| CLR                  | [Faith2007]_ [Daub2004]_ | MI + CLR                 | MATLAB / C / C++ | C++                       | No             | Yes            |
+----------------------+--------------------------+--------------------------+------------------+---------------------------+----------------+----------------+
| Elastic Net ensemble | [Ruyssinck2014]_         | Elastic Net Regression   | R (glmnet)       | C++ (glmnet)              | No             | Yes            |
+----------------------+--------------------------+--------------------------+------------------+---------------------------+----------------+----------------+
| GENIE3               | [Huynh-Thu2010]_         | Random Forest Regression | R (randomForest) | C++ (ranger)              | No             | Yes            |
+----------------------+--------------------------+--------------------------+------------------+---------------------------+----------------+----------------+
| NARROMI              | [Zhang2013]_             | MI + Linear Programming  | MATLAB           | C++ (glpk)                | No             | Yes            |
+----------------------+--------------------------+--------------------------+------------------+---------------------------+----------------+----------------+
| Partial Correlation  | [Schafer2005]_           | Correlation              | R                | C++                       | No             | No             |
+----------------------+--------------------------+--------------------------+------------------+---------------------------+----------------+----------------+
| Pearson Correlation  | NA                       | Correlation              | NA               | C++                       | No             | No             |
+----------------------+--------------------------+--------------------------+------------------+---------------------------+----------------+----------------+
| PLSNET               | [Guo2016]_               | PLS                      | MATLAB           | C++                       | No             | Yes            |
+----------------------+--------------------------+--------------------------+------------------+---------------------------+----------------+----------------+
| Spearman Correlation | NA                       | Correlation              | NA               | C++                       | No             | No             |
+----------------------+--------------------------+--------------------------+------------------+---------------------------+----------------+----------------+
| SVM ensemble         | [Ruyssinck2014]_         | SVM regression           | R (libsvm) / C   | C++ (libsvm or liblinear) | No             | Yes            |
+----------------------+--------------------------+--------------------------+------------------+---------------------------+----------------+----------------+
| TIGRESS              | [Haury2012]_             | LASSO Regression         | MATLAB / R       | C++ (glmnet)              | No             | Yes            |
+----------------------+--------------------------+--------------------------+------------------+---------------------------+----------------+----------------+

Downstream
^^^^^^^^^^

Once you have a network, you probably want to explore it. To that end we provide 
some utilities to investigate the networks and to prepare them for input into
other software.


References
^^^^^^^^^^

.. [Marbach2012] Marbach, D., Costello, J. C., Küffner, R., Vega, N. M., Prill, R. J., Camacho, D. M., … Stolovitzky, G. (2012). Wisdom of crowds for robust gene network inference. Nature Methods, 9(8), 796–804. https://doi.org/10.1038/nmeth.2016
.. [Kueffner2012] Küffner, R., Petri, T., Tavakkolkhah, P., Windhager, L., & Zimmer, R. (2012). Inferring gene regulatory networks by ANOVA. Bioinformatics, 28(10), 1376–1382. https://doi.org/10.1093/bioinformatics/bts143
.. [Margolin2006] Margolin, A. A., Nemenman, I., Basso, K., Wiggins, C., Stolovitzky, G., Dalla Favera, R., & Califano, A. (2006). ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. BMC Bioinformatics, 7 Suppl 1, S7. https://doi.org/10.1186/1471-2105-7-S1-S7
.. [Faith2007] Faith, J. J., Hayete, B., Thaden, J. T., Mogno, I., Wierzbowski, J., Cottarel, G., … Gardner, T. S. (2007). Large-scale mapping and validation of Escherichia coli transcriptional regulation from a compendium of expression profiles. PLoS Biology, 5(1), e8. https://doi.org/10.1371/journal.pbio.0050008
.. [Daub2004] Daub, C. O., Steuer, R., Selbig, J., & Kloska, S. (2004). Estimating mutual information using B-spline functions--an improved similarity measure for analysing gene expression data. BMC Bioinformatics, 5, 118. https://doi.org/10.1186/1471-2105-5-118
.. [Ruyssinck2014] Ruyssinck, J., Huynh-Thu, V. A., Geurts, P., Dhaene, T., Demeester, P., & Saeys, Y. (2014). NIMEFI: Gene regulatory network inference using multiple ensemble feature importance algorithms. PLoS ONE, 9(3). https://doi.org/10.1371/journal.pone.0092709
.. [Huynh-Thu2010] Huynh-Thu, V. A., Irrthum, A., Wehenkel, L., & Geurts, P. (2010). Inferring regulatory networks from expression data using tree-based methods. PLoS ONE, 5(9), 1–10. https://doi.org/10.1371/journal.pone.0012776
.. [Zhang2013] Zhang, X., Liu, K., Liu, Z. P., Duval, B., Richer, J. M., Zhao, X. M., … Chen, L. (2013). NARROMI: A noise and redundancy reduction technique improves accuracy of gene regulatory network inference. Bioinformatics, 29(1), 106–113. https://doi.org/10.1093/bioinformatics/bts619
.. [Guo2016] Guo, S., Jiang, Q., Chen, L., & Guo, D. (2016). Gene regulatory network inference using PLS-based methods. BMC Bioinformatics, 17(1), 545. https://doi.org/10.1186/s12859-016-1398-6
.. [Schafer2005] Schäfer, J., & Strimmer, K. (2005). A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics. Statistical Applications in Genetics and Molecular Biology, 4(1). https://doi.org/10.2202/1544-6115.1175
.. [Haury2012] Haury, A.-C., Mordelet, F., Vera-Licona, P., & Vert, J.-P. (2012). TIGRESS: Trustful Inference of Gene REgulation using Stability Selection. BMC Systems Biology, 6(1), 145. https://doi.org/10.1186/1752-0509-6-145