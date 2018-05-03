.. _aracne-label:

ARACNE
==========

ARACNE ([Margolin2006]_) is an inference algorithm based on mutual information
and applies data processing inequality to delete most indirect edges.

Our implementation differs to the original in that it estimates the initial
mutual information using a B-spline approach as described in [Daub2004]_ .

Running ARACNE
^^^^^^^^^^^^^^^^^^



References
^^^^^^^^^^

.. [Margolin2006] Margolin, A. A., Nemenman, I., Basso, K., Wiggins, C., Stolovitzky, G., Dalla Favera, R., & Califano, A. (2006). ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. BMC Bioinformatics, 7 Suppl 1, S7. https://doi.org/10.1186/1471-2105-7-S1-S7
.. [Daub2004] Daub, C. O., Steuer, R., Selbig, J., & Kloska, S. (2004). Estimating mutual information using B-spline functions--an improved similarity measure for analysing gene expression data. BMC Bioinformatics, 5, 118. https://doi.org/10.1186/1471-2105-5-118
