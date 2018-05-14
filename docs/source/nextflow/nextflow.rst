.. _nextflow-label:

Running seidr with nextflow
===========================

If you are looking for a convenient way to run seidr, you can consider running it off 
`nextflow <https://www.nextflow.io/>`_ . We provide an example configuration 
(``nextflow.config``) and nextflow pipeline ``vala.nf`` in the ``nextflow``
directory of the project root.

Currently, running on a local machine, and on a SLURM cluster are supported. Modify
the relevant entries in ``nextflow.config`` and then run::

  nextflow run vala.nf

