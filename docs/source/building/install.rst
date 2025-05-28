.. _install-label:

Installing Seidr
================

Conda
^^^^^

The easiest way to install seidr is via ``conda``. The ``conda`` package comes in three flavours. One without MPI support, and one package each with support for OpenMPI and MPICH respectively::

  # No MPI support (run on one node/computer)
  conda create -n seidr -c bioconda -c conda-forge 'seidr=*=nompi*'

  # OpenMPI (multi-node)
  conda create -n seidr -c bioconda -c conda-forge 'seidr=*=mpi_openmpi*'
  # MPICH
  conda create -n seidr -c bioconda -c conda-forge 'seidr=*=mpi_mpich*'

Building from source
^^^^^^^^^^^^^^^^^^^^

See :ref:`building-label`.