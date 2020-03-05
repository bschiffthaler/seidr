.. _mpirun-label:

Using multiple processors to infer networks
===========================================

A number of computationally intensive network inference algorithms in ``seidr``
are written using a hybrid MPI/OpenMP approach. This allows for sahred memory
parallelism on a single computer or across many nodes in a cluster. Some inference
algorithms in ``seidr`` have been run on hundreds of CPUs across many nodes in
a high performance compute cluster.

Running in OMP mode
^^^^^^^^^^^^^^^^^^^

By default, if your computer has multiple CPU cores availble, ``seidr`` will use
as many as it can. If the subprogram has parallel processing support, you can
control the extent of the parallelization with the ``-O,--threads`` option.

Example::
  # Use all available threads by default:
  seidr import ...

  # Use two threds
  seidr import -O 2 ...

  # Use environment variables to control the number of threads
  export OMP_NUM_THREADS=2
  seidr import ..

Running in MPI mode
^^^^^^^^^^^^^^^^^^^

By default all inference algorithms will use all cores to process data. Let's
use ``CLR`` as an example::

  mi -m CLR -i expr_mat.tsv -g genes.txt

This will spawn eight compute threads (on my laptop) to process the data.
In order to control the allocated number of CPUs, we can use the ``-O`` flag
of the ``mi`` program::

  mi -O 4 -m CLR -i expr_mat.tsv -g genes.txt

This will use 4 compute threads.

If we want to use multiple nodes, we can use  we can run the same command as a child of the ``mpirun``
program. You should first define a `hostfile <https://www.open-mpi.org/doc/current/man1/mpirun.1.php#sect6>`_.::

  mpirun -hostfile myhostfile.cfg mi -m CLR -i expr_mat.tsv -g genes.txt

This will spawn a distributed versiob of the MI inference, running the maximum
amount of OpenMP threads. You can combine ``mpirun`` and the program's ``-O``
argument to control the number of compute threads each MPI worker spawns.

A special note on MPI rank order: the highest memory node on the cluster you are
using should always be rank 0. If there are any high memory tasks, Seidr will
assign them to this MPI worker.

For more info on running MPI jobs (including running them on several nodes), please
refer to the `OpenMPI webpage <https://www.open-mpi.org/faq/?category=running>`_

.. _batchsize-label:

The batchsize argument
^^^^^^^^^^^^^^^^^^^^^^

All MPI enabled inference algorithms in ``seidr`` have a ``--batch-size`` argument.
This is the number of genes a compute thread will process at once before requesting
more from the master thread. Lower batch sizes will lead to more time spent in I/O
operations and more temporary files, but setting it too high might leave compute
threads without work for portions of the run. A good rule of thumb is to set this
to :math:`\frac{n_{genes}}{n_{nodes}}`. As an example, if I am estimating the
network for 25,000 genes using a five nodes, I set ``--batch-size`` to :math:`\frac{25000}{5} = 5000`. In general, it is safe to let ``seidr`` decide on the batch size.