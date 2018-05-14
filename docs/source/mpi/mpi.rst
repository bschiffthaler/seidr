.. _mpirun-label:

Using multiple processors to infer networks
===========================================

A number of computationally intensive network inference algorithms in ``seidr``
are written within a message parsing framework (OpenMPI 3). This allows for
parallelism on a single computer or across many nodes in a cluster. Some inference
algorithms in ``seidr`` have been run on hundreds of CPUs across many nodes in
a high performance compute cluster.

Running in MPI mode
^^^^^^^^^^^^^^^^^^^

By default all inference algorithms will use a single core to process data. Let's
use ``CLR`` as an example::

  mi -m CLR -i expr_mat.tsv -g genes.txt

This will spawn a single compute thread to process the data. In order to use 
multiple processors, we can run the same command as a child of the ``mpirun``
program. The ``-np`` option in ``mpirun`` defines the number of CPU cores.
In this case, I will choose 8 CPU cores::

  mpirun -np 8 mi -m CLR -i expr_mat.tsv -g genes.txt

This will spawn 7 compute threads and a single master thread. The master thread
will handle communication with all other threads and final data processing of the
result matrix and will require more memory than the compute threads. In a 
multi-node run, the master process (rank 0) should be on high memory node.

As the master thread is not very active during most of the computation, it doesn't
make sense to reserve an entire compute core for it, but ``mpirun`` won't let
us use more threads than compute cores by default. We can get around that by
using the ``--oversubscribe`` option to ``mpirun``. I have 8 threads available in
my computer, therefore I want to use 8 compute threads and 1 master thread::

  mpirun -np 9 --oversubscribe mi -m CLR -i expr_mat.tsv -g genes.txt

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
to :math:`\frac{n_{genes}}{n_{cpus} * 10}`. As an example, if I am estimating the
network for 25,000 genes using 64 compute threads, I set ``--batch-size`` to :math:`\frac{25000}{64 * 10} = 39`.