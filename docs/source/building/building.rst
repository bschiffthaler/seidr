.. _building-label:

Building Seidr
==============

Supported OSs
^^^^^^^^^^^^^

``seidr`` should build fine on most Linux distributions.

Test builds of ``seidr`` are created on Ubuntu 18.04 and Fedora 31. It is possible to
build on Mac OS X (with some effort). Microsoft Windows is currently not supported.

Basic Build
^^^^^^^^^^^

Currently, ``seidr`` has the following dependencies (examplary dnf packages on Fedora):

* ``gcc`` 
* ``gcc-c++`` 
* ``gcc-gfortran`` 
* ``cmake`` 
* ``git``
* ``boost-devel`` 
* ``glpk-devel`` or COIN-OR CLP (see :ref:`clp-label`) 
* ``armadillo-devel``
* ``zlib-devel``

Once the dependencies are satisfied, build with::

  git clone --recursive https://github.com/bschiffthaler/seidr
  cd seidr
  mkdir build
  cmake -DCMAKE_BUILD_TYPE=Release ..
  make

If you have multiple CPU cores, run make as ``make -j <ncpus>`` to speed up building.

Building with MPI
^^^^^^^^^^^^^^^^^

If you have access to a compute cluster, you might want to build ``seidr`` with 
MPI support. If you have the MPI libraries installed (e.g.: ``openmpi-devel`` on Fedora) add::

  cmake -DSEIDR_WITH_MPI=ON ..

to the CMake build options.

Building Parallel STL (PSTL)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have Intel TBB and PSTL availble, you can build ``seidr`` with support for
parallel STL algorithms, which can speed some operations. To do that, add::

  cmake -DSEIDR_PSTL=ON ..

to the CMake build options.

.. _clp-label:

A note on CLP and GLPK
^^^^^^^^^^^^^^^^^^^^^^

The ``narromi`` algorithm uses linear programming routines, which in ``seidr`` is
implemented via either GLPK or CLP backends. GLPK is widely available, but not
safe to use in an OpenMP context, you will therefore be limited to a single OMP thread.
CLP is safer, but packages are less widely available (you might need to build form source).
If you want to build ``seidr`` with the CLP backend add::

  cmake -DNARROMI_USE_CLP=ON ..

