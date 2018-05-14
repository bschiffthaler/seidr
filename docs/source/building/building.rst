.. _building-label:

Building Seidr
==============

Currently, seidr has the following dependencies (examplary dnf packages on Fedora):

* gcc 
* gcc-c++ 
* gcc-gfortran 
* cmake 
* git
* openmpi-devel 
* tclap 
* boost-devel 
* glpk-devel 
* armadillo-devel 
* zlib-devel
* gsl-devel 

It also needs `htslib <https://github.com/samtools/htslib>`_

Once the dependencies are satisfied, build with::

  git clone --recursive https://github.com/bschiffthaler/seidr
  cd seidr
  mkdir build
  cmake -DCMAKE_BUILD_TYPE=Release ..
  make

If you have multiple CPU cores, run make as ``make -j <ncpus>`` to speed up building.