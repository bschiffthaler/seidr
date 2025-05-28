<p align=center>
<img src="https://raw.githubusercontent.com/bschiffthaler/seidr/devel/docs/resources/logo.svg" width=150px height=150px>
</p>

# Seiðr

Harness the wisdom of the crowd for gene network inference.

# Questions?

Talk to us on Gitter: [![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/seidr-networks/community)

# Build Status

* Master branch : [![CircleCI](https://circleci.com/gh/bschiffthaler/seidr/tree/master.svg?style=svg)](https://circleci.com/gh/bschiffthaler/seidr/tree/master)
* Devel branch : [![CircleCI](https://circleci.com/gh/bschiffthaler/seidr/tree/devel.svg?style=svg)](https://circleci.com/gh/bschiffthaler/seidr/tree/devel)

# Installation

## From Source

[https://seidr.readthedocs.io/en/latest/source/building/building.html](https://seidr.readthedocs.io/en/latest/source/building/building.html)

## Bioconda

### Non-MPI version

```bash
conda create -n seidr -c bioconda -c conda-forge 'seidr=*=nompi*'
```

### OpenMPI

```bash
conda create -n seidr -c bioconda -c conda-forge 'seidr=*=mpi_openmpi*'
```

### MPICH

```bash
conda create -n seidr -c bioconda -c conda-forge 'seidr=*=mpi_mpich*'
```

# Documentation

Please find the full documentation at [http://seidr.readthedocs.io/](http://seidr.readthedocs.io/)
