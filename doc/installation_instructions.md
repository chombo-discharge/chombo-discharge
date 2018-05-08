Installation
============

Prerequisites
----------------------

The numerical code is written in C++ using the [MFChombo](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations) library and can be compiled on UNIX systems (both workstations and clusters). 
Compilation requires the following prerequisites to be in place:

1. A C++ compiler (e.g. g++, icpc)
2. A Fortran compiler (e.g. gfortran, ifort)
3. An MPI compiler (e.g. OpenMPI)
4. A working HDF5 installations, both parallel and serial. 

Compilers can be installed straightforwardly on UNIX systems; for users of Ubuntu-like systems we recommend that, as a first step, you install g++ (and gcc), gfortran, and OpenMPI. In general, compilation of Chombo requires a pre-compiled, parallel version of HDF5, and modification of Chombo's most important configuration file Make.defs.local which resides in mf-chomb/lib/mk/Make.defs.local. This file contains the compiler configuration (e.g. Intel or GNU compilers), dimensionality, optimization flags etc. This file is system-dependent and must be adapted to your own system.

For clusters on which we have successfully run the executables, the process of building and running is more straightforward, so you may find the configuration file attached below (see @ref cluster-installations). 

* [Local installations](installation_local.md)
* [Cluster installations](@ref cluster-installations)
