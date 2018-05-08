Compilation on fram.sigma2.no {#installation-fram}
=============================

TLDR
----
Load the modules by

    module load intel/2017a

The configuration file for fram.sigma2.no can be found here. Store it as <chombo_root>/lib/mk/Make.defs.local

Loading the modules
-------------------
To compile the system on the *fram* supercomputer, you need to load compilers and libraries before compiling the code. The intel compiler toochain is loaded by

    module load intel/2017a

which will load Intel compilers for Fortran and C/C++, as well as Intel MPI.

Modifying the configuration file
--------------------------------
The modify the configuration file, find the HDF5 path by

    module disp HDF5/1.10.1-intel-2017a

which will display the library paths to HDF5.

To use Intel's Fortran and C/C++ compilers, please specify

     MPICXX = mpiicpc
     FC     = ifort
     CXX    = icpc

in Make.defs.local. Likewise, update Make.defs.local so that it finds the HDF5 libraries and include files. 

We also recommend that you compile with -xCORE-AVX2 (put it in cxxoptflags and foptflags) to enable advanced vector extensions. 