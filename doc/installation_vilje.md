Compilation on vilje.hpc.ntnu.no {#installation-vilje}
================================

TLDR
----
Load modules by

    module load intelcomp/15.0.1
    module load gcc/4.9.1
    module load mpt/2.10

The configuration file for vilje is found here. Store it as <chombo_root>/lib/mk/Make.defs.local

Loading the modules
-------------------
To compile the system on the *vilje* supercomputer, you need to load compilers and libraries before compiling the code. We recommend that you use intel compilers by loading

    module load intelcomp/15.0.1

which will load Intel compilers for Fortran and C/C++ (version 15.0.1). An updated version of gcc is loaded by

    module load gcc/4.9.1
    
MPI libraries are loaded by 

    module load mpt/2.10

and HDF5 is loaded by 

    module load hdf5/1.8.14

Modifying the configuration file
--------------------------------

It is likely possible to load newer versions of compilers, MPI, and HDF5 (but we haven't tried). Once the libraries are loaded, the code may be compiled as usual. It will be necessary to update CHOMBO_HOME and supply the library path to the Chombo makefiles. 
CHOMBO_HOME is simply updated to the folder in which Chombo is installed. HDF5 paths are found by typing

    module disp hdf5/1.8.14

which will display the library paths to HDF5. The typical output is

    /sw/sdev/modulefiles/hdf5/1.8.14:
    
    module-whatis    general purpose library and file format for storing scientific data 
    prereq   intelcomp/15.0.1 
    prereq   mpt/2.10 
    prepend-path     PATH /sw/sdev/Modules/hdf5/hdf5-1.8.14/bin 
    prepend-path     CPATH /sw/sdev/Modules/hdf5/hdf5-1.8.14/include 
    prepend-path     FPATH /sw/sdev/Modules/hdf5/hdf5-1.8.14/include 
    prepend-path     LIBRARY_PATH /sw/sdev/Modules/hdf5/hdf5-1.8.14/lib 
    prepend-path     LD_LIBRARY_PATH /sw/sdev/Modules/hdf5/hdf5-1.8.14/lib 
    prepend-path     PYTHONPATH /sw/sdev/Modules/hdf5/hdf5-1.8.14/lib/python2.7/site-packages 
    setenv           HDF5_DIR /sw/sdev/Modules/hdf5/hdf5-1.8.14 
        
    Loads HDF5 version 1.8.14

To use Intel's Fortran and C/C++ compilers, please specify

     FC  = ifort
     CXX = icpc

in Make.defs.local. Likewise, update Make.defs.local so that it finds the HDF5 libraries and include files. For this version of HDF5, they are located in

     /sw/sdev/Modules/hdf5/hdf5-1.8.14/include 
     /sw/sdev/Modules/hdf5/hdf5-1.8.14/lib 

With these paths, the appropriate flags in Make.defs.local are

    HDFINCFLAGS    = -I/sw/sdev/Modules/hdf5/hdf5-1.8.14/include 
    HDFLIBFLAGS    = -L/sw/sdev/Modules/hdf5/hdf5-1.8.14/lib -lhdf5 -lz 

    HDFMPIINCFLAGS = -I/sw/sdev/Modules/hdf5/hdf5-1.8.14/include 
    HDFMPILIBFLAGS = -L/sw/sdev/Modules/hdf5/hdf5-1.8.14/lib -lhdf5 -lz