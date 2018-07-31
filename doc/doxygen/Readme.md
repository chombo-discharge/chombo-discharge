Installation
============

Building the documentation
==========================

Before proceeding further, please install doxygen and build the code documentation. This requires that doxygen is installed on your computer, but it does not require a working installation of the code. To build the documentation, please type

    doxygen doxyfile 

in your shell (from the main folder). This will compile the documentation of the code in which you will find installation instructions. After building the documentation, please proceed ../chombo-streamer/doc/doxygen/html/index.html and navigate to installation instructions. 


Installation instructions
=========================

Prerequisites
----------------------

The numerical code is written in C++ using the [MFChombo](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations) library and can be compiled on UNIX systems (both workstations and clusters). 
Compilation requires the following prerequisites to be in place:

1. A C++ compiler (e.g. g++, icpc)
2. A Fortran compiler (e.g. gfortran, ifort)
3. An MPI compiler (e.g. OpenMPI)
4. A working HDF5 installations, both parallel and serial. 

Compilers can be installed straightforwardly on UNIX systems; for users of Ubuntu-like systems we recommend that you install g++ (and gcc), gfortran, and OpenMPI. 


HDF5 installation
-----------------

Before attempting to install HDF5, please make sure that C++, Fortran, and MPI compilers are in place. 

HDF5 (Hierarchical Data Format) is a portable, binary data file with parallel support and is used to reduce I/O bottlenecks in simulations. HDF5 can be obtained from [here](https://support.hdfgroup.org/HDF5/). 

While it is possible to compile our code for serial execution only, we generally recommend that the user installs both serial and parallel versions of HDF5. The two versions should *not* be installed in the same directory.
For versions of HDF5 newer than 1.6, we require that the user configures HDF5 with the 1.6 API version.

1. Compile and install zlib (a compression library). zlib can be obtained through

        sudo apt-get install zlib1g-dev

2. Download and install HDF5 for *serial execution*. 

        ./configure --prefix=/usr/local/hdf5 --enable-production --enable-cxx --enable-fortran --with-default-api-version=v16
        make
        make install

 This will install a serial version of HDF5 in /usr/local/hdf5 (if you want to install 

3. Install HDF5 for *parallel execution*

        ./configure --prefix=/usr/local/phdf5 --enable-production --enable-fortran --with-default-api-version=v16 --enable-parallel
        make
        make install

 This will install a parallel version of HDF5 in /usr/local/phdf5


Building the source code
------------------------

Once sequential and parallel versions of HDF5 are installed, the source code may be built. The source is built by inheriting a set of makefiles from Chombo. By default, a version of Chombo-3.2 is shipped with this code.
To build the source code, you must ensure that the makefile finds the Chombo library by modifying the GNUmakefile which lies at the top directory in the code. Modify this file so that `CHOMBO_HOME` points to the correct location (a full path specification may be necessary). 
For example, if the code is in the home folder of a user `exampleUser`, `CHOMBO_HOME` should be

    CHOMBO_HOME := /home/exampleUser/FluidAMRStreamer/Chombo-3.2/lib

Specification of default settings for Chombo are necessary for compiling both the Chombo library and the source code. You should modify Make.defs.local which lies in Chombo-3.2/lib/mk/ for specification of compilers and path to the serial and parallel HDF5 installations. If you installed HDF5 in /usr/local/hdf5 and /usr/local/phdf5, the default settings in Make.defs.local are (probably) sufficient. 

The streamer code can be compiled for two-dimensional and three-dimensional executation on Cartesian grids (default is 2D). To compile the code, it should be sufficient to navigate to one of the plasma models and then execute

    make "model"

where "model" is replaced by the model name. If you want to silence GNUmake during the compilation stage, include a flag -s. I.e.

    make -s "model"

This will inherit a set of makefiles: Chombo source codes are compiled first and the streamer code is compiled second. To execute for three-dimensional execution, either make appropriate changes to Make.defs.local (DIM=3), or run

    make DIM=3 "model"

By default, the executable will be named according to the specifications in Make.defs.local. For example, if the code was compiled for optimized three-dimensional execution using MPI (mpiCC as compiler) gfortran on a Linux system, the executable is named

    plasma3d.Linux.64.mpiCC.gfortran.OPT.MPI.ex

The code has compiled successfully if this file is created. 

Compilation on vilje
====================================

To compile the system on the *vilje* supercomputer, you need to load compilers and libraries before compiling the code. We recommend that you use intel compilers by loading

    module load intelcomp/15.0.1

which will load Intel compilers for Fortran and C/C++ (version 15.0.1). MPI libraries are loaded by 

    module load mpt/2.10

and HDF5 is loaded by 

    module load hdf5/1.8.14

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

Compilation on fram
====================================
To compile the system on the *fram* supercomputer, you need to load compilers and libraries before compiling the code. The intel compiler toochain is loaded by

    module load intel/2017a

which will load Intel compilers for Fortran and C/C++, as well as Intel MPI. The code may be compiled as usual. It will be necessary to update CHOMBO_HOME and supply the library path to the Chombo makefiles. 
CHOMBO_HOME is simply updated to the folder in which Chombo is installed. HDF5 paths are found by typing

    module disp HDF5/1.10.1-intel-2017a

which will display the library paths to HDF5.

To use Intel's Fortran and C/C++ compilers, please specify

     MPICXX = mpiicpc
     FC     = ifort
     CXX    = icpc

in Make.defs.local. Likewise, update Make.defs.local so that it finds the HDF5 libraries and include files. 

We also recommend that you compile with -xCORE-AVX2 (put it in cxxoptflags and foptflags) to enable advanced vector extensions. 
