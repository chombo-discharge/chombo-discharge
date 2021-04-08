Installation on local machines {#installation-local}
==============================


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