.. _Chap:GettingStarted:

Getting started
===============



To get started using PlasmaC, we recommend that you follow the :ref:`Chap:WorkedExample` chapter. This example will provide the entire workflow from generating a geometry, setting up your mini-app, and running it. For information on how to set up your environment so that PlasmaC can be compiled, see :ref:`Chap:Compiling`

If you are looking for a detailed description on the PlasmaC interface, you can find it in the :ref:`Chap:Interface` chapter.

If you are already using PlasmaC, and want to check out what various parameters in your input script do, check out the :ref:`Chap:ImportantClasses` chapter. 

.. _Chap:Compiling:

Compiling PlasmaC
-----------------

From the ground up, PlasmaC is built on top of the `Chombo <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>`_ framework. To compile PlasmaC, you must therefore have the following in place:

* A Fortran compiler, typically gfortran or Intel Fortran
* A C++ compiler, typically g++ or Intel C++
* An MPI installation
* A parallel HDF5 installation

Typically, local clients (laptops and desktops) already have appropriate Fortran and C++ compilers installed, as well as a version of MPI. On clusters, HDF5 is also preinstalled (usually), and in this case, it will be sufficient to modify the Chombo build files in order to compile PlasmaC. If you already have HDF5 installed, you may skip directly to :ref:`Chap:MakeDefsLocal`.

.. _Chap:HDF5:

Installing HDF5
_______________

If you do not have HDF5 installed, you may do the following:

1. Compile and install zlib, which is a compression library used by HDF5. zlib can be installed by
   
   .. code-block:: c++
		
		sudo apt-get install zlib1g-dev

2. Download HDF5 (version 1.8 or newer) and install it for parallel execution

      .. code-block:: c++
		
		      ./configure --prefix=/usr/local/hdf5 --enable-production --enable-fortran --enable-parallel
		      make
		      make install

   This will install HDF5 in /usr/local/hdf5. 
   
.. _Chap:MakeDefsLocal:

Setting up your environment
___________________________

In Chombo,the system information is supplied through a file known as Make.defs.local, which resides in the Chombo library itself. This file contains a number of build settings, such as dimension, compilers, paths to HDF5 and so on. The file itself is in Chombo/lib/mk/Make.defs.local. Strictly speaking, it is not necessary to modify this file since all build parameters can be controlled through the command line, but it usually pays off.

Here are what configuration variables that we use on the fram supercomputer

.. code-block:: c++

		DIM           = 3
		DEBUG         = FALSE
		OPT           = TRUE
		#PRECISION     =
		#PROFILE       =
		CXX            = icpc
		FC             = ifort
		MPI            = TRUE
		## Note: don't set the MPICXX variable if you don't have MPI installed
		MPICXX         = mpicxx
		#OBJMODEL      =
		#XTRACONFIG    =
		## Optional features
		#USE_64        =
		#USE_COMPLEX   =
		USE_EB         = TRUE
		#USE_CCSE      =
		USE_HDF        = TRUE
		HDFINCFLAGS    = -I/cluster/software/HDF5/1.10.1-intel-2017a/include
		HDFLIBFLAGS    = -L/cluster/software/HDF5/1.10.1-intel-2017a/lib -lhdf5 -lz
		## Note: don't set the HDFMPI* variables if you don't have parallel HDF installed
		HDFMPIINCFLAGS = -I/cluster/software/HDF5/1.10.1-intel-2017a/include
		HDFMPILIBFLAGS = -L/cluster/software/HDF5/1.10.1-intel-2017a/lib -lhdf5 -lz
		USE_MF         = TRUE
		#USE_MT        =
		#USE_SETVAL    =
		#CH_AR         =
		#CH_CPP        =
		#DOXYGEN       =
		#LD            =
		#PERL          =
		#RANLIB        =
		#cppdbgflags   =
		#cppoptflags   =
		#cxxcppflags   =
		#cxxdbgflags   =
		cxxoptflags    = -O2 -xCORE-AVX2
		#cxxprofflags  =
		#fcppflags     =
		#fdbgflags     =
		foptflags      = -O2 -xCORE-AVX2
		#fprofflags    =
		flibflags      = -lblas -llapack
		#lddbgflags    =
		#ldoptflags    =
		#ldprofflags   =
		syslibflags    = -ldl -lm -lz


We also recommend that you create environment variables that hold the path to your Chombo and PlasmaC libraries. For example,

.. code-block:: c++

		CHOMBO_HOME=/usr/local/Chombo-3.2
		PLASMAC_HOME=/home/foo/plasmac

These two environment variables are used in the PlasmaC makefile system so that our makefiles can find Chombo and PlasmaC. Strictly speaking, you don't HAVE to set these as environment variables. However, both variables are used in the mini-application makefiles so if you don't use environment variables, you will need to specify them directly in your makefile. 
