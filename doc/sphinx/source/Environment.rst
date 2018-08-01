.. _Chap:Environment:

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

These two environment variables are used in the PlasmaC makefile system so that our makefiles can find Chombo and PlasmaC. Strictly speaking, you don't *have* to set these as environment variables. However, both variables are used in the mini-application makefiles so if you don't use environment variables, you will need to specify them directly in your makefile. 
