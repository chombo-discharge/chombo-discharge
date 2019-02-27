.. _Chap:GettingStarted:

Getting started
===============

This chapter discusses what `PlasmaC` really is, how you may obtain `PlasmaC` and compile it.

.. _Chap:Obtaining:

Obtaining `PlasmaC`
---------------------

`PlasmaC` is obtained by cloning the following repository:

.. code-block:: bash

   git clone ssh://git@git.code.sintef.no/~robertm/plasmac

You will also need a version of Chombo, which you may download online.

.. _Chap:Prerequisites:

Prerequisites
-------------

From the ground up, `PlasmaC` is built on top of the `Chombo <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>`_ framework. To compile `PlasmaC`, you must therefore have the following in place:

* A Fortran compiler, usually gfortran or Intel Fortran
* A C++ compiler, usually g++ or Intel C++
* An MPI installation
* A parallel HDF5 installation
* A Chombo library

Usually, laptops and desktops already have appropriate Fortran and C++ compilers installed, as well as a version of MPI. On clusters, HDF5 is (usually) preinstalled, and in this case, it will be sufficient to modify the Chombo build files in order to compile `PlasmaC`. If you already have HDF5 installed, you may skip directly to :ref:`Chap:Environment`.

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

   .. _Chap:Environment:

Setting up your environment
___________________________

In Chombo,the system information is supplied through a file known as Make.defs.local, which resides in the Chombo library itself. This file contains a number of build settings, such as dimension, compilers, paths to HDF5 and so on. The file itself is in Chombo/lib/mk/Make.defs.local. The build parameters can also be controlled through the command line. 

Here are what configuration variables that are used on the ``fram`` supercomputer

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


We also recommend that you create environment variables that hold the path to your Chombo and `PlasmaC` libraries. For example,

.. code-block:: c++

		CHOMBO_HOME=/usr/local/Chombo-3.2
		PLASMAC_HOME=/home/foo/plasmac

These two environment variables are used in the `PlasmaC` makefile system so that our makefiles can find both Chombo and `PlasmaC`.

.. _Chap:Compiling:

Compiling `PlasmaC`
---------------------

Once an application has been set up, compiling is done by

.. code-block:: bash

   make -s -j 16 DIM=2 <application_name>

Compiling must be performed from the folder which houses your makefile. 

Currently, all of `PlasmaC` is compiled into your mini-applications. While this is something that we are working on improving, this means that there is no separate build for the `PlasmaC` source code and your application files.

.. _Chap:Visualization:

Visualization
-------------

PlasmaC writes it's output files to HDF5. In general, user's can decide what data to output, as well as restrict plot depth to a certain AMR level. There are also options for including ghost cells in the output files.

Our favorite tool for visualization is `VisIt <https://wci.llnl.gov/codes/visit/>`_, which can be freely downloaded. Our experience is that client-server visualization is beneficial, especially for in three dimensions. For information on how to set up host profiles for VisIt, please contact your local guru or refer to the `VisIt documentation <http://visit-sphinx-user-manual.readthedocs.io/en/latest/index.html>`_. 
