.. _Chap:GettingStarted:

Getting started
===============

This chapter discusses how you may obtain `PlasmaC` and compile it.

.. _Chap:Obtaining:

Obtaining `PlasmaC`
---------------------

`PlasmaC` is obtained by cloning the following repository:

.. code-block:: bash

   git clone ssh://git@git.code.sintef.no/~robertm/plasmac



.. _Chap:Prerequisites:

Prerequisites
-------------

From the ground up, `PlasmaC` is built on top of the `Chombo <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>`_ framework.
To compile `PlasmaC`, you must first install the following:

* A Fortran compiler, usually gfortran or Intel Fortran
* A C++ compiler, usually g++ or Intel C++
* An MPI installation
* A parallel HDF5 installation
* LAPACK and BLAS

Usually, laptops and desktops already have appropriate Fortran and C++ compilers installed, as well as a version of MPI.
On clusters, HDF5 is (usually) preinstalled, and in this case, it will be sufficient to modify the `Chombo` build files in order to compile `PlasmaC`.
Following some changes that we've made to the way `Chombo` generates its embedded boundaries, `Chombo` is released together with `PlasmaC`. 

If you already have HDF5 installed, you may skip directly to :ref:`Chap:Environment`. 

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

In `Chombo`,the system information is supplied through a file known as Make.defs.local, which resides in the `Chombo` library itself.
This file contains a number of build settings, such as dimension, compilers, paths to HDF5 and so on.
The configuration file is /lib/mk/Make.defs.local in your `Chombo` folder.
The build parameters can also be controlled through the command line. 

Here are the configuration variables that have been used used on the `fram <https://www.top500.org/system/179072>`_ supercomputer

.. code-block:: c++

		DIM           = 3
		DEBUG         = FALSE
		OPT           = HIGH
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
		cxxoptflags    = -Ofast -xCORE-AVX2 -march=native
		#cxxprofflags  =
		#fcppflags     =
		#fdbgflags     =
		foptflags      = -Ofast -xCORE-AVX2
		#fprofflags    =
		flibflags      = -lblas -llapack
		#lddbgflags    =
		#ldoptflags    =
		#ldprofflags   =
		syslibflags    = -ldl -lm -lz


We also recommend that you create environment variables that hold the path to your `PlasmaC` version. For example,

.. code-block:: c++

		PLASMAC_HOME=/home/foo/plasmac

This environment variables is used in the `PlasmaC` makefile system so that our makefiles can find `PlasmaC` and `Chombo`.

.. _Chap:Compiling:

Compiling `PlasmaC`
---------------------

In `PlasmaC`, each problem is compiled as a mini-application into a subfolder.
Mini-apps are usually set up through a Python pre-compilation script that generates the required source code, makefiles, and simulation parameters.
There is no separate build for the `PlasmaC` source code and your own application files and you will *not* be able to install `PlasmaC` as a separate library.

Once an application has been set up, the makefile system tracks the necessary `Chombo` and `PlasmaC` source files.
Compiling is done in the subfolder that houses your mini-app:

.. code-block:: bash

   make -s -j8 DIM=2 OPT=HIGH <application_name>

We generally recommend that you compile with ``OPT=HIGH`` for performance reasons. 

.. _Chap:Visualization:

Visualization
-------------

PlasmaC writes output files to HDF5. Users can decide what data to output, as well as restrict plot depth to a certain grid levels level. There are also options for including ghost cells in the output files.

Our favorite tool for visualization is `VisIt <https://wci.llnl.gov/codes/visit/>`_, which can be freely downloaded. Our experience is that client-server visualization is beneficial for visualization of three-dimensional simulation data. For information on how to set up host profiles for VisIt, please contact your local guru or refer to the `VisIt documentation <http://visit-sphinx-user-manual.readthedocs.io/en/latest/index.html>`_. 

.. _Chap:MyFirstCompilation:

My first compilation
--------------------

Before moving on with more complex descriptions of `PlasmaC`, we will try to compile a test problem which simply advects a scalar.
The application we will use is a part of the regression testing in `PlasmaC`.

To run this application, navigate to :file:`/regression/advection_diffusion` and compile with

.. code-block:: bash

   make -s -j4 DIM=2 main

where ``-j4`` is the number of cores used for the compilation. If you want to compile this example in 3D, you should put DIM=3.
If the application compiles successfully, you will see a file called :file:`main2d.<bunch_of_options>.ex`.
If you see this file, you will be able to compile all of `PlasmaC`. If you don't, you won't be able to compile any of it.
Before moving on further, please make sure that your model compiles.

Once we have compiled our application, we are ready to run it.
The example that we will run is a very simple setup of scalar advection and diffusion of a rotating flow, where the base code is provided in :file:`/physics/advection_diffusion`.
To run the example, you can do

.. code-block:: bash

   mpirun -np 4 main2d.<bunch_of_options>.ex regression2d.inputs

Output files should now appear in :file:`/regression/advection_diffusion/plt`. 

Troubleshooting
---------------

If the prerequisites are in place, compilation of `PlasmaC` is usually straightforward.
However, due to dependencies on `Chombo` and HDF5, compilation can sometimes be an issue.
Our experience is that if `Chombo` compiles, so does `PlasmaC`.
For that reason we refer you to the `Chombo` user guide for troubleshooting.

Using this documentation
------------------------

This documentation was built using `reStructuredText` with `Sphinx`. If you want to build a PDF version of this documentation, please navigate to :file:`plasmac/doc/sphinx` and execute

.. code-block:: bash
   
   make latexpdf

A PDF version of this documentation named :file:`PlasmaC.pdf` will appear in :file:`plasma/doc/sphinx/build/latex`. 
