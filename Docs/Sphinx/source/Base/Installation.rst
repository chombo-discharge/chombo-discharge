.. _Chap:Installation:

Installation
============

Obtaining ``chombo-discharge``
------------------------------

``chombo-discharge`` can be freely obtained from `<https://github.com/chombo-discharge/chombo-discharge>`_.
The following packages are *required*:

* ``Chombo``, which is supplied with ``chombo-discharge``.
* The C++ JSON file parser `<https://github.com/nlohmann/json>`_.
* The ``EBGeometry`` package, see `<https://github.com/rmrsk/EBGeometry>`_.
* LAPACK and BLAS

The ``Chombo``, ``nlohmann/json``, and ``EBGeometry`` dependencies are automatically handled by ``chombo-discharge`` through git submodules.

.. warning::
   Our version of ``Chombo`` is hosted at `<https://github.com/chombo-discharge/Chombo-3.3.git>`_. 
   ``chombo-discharge`` has made substantial changes to the embedded boundary generation in ``Chombo``.
   It will not compile with other versions of ``Chombo`` than the one above.  

Optional packages are

* A serial or parallel version of HDF5, which is used for writing plot and checkpoint files.
* An MPI installation, which is used for parallelization.
* `<https://visit-dav.github.io/visit-website/>`_ visualization, which used for visualization.


Cloning ``chombo-discharge``
----------------------------

When compiling ``chombo-discharge``, makefiles must be able to find both ``chombo-discharge`` and ``Chombo``.
In our makefiles the paths to these are supplied through the environment variables

* ``DISCHARGE_HOME``, pointing to the ``chombo-discharge`` root directory.
* ``CHOMBO_HOME``, pointing to your ``Chombo`` library  

Note that ``DISCHARGE_HOME`` must point to the root folder in the ``chombo-discharge`` source code, while ``CHOMBO_HOME`` must point to the :file:`lib/` folder in your ``Chombo`` root directory.
When cloning with submodules, both ``Chombo`` and ``nlohmann/json`` will be placed in the :file:`Submodules` folder in ``$DISCHARGE_HOME``.  

.. tip::
   
   To clone ``chombo-discharge`` directly to ``$DISCHARGE_HOME``, set the environment variables and clone (using ``--recursive`` to fetch submodules):

   .. code-block:: bash

      export DISCHARGE_HOME=/home/foo/chombo-discharge
      export CHOMBO_HOME=$DISCHARGE_HOME/Submodules/Chombo-3.3/lib
		
      git clone --recursive git@github.com:chombo-discharge/chombo-discharge.git ${DISCHARGE_HOME}

   Alternatively, if cloning using https:

   .. code-block:: bash

      export DISCHARGE_HOME=/home/foo/chombo-discharge
      export CHOMBO_HOME=$DISCHARGE_HOME/Submodules/Chombo-3.3/lib
		
      git clone --recursive https://github.com/chombo-discharge/chombo-discharge.git ${DISCHARGE_HOME}   

``chombo-discharge`` is built using a configuration file supplied to ``Chombo``.
This file must reside in ``$CHOMBO_HOME/mk``.
Some standard configuration files are supplied with ``chombo-discharge``, and reside in ``$DISCHARGE_HOME/Lib/Local``.
These files may or may not work right off the bat. 

Test build
----------

For a quick compilation test the user can use the GNU configuration file supplied with ``chombo-discharge`` by following the steps below.

#. Copy the GNU configuration file

   .. code-block:: text

      cp $DISCHARGE_HOME/Lib/Local/Make.defs.GNU $CHOMBO_HOME/mk/Make.defs.local

#. If you do not have the GNU compiler suite, install it by

   .. code-block::
   
      sudo apt install csh gfortran g++ libblas-dev liblapack-dev

   This will install
   
      * LAPACK and BLAS
      * GNU compilers for Fortran and C++

#. Compile ``chombo-discharge`` 

   .. code-block:: text

      cd $DISCHARGE_HOME
      make -s -j4

This will compile the ``chombo-discharge`` source code in serial and without HDF5 (using four cores for the compilation).
If successful, ``chombo-discharge`` libraries will appear in ``$DISCARGE_HOME/Lib``.

.. _Chap:AdvancedConfig:

Full configuration
------------------

``chombo-discharge`` is compiled using GNU Make, following the ``Chombo`` configuration methods.

.. important::

   Compilers, libraries, and configuration options are defined in a file ``Make.defs.local`` which resides in ``$CHOMBO_HOME/mk``.
   Users need to supply this file in order to compile ``chombo-discharge``.
   
Typically, a full configuration consists of specifying

* Fortran and C++ compilers
* Specifying configurations. E.g., serial or parallel builds, and compiler flags. 
* Library paths (in particular for HDF5).

Main settings
_____________

The main variables that the user needs to set are

* ``DIM = 2/3`` The dimensionality (must be 2 or 3). 
* ``DEBUG = TRUE/FALSE``
  This enables or disables debugging flags and code checks/assertions.
* ``OPT = FALSE/TRUE/HIGH``.
  Setting ``OPT=TRUE/HIGH`` enables optimization flags that will speed up ``Chombo`` and ``chombo-discharge``.
* ``PRECISION = DOUBLE``
  Currently, ``chombo-discharge`` has not been wetted with single precision.
  Many algorithms (like conjugate gradient) depend on the use of double precision.
* ``CXX = <C++ compiler>``
* ``FC = <Fortran compiler>``
* ``MPI = TRUE/FALSE``
  This enables/disables MPI.
* ``MPICXX = <MPI compiler>`` This sets the MPI compiler.
* ``CXXSTD = 14`` For specifying the C++ standard. We are currently at C++14.
  Sets the C++ standard - we are currently at C++14.
* ``USE_EB=TRUE``
  Configures ``Chombo`` with embedded boundary functionality.
  This is a requirement. 
* ``USE_MF=TRUE``
  Configures ``Chombo`` with multifluid functionality.
  This is a requirement.
* ``USE_HDF5 = TRUE/FALSE``
  This enables and disables HDF5 code.
* ``OPENMPCC = TRUE/FALSE``
  Turn on/off OpenMP threading. 
  

MPI
___

To enable MPI, make sure that ``MPI`` is set to true and that the ``MPICXX`` compiler is set.
For GNU installations, one will usually have ``MPICXX = mpicxx`` or ``MPICXX = mpic++``, while for Intel builds one will usually have ``MPICXX = mpiicpc``.

.. note::

   The MPI layer distributes grid patches among processes, i.e. uses *domain decomposition*.

HDF5
____

If using HDF5, one must also set the following flags:

* ``HDFINCFLAGS      = -I<path to hdf5-serial>/include`` (for serial HDF5). 
* ``HDFLIBFLAGS      = -L<path to hdf5-serial>/lib -lhdf5 -lz`` (for serial HDF5)
* ``HDFMPIINCFLAGS   = -I<path to hdf5-parallel>/include`` (for parallel HDF5)
* ``HDFMPILIBFLAGS   = -L<path to hdf5-parallel>/lib -lhdf5 -lz`` (for parallel HDF5).

.. warning::

   ``Chombo`` only supports HDF5 APIs at version 1.10 and below.
   To use a newer version of HDF5 together with the 1.10 API, add ``-DH5_USE_110_API`` to the HDFINC flags.

OpenMP
______

To turn on OpenMP threading one can set the ``OPENMPCC`` to ``TRUE``.

.. note::

   The OpenMP parallelization layer uses threaded loops over grid patches (potentially within each MPI rank).


Compiler flags
______________

Compiler flags are set through

* ``cxxoptflags  = <C++ compiler flags``
* ``foptflags    = <Fortran compiler flags``
* ``syslibflags  = <system library flags>``

Note that LAPACK and BLAS are requirements in ``chombo-discharge``.
Linking to these can often be done using

* ``syslibflag = -llapack -lblas`` (for GNU compilers)
* ``syslibflag = -mkl=sequential`` (for Intel compilers)  
  

Pre-defined configuration files
_______________________________

Some commonly used configuration files are found in ``$DISCHARGE_HOME/Lib/Local``.
``chombo-discharge`` can be compiled in serial or with MPI, and with or without HDF5.
The user need to configure the ``Chombo`` makefile to ensure that the ``chombo-discharge`` is properly configured.
Below, we include brief instructions for compilation on a Linux workstation and for a cluster. 


GNU configuration for workstations
__________________________________

Here, we provide a more complete installation example using GNU compilers for a workstation.
These steps are intended for users that do not have MPI or HDF5 installed.
If you already have installed MPI and/or HDF5, the steps below might require modifications.

#. Ensure that ``$DISCHARGE_HOME`` and ``$CHOMBO_HOME`` point to the correct locations:

   .. code-block:: bash
		   
      echo $DISCHARGE_HOME
      echo $CHOMBO_HOME

#. Install GNU compiler dependencies by

   .. code-block::
   
      sudo apt install csh gfortran g++ libblas-dev liblapack-dev

   This will install

      * LAPACK and BLAS
      * GNU compilers for Fortran and C++   

#. To also install OpenMPI and HDF5:

   .. code-block::

      sudo apt install libhdf5-dev libhdf5-openmpi-dev openmpi-bin

   This will install

      * OpenMPI
      * HDF5, both serial and parallel.

   Both serial and parallel HDF5 will be installed, and these are *usually* found in folders

     * ``/usr/lib/x86_64-linux-gnu/hdf5/serial/`` for serial HDF5
     * ``/usr/lib/x86_64-linux-gnu/hdf5/openmpi/`` for parallel HDF5 (using OpenMPI). 
     
   Before proceeding further, the user need to locate the HDF5 libraries (if building with HDF5). 

#. After installing the dependencies, copy the desired configuration file to ``$CHOMBO_HOME/mk``:

   * **Serial build without HDF5**:

     .. code-block:: text

	cp $DISCHARGE_HOME/Lib/Local/Make.defs.GNU $CHOMBO_HOME/mk/Make.defs.local

   * **Serial build with HDF5**:

     .. code-block:: text

	cp $DISCHARGE_HOME/Lib/Local/Make.defs.HDF5.GNU $CHOMBO_HOME/mk/Make.defs.local

   * **MPI build without HDF5**:

     .. code-block:: text

	cp $DISCHARGE_HOME/Lib/Local/Make.defs.MPI.GNU $CHOMBO_HOME/mk/Make.defs.local

   * **MPI build with HDF5**:

     .. code-block:: text

	cp $DISCHARGE_HOME/Lib/Local/Make.defs.MPI.HDF5.GNU $CHOMBO_HOME/mk/Make.defs.local

#. Compile the ``chombo-discharge``

   .. code-block:: text

      cd $DISCHARGE_HOME
      make -s -j4 discharge-lib

This will compile the ``chombo-discharge`` source code using the configuration settings set by the user.
To compile ``chombo-discharge`` in 3D, do ``make -s -j4 DIM=3 discharge-lib``.
If successful, ``chombo-discharge`` libraries will appear in ``$DISCARGE_HOME/Lib``.

Configuration on clusters
_________________________

To configure ``chombo-discharge`` for executation on a cluster, use one of the makefiles supplied in ``$DISCHARGE_HOME/Lib/Local`` if it exists for your computer.
Alternatively, copy ``$DISCHARGE_HOME/Lib/Local/Make.defs.local.template`` to ``$CHOMBO_HOME/mk/Make.defs.local`` and set the compilers, optimization flags, and paths to HDF5 library.

On clusters, MPI and HDF5 are usually already installed, but must usually be loaded (e.g. as modules) before compilation.

Configuration files for GitHub
______________________________

``chombo-discharge`` uses GitHub actions for continuous integration and testing.
These tests run on Linux for a selection of GNU and Intel compilers.
The configuration files are located in :file:`$DISCHARGE_HOME/Lib/Local/GitHub`. 

Running example applications
----------------------------

In ``chombo-discharge``, applications are set up so that they use the ``chombo-discharge`` source code and one ``chombo-discharge`` physics module.
These are normally set up through Python interfaces accompanying each module. 
Several example applications are given in :file:`$DISCHARGE_HOME/Exec`, which are organized by example type (e.g., plasma simulation, electrostatics, radiative transfer, etc).
If ``chombo-discharge`` built successfully, it will usually be sufficient to compile the example by navigating to the folder containing the program file (:file:`program.cpp`) and compiling it:

.. code-block:: text

   make -s -j4 program

To see how these programs are run, see :ref:`Chap:Control`.   

Positive streamer in air
________________________

To run one of the applications that use a particular ``chombo-discharge`` physics module, we will run a simulation of a positive streamer (in air). 

The application code is located in ``$DISCHARGE_HOME/Exec/Examples/CdrPlasma/DeterministicAir`` and it uses the convection-diffusion-reaction plasma module (located in ``$DISCHARGE_HOME/Physics/CdrPlasma``).

First, compile the application by

.. code-block:: text

   cd $DISCHARGE_HOME/Exec/Examples/CdrPlasma/DeterministicAir
   make -s -j4 DIM=2 program

This will provide an executable named ``program2d.<bunch_of_options>.ex``.
If one compiles for 3D, i.e. ``DIM=3``, the executable will be named ``program3d.<bunch_of_options>.ex``.

To run the application do:

**Serial build**

.. code-block:: text

   ./program2d.<bunch_of_options>.ex positive2d.inputs

**Parallel build**
  
.. code-block:: text

   mpirun -np 8 program2d.<bunch_of_options>.ex positive2d.inputs   

If the user also compiled with HDF5, plot files will appear in the subfolder ``plt``.

.. tip::

   One can track the simulation progress through the :file:`pout.*` files, see :ref:`Chap:pout`.

.. _Chap:TroubleShooting:

Troubleshooting
---------------

If the prerequisites are in place, compilation of ``chombo-discharge`` is usually straightforward.
However, due to dependencies on ``Chombo`` and HDF5, compilation can sometimes be an issue.
Our experience is that if ``Chombo`` compiles, so does ``chombo-discharge``.

If experiencing issues, try cleaning ``chombo-discharge`` by

.. code-block:: bash

   cd $DISCHARGE_HOME
   make pristine

.. note::

   Do not hesitate to contact us at `GitHub <https://github.com/chombo-discharge/chombo-discharge>`_ regarding installation issues.

Recommended configurations
__________________________

Production runs
^^^^^^^^^^^^^^^

For production runs, we generally recommend that the user compiles with ``DEBUG=FALSE`` and ``OPT=HIGH``.
These settings can be set directly in :file:`Make.defs.local`.
Alternatively, they can be included directly on the command line when compiling problems.

Debugging
^^^^^^^^^

If you believe that there might be a bug in the code, one can compile with ``DEBUG=TRUE`` and ``OPT=TRUE``.
This will turn on some assertions throughout ``Chombo`` and ``chombo-discharge``.    

Common problems
_______________

* Missing library paths:

   On some installations the linker can not find the HDF5 library.
   To troubleshoot, make sure that the the environment variable ``LD_LIBRARY_PATH`` can find the HDF5 libraries:

   .. code-block:: bash

      echo $LD_LIBRARY_PATH

   If the path is not included, it can be defined by:

   .. code-block:: bash

      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/<path_to_hdf5_installation>/lib

* Incomplete perl installations.

  ``Chombo`` may occasionally complain about incomplete perl modules.
  These error messages are unrelated to ``Chombo`` and ``chombo-discharge``, but the user may need to install additional perl modules before compiling ``chombo-discharge``.

  
