chombo-discharge
----------------

This is chombo-discharge, a multiphysics code which uses Chombo for plasma simulations with adaptive mesh refinement (AMR) on embedded boundary grids.

A modified version of Chombo is distributed together with this code.
chombo-discharge only uses Chombo; it is not affiliated nor endorsed by LBNL.

Installation instructions are found below, and also included in the full user documentation which is available at https://chombo-discharge.github.io.

See LICENSE and Copyright.txt for redistribution rights. 

Dependencies
------------

chombo-discharge is dependent on the following packages:

* Chombo (supplied with chombo-discharge)
* LAPACK and BLAS (required). 
* HDF5 (optional, used for writing plot and checkpoint file).
* MPI (optional, used for parallelization).
* VisIt visualization (optional, but frequently used for visualization). 

Obtaining chombo-discharge
--------------------------

To clone chombo-discharge, set the environment variable ``$DISCHARGE_HOME`` as follows

.. code-block:: text
		
   export DISCHARGE_HOME=<Location for chombo-discharge>

There are two ways of cloning ``chombo-discharge``: 

Using submodules
________________

``Chombo`` can be included as a submodule in ``chombo-discharge``.
To clone the ``chombo-discharge`` code and ``Chombo``, do

.. code-block:: text
		
		git clone --recursive git@github.com:chombo-discharge/chombo-discharge.git ${DISCHARGE_HOME}

Next, set the ``Chombo`` environment variable ``$CHOMBO_HOME`` to ``$DISCHARGE_HOME/Chombo-3.3/lib``, i.e.

.. code-block:: text

		export CHOMBO_HOME=$DISCHARGE_HOME/Chombo-3.3/lib

Maintaining  ``Chombo`` separately.
___________________________________

First clone ``chombo-discharge`` *without* submodules by

.. code-block:: text
		
		git clone git@github.com:chombo-discharge/chombo-discharge.git ${DISCHARGE_HOME}

Next, set the ``$CHOMBO_HOME`` environment variable and clone ``Chombo`` there, i.e.

.. code-block:: text

		export CHOMBO_HOME=<Location for Chombo>
		git clone git@github.com:chombo-discharge/Chombo-3-3.git ${CHOMBO_HOME}

Configuration
-------------

``chombo-discharge`` is compiled using Make, following the ``Chombo`` makefiles.
Compilers, libraries, and configuration options are defined in a file ``Make.defs.local`` which resides in ``$CHOMBO_HOME/mk``.
Users need to supply this file in order to compile ``chombo-discharge``.
Typically, these steps include

* Specifying Fortran and C++ compilers
* Specifying configurations. E.g., serial or parallel builds, and compiler flags. 
* Specifying library paths (in particular for HDF5).

Some commonly used configuration files are found in :file:`$DISCHARGE_HOME/Local`. 

Test build
----------

chombo-discharge can be compiled in serial or with MPI, and with or without HDF5.
The user need to configure the Chombo makefile to ensure that the chombo-discharge is properly configured.

An example configuration process for workstations is included below, but for a quick test the user can use the GNU configuration file as follows:

.. code-block:: text

   cp $DISCHARGE_HOME/Local/Make.defs.GNU $CHOMBO_HOME/lib/mk/Make.defs.local

Next, compile ``chombo-discharge`` by

.. code-block:: text

   cd $DISCHARGE_HOME
   make -s -j4 lib

This will compile the ``chombo-discharge`` source code in serial and without HDF5 (using four cores for the compilation).
If successful, ``chombo-discharge`` libraries will appear in :file:`$DISCARGE_HOME/Lib`.

GNU configuration for workstations
----------------------------------

Here, we provide an example installation process using GNU compilers for a workstation.

First, one can currently install the GNU compiler suite by

.. code-block::
   
   sudo apt install csh gfortran g++ libblas-dev liblapack-dev

To also install OpenMPI and HDF5:

.. code-block::

   sudo apt install libhdf5-dev libhdf5-openmpi-dev openmpi-bin

This will install

* LAPACK and BLAS
* GNU compilers for Fortran and C++
* OpenMPI
* HDF5, both serial and parallel.
  
Equivalent steps for Intel compilers will differ slightly.

Both serial and parallel HDF5 will be installed, and these are *usually* found in folders ``/usr/lib/x86_64-linux-gnu/serial/`` and ``/usr/lib/x86_64-linux-gnu/parallel/``.
Before proceeding further, the user should ensure that he can locate both the serial and parallel HDF5 libraries.

After installing the dependencies, copy the desired configuration file to ``$CHOMBO_HOME/lib/mk``:

* For serial build without HDF5:

  .. code-block:: text

     cp $DISCHARGE_HOME/Local/Make.defs.GNU $CHOMBO_HOME/lib/mk/Make.defs.local

* For serial build with HDF5:

  .. code-block:: text

     cp $DISCHARGE_HOME/Local/Make.defs.HDF5.GNU $CHOMBO_HOME/lib/mk/Make.defs.local

* For MPI build without HDF5:

  .. code-block:: text

     cp $DISCHARGE_HOME/Local/Make.defs.MPI.GNU $CHOMBO_HOME/lib/mk/Make.defs.local

* For MPI build with HDF5:

  .. code-block:: text

     cp $DISCHARGE_HOME/Local/Make.defs.MPI.HDF5.GNU $CHOMBO_HOME/lib/mk/Make.defs.local               

After that, compile one of the applications by

.. code-block:: text

   cd $DISCHARGE_HOME/Regression/AdvectionDiffusion/Godunov
   make -s -j4 main

To run the application, do

.. code-block:: text

   ./main2d.<bunch_of_options>.ex regression2d.inputs

if using a serial build, and

.. code-block:: text

   mpirun -np 4 main2d.<bunch_of_options>.ex regression2d.inputs

if using a parallel build.
If compiling with HDF5, plot files will then appear in ``$DISCHARGE_HOME/Regression/AdvectionDiffusion/Godunov/plt``. 

Configuration on clusters
-------------------------

To configure chombo-discharge for executation on a cluster, use one of the makefiles supplied in ``$DISCHARGE_HOME/Local`` if it exists for your computer.
Alternatively, copy ``$DISCHARGE_HOME/Local/Make.defs.local.template`` to ``$CHOMBO_HOME/lib/mk/Make.defs.local`` and set the compilers, optimization flags, and paths to HDF5 library.

On clusters, MPI and HDF5 are usually already installed, but must usually be loaded (e.g. as modules) before compilation.

Troubleshooting
---------------

Compilation is normally straightforward, but if experiencing problems, try cleaning ``Chombo`` by

.. code-block:: text

   cd $CHOMBO_HOME
   make realclean

Likewise, when compiling applications, compile with ``make clean`` rather than just ``make``.
More tips and tricks are given in the documentation at https://chombo-discharge.github.io. 

Contributing
------------
We welcome feedback, bug reports, or code contributions. Use the github issue tracker and pull request system for code contributions
See code documentation for coding style and review system. 


