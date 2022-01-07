chombo-discharge
----------------

This is ``chombo-discharge``, a multiphysics code which uses ``Chombo`` for discharge simulations with adaptive mesh refinement (AMR) on embedded boundary grids.

A modified version of ``Chombo`` is distributed together with this code.
``chombo-discharge`` only uses ``Chombo``; it is not affiliated nor endorsed by LBNL.

Installation instructions are found below and in the user documentation which is available at https://chombo-discharge.github.io. 

See LICENSE and Copyright.txt for redistribution rights. 


Obtaining chombo-discharge
--------------------------

To clone ``chombo-discharge``, set the environment variable ``$DISCHARGE_HOME`` as follows

.. code-block:: text
		
   export DISCHARGE_HOME=<Location for chombo-discharge>

There are two ways of cloning ``chombo-discharge``, using submodules or maintaining ``Chombo`` separately.
The use of submodules is recommended.

Using submodules
________________

``Chombo`` can be included as a submodule in ``chombo-discharge``.
To clone both, do

.. code-block:: text
		
   git clone --recursive git@github.com:chombo-discharge/chombo-discharge.git ${DISCHARGE_HOME}

Next, set the ``Chombo`` environment variable ``$CHOMBO_HOME`` to ``$DISCHARGE_HOME/Chombo-3.3/lib``, i.e.

.. code-block:: text

   export CHOMBO_HOME=$DISCHARGE_HOME/Chombo-3.3/lib

Maintaining  ``Chombo`` separately.
___________________________________

When maintaining ``Chombo`` separately, first clone ``chombo-discharge`` *without* submodules by

.. code-block:: text
		
   git clone git@github.com:chombo-discharge/chombo-discharge.git ${DISCHARGE_HOME}

Next, set the ``$CHOMBO_HOME`` environment variable and clone ``Chombo`` there, i.e.

.. code-block:: text

   export CHOMBO_HOME=<Location for Chombo>
   git clone git@github.com:chombo-discharge/Chombo-3-3.git ${CHOMBO_HOME}

Installation
-------------

Dependencies
____________

``chombo-discharge`` is dependent on the following packages:

* ``Chombo`` (required, and supplied with ``chombo-discharge``)
* LAPACK and BLAS (required). 
* HDF5 (optional, used for writing plot and checkpoint file).
* MPI (optional, used for parallelization).
* VisIt visualization (optional, but frequently used for visualization). 

Configuration
_____________

``chombo-discharge`` is compiled using GNU Make, following the ``Chombo`` makefiles.
Compilers, libraries, and configuration options are defined in a file ``Make.defs.local`` which resides in ``$CHOMBO_HOME/mk``.
Users need to supply this file in order to compile ``chombo-discharge``.
Typically, these steps include

* Specifying Fortran and C++ compilers
* Specifying configurations. E.g., serial or parallel builds, and compiler flags. 
* Specifying library paths (in particular for HDF5).

Some commonly used configuration files are found in ``$DISCHARGE_HOME/Local``.

Before proceeding further, ensure that both ``$DISCHARGE_HOME`` and ``$CHOMBO_HOME`` are defined and point to correct locations.
E.g.,

.. code-block:: txt

   echo $DISCHARGE_HOME
   echo $CHOMBO_HOME

``chombo-discharge`` can be compiled in serial or with MPI, and with or without HDF5.
The user need to configure the ``Chombo`` makefile to ensure that the ``chombo-discharge`` is properly configured.

Below, we include brief instructions for compilation in serial without HDF5 (for testing purposes), full configuration with MPI and HDF5 on a workstation, and for a cluster. 

Test build
__________

For a quick test the user can use the GNU configuration file supplied with ``chombo-discharge`` by following the steps below:

#. Copy the GNU configuration file

   .. code-block:: text

      cp $DISCHARGE_HOME/Local/Make.defs.GNU $CHOMBO_HOME/lib/mk/Make.defs.local

#. Install dependencies

   .. code-block::
   
      sudo apt install csh gfortran g++ libblas-dev liblapack-dev

   This will install
   
   * LAPACK and BLAS
   * GNU compilers for Fortran and C++

#. Optionally, compile ``Chombo`` before compiling ``chombo-discharge``

   .. code-block:: text

      cd $DISCHARGE_HOME
      make -s -j4 chombo   

#. Compile ``chombo-discharge`` 

   .. code-block:: text

      cd $DISCHARGE_HOME
      make -s -j4 discharge-lib

This will compile the ``chombo-discharge`` source code in serial and without HDF5 (using four cores for the compilation).
If successful, ``chombo-discharge`` libraries will appear in ``$DISCARGE_HOME/Lib``.

GNU configuration for workstations
__________________________________

Here, we provide a more complete installation example using GNU compilers for a workstation.

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

   Both serial and parallel HDF5 will be installed, and these are *usually* found in folders ``/usr/lib/x86_64-linux-gnu/hdf5/serial/`` and ``/usr/lib/x86_64-linux-gnu/hdf5/parallel/``.
   Before proceeding further, the user need to locate the HDF5 libraries (if building with HDF5). 

#. After installing the dependencies, copy the desired configuration file to ``$CHOMBO_HOME/lib/mk``:

   * **Serial build without HDF5**:

     .. code-block:: text

	cp $DISCHARGE_HOME/Local/Make.defs.GNU $CHOMBO_HOME/lib/mk/Make.defs.local

   * **Serial build with HDF5**:

     .. code-block:: text

	cp $DISCHARGE_HOME/Local/Make.defs.HDF5.GNU $CHOMBO_HOME/lib/mk/Make.defs.local

   * **MPI build without HDF5**:

     .. code-block:: text

	cp $DISCHARGE_HOME/Local/Make.defs.MPI.GNU $CHOMBO_HOME/lib/mk/Make.defs.local

   * **MPI build with HDF5**:

     .. code-block:: text

	cp $DISCHARGE_HOME/Local/Make.defs.MPI.HDF5.GNU $CHOMBO_HOME/lib/mk/Make.defs.local

#. Optionally, compile ``Chombo`` before compiling ``chombo-discharge``

   .. code-block:: text

      cd $DISCHARGE_HOME
      make -s -j4 chombo

#. Compile the ``chombo-discharge``

   .. code-block:: text

      cd $DISCHARGE_HOME
      make -s -j4 discharge-lib

This will compile the ``chombo-discharge`` source code using the configuration settings set by the user.
To compile ``chombo-discharge`` in 3D, do ``make -s -j4 DIM=3 discharge-lib``.
If successful, ``chombo-discharge`` libraries will appear in ``$DISCARGE_HOME/Lib``.


Configuration on clusters
_________________________

To configure ``chombo-discharge`` for executation on a cluster, use one of the makefiles supplied in ``$DISCHARGE_HOME/Local`` if it exists for your computer.
Alternatively, copy ``$DISCHARGE_HOME/Local/Make.defs.local.template`` to ``$CHOMBO_HOME/lib/mk/Make.defs.local`` and set the compilers, optimization flags, and paths to HDF5 library.

On clusters, MPI and HDF5 are usually already installed, but must usually be loaded (e.g. as modules) before compilation.

Compiling physics modules
-------------------------

The ``chombo-discharge`` physics modules are maintained separately from the ``chombo-discharge`` source code. 
To compile the physics modules, navigate to ``$DISCHARGE_HOME`` and compile the physics modules by

.. code-block:: text

   cd $DISCHARGE_HOME
   make -s -j4 physics

This will compile all physics modules.
If successful, ``chombo-discharge`` libraries will appear in ``$DISCHARGE_HOME/Lib``. 

Running an example application
------------------------------

In ``chombo-discharge``, applications are set up so that they use the ``chombo-discharge`` source code and one ``chombo-discharge`` physics module. 
To run one of the applications that use a particular ``chombo-discharge`` physics module, we will run a simulation of a positive streamer (in air). 

The application code is located in ``$DISCHARGE_HOME/Examples/CdrPlasma/AirStreamer`` and it uses the convection-diffusion-reaction plasma module (located in ``$DISCHARGE_HOME/Physics/CdrPlasma``).

First, compile the application by

.. code-block:: text

   cd $DISCHARGE_HOME/Examples/CdrPlasma/AirStreamer
   make -s -j4 DIM=2 program

This will provide an executable named ``program2d.<bunch_of_options>.ex``.
If one compiles for 3D, i.e. ``DIM=3``, the executable will be named ``program2d.<bunch_of_options>.ex``.

To run the application do:

* **Serial build**

  .. code-block:: text

     ./program2d.<bunch_of_options>.ex positive2d.inputs

* **Parallel build**
  
  .. code-block:: text

     ./program2d.<bunch_of_options>.ex positive2d.inputs   

If the user also compiled with HDF5, plot files will appear in the subfolder ``plt``. 


Troubleshooting
---------------

Compilation is normally straightforward, but if experiencing problems, try cleaning ``Chombo`` and ``chombo-discharge`` by

.. code-block:: text

   cd $CHOMBO_HOME
   make realclean

   cd $DISCHARGE_HOME
   make pristine

Likewise, when compiling applications, compile with ``make clean <target>`` rather than just ``make <target>``.
More tips and tricks are given in the documentation at https://chombo-discharge.github.io. 

Contributing
------------
We welcome feedback, bug reports, or code contributions. Use the github issue tracker and pull request system for code contributions
See code documentation for coding style and review system. 
