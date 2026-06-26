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
* VisIt (`<https://visit-dav.github.io/visit-website>`_), which used for visualization and analysis.


Cloning ``chombo-discharge``
----------------------------

``chombo-discharge`` is compiled using GNUmake.
When compiling ``chombo-discharge``, the makefiles must be able to find both ``chombo-discharge`` and ``Chombo``.
In our makefiles the paths to these are supplied through the environment variables

* ``DISCHARGE_HOME``, pointing to the ``chombo-discharge`` root directory.
* ``CHOMBO_HOME``, pointing to your ``Chombo`` library.

Note that ``DISCHARGE_HOME`` must point to the root folder in the ``chombo-discharge`` source code, while ``CHOMBO_HOME`` must point to the :file:`lib/` folder in your ``Chombo`` root directory.
When cloning recursively with submodules, both ``Chombo`` and ``nlohmann/json`` will be placed in the :file:`Submodules` folder in ``$DISCHARGE_HOME``.

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

See :ref:`Chap:BuildConfig` below for how to supply ``Make.defs.local`` before building.

.. _Chap:BuildConfig:

Build configuration
-------------------

``chombo-discharge`` inherits the ``Chombo`` build system, which reads compiler settings,
library paths, and build flags from a file called ``Make.defs.local``.
This file must exist at:

.. code-block:: text

   $DISCHARGE_HOME/Submodules/Chombo-3.3/lib/mk/Make.defs.local

For convenience, ``$DISCHARGE_HOME`` contains a symlink:

.. code-block:: text

   $DISCHARGE_HOME/Make.defs.local -> Submodules/Chombo-3.3/lib/mk/Make.defs.local

You can therefore create and edit ``Make.defs.local`` directly from the repository root
without navigating into the ``Chombo`` submodule.
The file is gitignored and never committed — each machine keeps its own copy.

Pre-made configuration files
____________________________

Several ready-to-use configuration files are provided in ``$DISCHARGE_HOME/Lib/Local``:

.. list-table::
   :header-rows: 1
   :widths: 40 20 10 10

   * - File
     - Compilers
     - MPI
     - HDF5
   * - ``Make.defs.GNU``
     - GCC
     - No
     - No
   * - ``Make.defs.HDF5.GNU``
     - GCC
     - No
     - Yes
   * - ``Make.defs.MPI.GNU``
     - GCC
     - Yes
     - No
   * - ``Make.defs.MPI.HDF5.GNU``
     - GCC
     - Yes
     - Yes

Copy the one that matches your setup via the symlink:

.. code-block:: bash

   cp $DISCHARGE_HOME/Lib/Local/Make.defs.GNU $DISCHARGE_HOME/Make.defs.local

For a fully custom configuration, start from the annotated template instead:

.. code-block:: bash

   cp $DISCHARGE_HOME/Lib/Local/Make.defs.local.template $DISCHARGE_HOME/Make.defs.local

Then open ``$DISCHARGE_HOME/Make.defs.local`` and uncomment and set the variables you need.
See :ref:`Chap:AdvancedConfig` for a description of the available variables.

Quick start on a Linux workstation
___________________________________

The steps below install all required dependencies and compile ``chombo-discharge`` in 2D
without MPI or HDF5.

#. Install dependencies:

   .. code-block:: bash

      sudo apt install csh gfortran g++ libblas-dev liblapack-dev

#. Optionally, install OpenMPI and HDF5:

   .. code-block:: bash

      sudo apt install libhdf5-dev libhdf5-openmpi-dev openmpi-bin

   Serial and parallel HDF5 are normally installed in:

   * ``/usr/lib/x86_64-linux-gnu/hdf5/serial/`` (serial)
   * ``/usr/lib/x86_64-linux-gnu/hdf5/openmpi/`` (parallel, OpenMPI)

#. Copy the desired configuration file:

   * **Serial, no HDF5**:

     .. code-block:: bash

        cp $DISCHARGE_HOME/Lib/Local/Make.defs.GNU $DISCHARGE_HOME/Make.defs.local

   * **Serial with HDF5**:

     .. code-block:: bash

        cp $DISCHARGE_HOME/Lib/Local/Make.defs.HDF5.GNU $DISCHARGE_HOME/Make.defs.local

   * **MPI, no HDF5**:

     .. code-block:: bash

        cp $DISCHARGE_HOME/Lib/Local/Make.defs.MPI.GNU $DISCHARGE_HOME/Make.defs.local

   * **MPI with HDF5**:

     .. code-block:: bash

        cp $DISCHARGE_HOME/Lib/Local/Make.defs.MPI.HDF5.GNU $DISCHARGE_HOME/Make.defs.local

#. Compile:

   .. code-block:: bash

      cd $DISCHARGE_HOME
      make -s -j4 discharge-lib

   For a 3D build: ``make -s -j4 DIM=3 discharge-lib``.

If successful, ``chombo-discharge`` libraries will appear in ``$DISCHARGE_HOME/Lib``.

.. _Chap:AdvancedConfig:

Full configuration
------------------

``chombo-discharge`` is compiled using GNU Make, following the ``Chombo`` configuration methods.

.. important::

   Compilers, libraries, and configuration options are defined in ``Make.defs.local``.
   See :ref:`Chap:BuildConfig` for where this file lives and how to create it.

.. _Chap:MainSettings:

Main settings
_____________

The main variables that the user needs to set are

* ``DIM = 2/3`` — number of spatial dimensions (must be 2 or 3).
* ``DEBUG = TRUE/FALSE`` — enables or disables debugging flags and assertions.
  ``chombo-discharge`` runs substantially slower with ``DEBUG=TRUE``.
* ``OPT = FALSE/TRUE/HIGH`` — optimization level.
* ``CXX = <C++ compiler>``
* ``FC = <Fortran compiler>``
* ``MPI = TRUE/FALSE`` — enables or disables MPI.
* ``MPICXX = <MPI C++ compiler>``
* ``CXXSTD = 17`` — C++ standard. ``chombo-discharge`` requires at least C++17.
* ``USE_EB = TRUE`` — enables embedded boundary support. This is required.
* ``USE_MF = TRUE`` — enables multifluid support. This is required.
* ``USE_MT = TRUE/FALSE`` — enables ``Chombo`` memory tracking.
  Not compatible with OpenMP.
* ``USE_HDF = TRUE/FALSE`` — enables HDF5 output.
* ``OPENMPCC = TRUE/FALSE`` — enables OpenMP threading.
* ``USE_PETSC = TRUE/FALSE`` — links to PETSc.
* ``PRECISION = DOUBLE`` — floating-point precision.
  Single precision is not supported.


MPI
___

To enable MPI, set ``MPI = TRUE`` and set the ``MPICXX`` compiler.
For GNU installations this is typically ``mpicxx`` or ``mpic++``; for Intel builds it is usually ``mpiicpc``.

.. note::

   The MPI layer distributes grid patches among processes (domain decomposition).

HDF5
____

If using HDF5, also set:

* ``HDFINCFLAGS      = -I<path>/include`` (serial HDF5)
* ``HDFLIBFLAGS      = -L<path>/lib -lhdf5 -lz`` (serial HDF5)
* ``HDFMPIINCFLAGS   = -I<path>/include`` (parallel HDF5)
* ``HDFMPILIBFLAGS   = -L<path>/lib -lhdf5 -lz`` (parallel HDF5)

.. warning::

   ``Chombo`` only supports HDF5 APIs at version 1.10 and below.
   To use a newer version with the 1.10 API, add ``-DH5_USE_110_API`` to the HDF include flags.

OpenMP
______

Set ``OPENMPCC = TRUE`` to enable OpenMP threading.
When compiled with OpenMP, all loops over grid patches use:

.. code-block:: c++

   #pragma omp parallel for schedule(runtime)
   for (int mybox = 0; mybox < nbox; mybox++) {

   }

.. warning::

   Memory tracking is not supported together with OpenMP.
   When using OpenMP, set ``USE_MT = FALSE``.

.. warning::

   OpenMP support (including hybrid OpenMP+MPI) is **experimental** and not yet validated for
   production. Several threaded code paths -- including some inside ``Chombo`` itself, which are
   only activated when it is compiled with ``OPENMPCC = TRUE`` -- contain data races that can
   cause intermittent heap corruption and crashes at higher thread counts. For production runs,
   use pure MPI (``OPENMPCC = FALSE``); OpenMP and OpenMP+MPI are intended for experimentation only.

PETSC
_____

Set ``USE_PETSC = TRUE`` to link against PETSc.
The user must first install PETSc and set ``PETSC_DIR`` and ``PETSC_ARCH``.


Compiler flags
______________

Compiler flags are set through:

* ``cxxoptflags  = <C++ compiler flags>``
* ``foptflags    = <Fortran compiler flags>``
* ``syslibflags  = <system library flags>``

LAPACK and BLAS are required. Typical values:

* ``syslibflags = -llapack -lblas`` (GNU compilers)
* ``syslibflags = -mkl=sequential`` (Intel compilers)

The ``cxxoptflags`` and ``foptflags`` variables apply to optimized builds (``OPT=TRUE/HIGH``).
Corresponding debug flags are ``cxxdbgflags`` and ``fdbgflags``.

The provided GNU configurations include ``-fno-math-errno`` in ``cxxoptflags``; this lets the compiler auto-vectorize loops that call ``sqrt``/``pow`` (it only drops the ``errno`` side-effect and does not change numerical results). The Intel configurations get the same effect via ``-Ofast``.

Configuration on clusters
_________________________

Use one of the provided makefiles in ``$DISCHARGE_HOME/Lib/Local`` if one exists for your cluster
(see ``$DISCHARGE_HOME/Lib/Local/betzy``, ``fram``, etc.).
Otherwise, copy ``Make.defs.local.template`` and set the compilers, optimization flags, and HDF5 paths.

On clusters, MPI and HDF5 are usually pre-installed but must be loaded as modules before compilation.

.. _Chap:TroubleShooting:

Troubleshooting
---------------

If the prerequisites are in place, compilation of ``chombo-discharge`` is usually straightforward.
However, due to dependencies on ``Chombo`` and HDF5, compilation can sometimes be an issue.
Our experience is that if ``Chombo`` compiles, so does ``chombo-discharge``.

If experiencing issues, try remove the ``chombo-discharge`` installation first by running

.. code-block:: bash

   cd $DISCHARGE_HOME
   make pristine

.. note::

   Do not hesitate to contact us at `GitHub <https://github.com/chombo-discharge/chombo-discharge>`_ regarding installation issues.

Recommended configurations
__________________________

Production runs
^^^^^^^^^^^^^^^

For production runs, compile with ``DEBUG=FALSE`` and ``OPT=HIGH``.
These can be set in :file:`Make.defs.local` or passed directly on the command line.

Debugging
^^^^^^^^^

To enable assertions throughout ``Chombo`` and ``chombo-discharge``, compile with ``DEBUG=TRUE`` and ``OPT=TRUE``.

Common problems
_______________

* Missing library paths:

   On some installations the linker cannot find the HDF5 library.
   To troubleshoot, check that ``LD_LIBRARY_PATH`` includes the HDF5 library path:

   .. code-block:: bash

      echo $LD_LIBRARY_PATH

   If not, add it:

   .. code-block:: bash

      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/<path_to_hdf5_installation>/lib

* Incomplete perl installations.

  ``Chombo`` may occasionally complain about incomplete perl modules.
  These error messages are unrelated to ``Chombo`` and ``chombo-discharge``, but additional perl modules may need to be installed before compiling.
