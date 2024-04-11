.. _Chap:Control:

Controlling ``chombo-discharge``
================================

In this chapter we give a brief overview of how to run a ``chombo-discharge`` simulation and control its behavior through input scripts or command line options.

Running ``chombo-discharge``
----------------------------

How one runs ``chombo-discharge`` depends on the type of parallelism one compiled with.
Below, we consider basic examples for serial and parallel execution.

Serial
______

If the application was compiled for serial execution one runs it with:

.. code-block:: bash

   ./<application_executable> <input_file>

where <input_file> is your input file.

Parallel with OpenMP
____________________

When running with OpenMP one must specify the number of threads, and possibly also the binding of threads.
``chombo-discharge`` is compiled with run-time thread scheduling (which defaults to static), which can be specified.
For example

.. code-block:: bash

   export OMP_NUM_THREADS=8
   export OMP_PLACES=cores
   export OMP_PROC_BIND=true
   export OMP_SCHEDULE="dynamic, 4"
   
   ./<application_executable> <input_file>


Parallel with MPI
_________________

If the executable was compiled with MPI, one executes with e.g. ``mpirun`` (or one of its aliases):

.. code-block:: bash
	     
   mpirun -np 8 <application_executable> <input_file>

On clusters, this is a little bit different and usually requires passing the above command through a batch system.
Normally, the MPI installation will map processes to cores.
With OpenMP one can use ``--report-bindings`` to verify the mapping.

Parallel with MPI+OpenMP
________________________

When running with both MPI and OpenMP the user must

#. Bind each MPI rank to a specified resource (e.g., a node, socket, or list of CPUs)
#. Bind OpenMP threads to resources available to each MPI rank.

For example, the following may not work as expected:

.. code-block:: bash

   export OMP_NUM_THREADS=4
   mpiexec -n 2 ./<application_executable> <input_file>

Each MPI rank may spawn threads on the same physical cores.
With e.g. OpenMPI one can map each rank to a specified number of CPUs, and bind threads to those CPUs.
For example, on a local workstation one might do

.. code-block:: bash

   export NRANKS=2		
   export OMP_NUM_THREADS=4
   export OMP_PLACES=cores
   export OMP_PROC_BIND=true
   
   mpiexec --bind-to core --map-by slot:PE=$OMP_NUM_THREADS -n $NRANKS ./<application_executable> <input_file>

MPI bindings are here bound to cores and will spawn threads only on the cores they are associated with.
It is often useful to verify this by reporting the bindings as follows:

.. code-block:: bash

   mpiexec --report-bindings --bind-to core --map-by slot:PE=$OMP_NUM_THREADS -n $NRANKS ./<application_executable> <input_file>


.. important::

   More sophisticated architectures (e.g., clusters with NUMA nodes) require careful specification of MPI and thread placement (e.g. binding of MPI ranks to sockets).


Simulation I/O
--------------

Simulation inputs
_________________

``chombo-discharge`` simulations take their input from a single simulation input file (possibly appended with overriding options on the command line).
Simulations may consist of several hundred possible switches for altering the behavior of a simulation, and all physics models in ``chombo-discharge`` are therefore equipped with Python setup tools that collect all such options into a single file when setting up a new application.
Generally, these input parameters are fetched from the options file of component that is used in a simulation. 
Simulation options usually consist of a prefix, a suffix, and a configuration value.
For example, the configuration options that adjusts the number of time steps that will be run in a simulation is

.. code-block:: none

   Driver.max_steps = 100

Likewise, for controlling how often plot are written:

.. code-block:: none

   Driver.plot_interval = 5

You may also pass input parameters through the command line. For example, running

.. code-block:: bash

   mpirun -np 32 <application_executable> <input_file> Driver.max_steps=10

will set the ``Driver.max_steps`` parameter to 10.
Command-line parameters override definitions in the input file.
Moreover, parameters parsed through the command line become static parameters, i.e. they are not run-time configurable (see :ref:`Chap:RuntimeConfig`).
Also note that if you define a parameter multiple times in the input file, the last definition is canon. 		

Simulation outputs
__________________

Mesh data from ``chombo-discharge`` simulations is by default written to HDF5 files, and if HDF5 is disabled ``chombo-discharge`` will not write any plot or checkpoint files. 
In addition to plot files, MPI ranks can output information to separate files so that the simulation progress can be tracked.

``chombo-discharge`` comes with controls for adjusting output.
Through the :ref:`Chap:Driver` class the user may adjust the option ``Driver.output_directory`` to specify where output files will be placed.
This directory is relative to the location where the application is run.
If this directory does not exist, ``chombo-discharge`` will create it. 
It will also create the following subdirectories given in :ref:`Tab:OutputDirectories`.

.. _Tab:OutputDirectories:
.. list-table:: Simulation output organization.
   :widths: 10 70
   :header-rows: 1

   * - Folder
     - Explanation
   * - :file:`chk`
     - Checkpoint files (these are used for restarting simulations from a specified time step). 
   * - :file:`crash`
     - Plot files written if a simulation crashes. 
   * - :file:`geo`
     - Plot files for geometries (if you run with ``Driver.geometry_only = true``). 
   * - :file:`mpi`
     - Information about individual MPI ranks, such as computational loads or memory consumption per rank. 
   * - :file:`plt`
     - All plot files.
   * - :file:`regrid`
     - Plot files written during regrids (if you run with ``Driver.write_regrid_files``).
   * - :file:`restart`
     - Plot files written during restarts (if you run with ``Driver.write_regrid_files``).

The reason for the output folder structure is that ``chombo-discharge`` can end up writing thousands of files per simulation and we feel that having a directory structure helps us navigate simulation data.  

Fundamentally, there are only two types of HDF5 files written:

1. Plot files, containing plots of simulation data.
2. Checkpoint files, which are binary files used for restarting a simulation from a given time step. 

The :ref:`Chap:Driver` class is responsible for writing output files at specified intervals, but the user is generally speaking responsible for specifying what goes into the plot files.
Since not all variables are always of interest, solver classes have options like ``plt_vars`` that specify which output variables in the solver will be written to the output file.
For example, one of our convection-diffusion-reaction solver classes have the following output options:

.. code-block:: text

   CdrGodunov.plt_vars = phi vel dco src ebflux # Plot variables. Options are 'phi', 'vel', 'dco', 'src', 'ebflux'

where ``phi`` is the state density, ``vel`` is the drift velocity, ``dco`` is the diffusion coefficient, ``src`` is the source term, and ``ebflux`` is the flux at embedded boundaries.
If you only want to plot the density, then you should put ``CdrGodunov.plt_vars = phi``.
An empty entry like ``CdrGodunov.plt_vars =`` may lead to run-time errors, so if you do not want a class to provide plot data you may put ``CdrGodunov.plt_vars = none``. 


.. _Chap:pout:

Parallel processor verbosity
____________________________

By default, ``Chombo`` will write a process output file *per MPI process* and this file will be named :file:`pout.n` where ``n`` is the MPI rank.
These files are written in the directory where you executed your application, and are *not* related to plot files or checkpoint files.
However, ``chombo-discharge`` prints information to these files as simulations advance (for example by displaying information of the current time step, or convergence rates for multigrid solvers).
To see information regarding the latest time steps, simply print a few lines in these files, e.g.

.. code-block:: bash

   tail -200 pout.0

While it is possible to monitor the evolution of ``chombo-discharge`` for each MPI rank, most of these files contain redundant information.
To adjust the number of files that will be written, ``Chombo`` can read an environment variable ``CH_OUTPUT_INTERVAL`` that determines which MPI ranks write :file:`pout.n` files. 
For example, if you only want the master MPI rank to write :file:`pout.0`, you would do

.. code-block:: bash

   export CH_OUTPUT_INTERVAL=999999999

.. important::
   
   If you run simulations at high concurrencies, you *should* turn off the number of process output files since they impact the performance of the file system. 
   
.. _Chap:RestartingSimulations:

Restarting simulations
______________________

Restarting simulations is done in exactly the same way as running simulations, although the user must set the ``Driver.restart`` parameter.
For example,

.. code-block:: bash

   mpirun -np 32 <application_executable> <input_file> Driver.restart=10

will restart from step 10.

Specifying anything but an integer is an error.
When a simulation is restarted, ``chombo-discharge`` will look for a checkpoint file with the ``Driver.output_names`` variable and the specified restart step.
It will look for this file in the subfolder :file:`/chk` relative to the execution directory.

If the restart file is not found, restarting will not work and ``chombo-discharge`` will abort.
You must therefore ensure that your executable can locate this file.
This also implies that you cannot change the ``Driver.output_names`` or ``Driver.output_directory`` variables during restarts, unless you also change the name of your checkpoint file and move it to a new directory.

.. note::

   If you set ``Driver.restart=0``, you will get a fresh simulation.

.. _Chap:RuntimeConfig:

Run-time configurations
_______________________

``chombo-discharge`` reads input parameters before the simulation starts, but also during run-time. 
This is useful when your simulation waited 5 days in the queue on a cluster before starting, but you forgot to tweak one parameter and don't want to wait another 5 days.

``Driver`` re-reads the simulation input parameters after every time step.
The new options are parsed by the core classes ``Driver``, ``TimeStepper``, ``AmrMesh``, and ``CellTagger`` through special routines ``parseRuntimeOptions()``.
Note that not all input configurations are suitable for run-time configuration.
For example, increasing the size of the simulation domain does not make sense but changing the blocking factor, refinement criteria, or plot intervals do.
To see which options are run-time configurable, see :ref:`Chap:Driver`, :ref:`Chap:AmrMesh`, or the :ref:`Chap:TimeStepper` and :ref:`Chap:CellTagger` that you use.

.. _Chap:Visualization:

Visualization
-------------

``chombo-discharge`` output files are always written to HDF5.
The plot files will reside in the ``plt`` subfolder where the application was run.

Currently, we have only used `VisIt <https://visit-dav.github.io/visit-website/>`_ for visualizing the plot files.
Learning how to use VisIt is not a part of this documentation; there are great tutorials on the `VisIt website <https://visit-dav.github.io/visit-website/>`_.
