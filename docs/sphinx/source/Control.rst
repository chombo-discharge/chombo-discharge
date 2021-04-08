.. _Chap:Control:

Controlling ``PlasmaC``
=======================

In this chapter we show how to run a ``PlasmaC`` simulation and control its behavior through input scripts or command line options.

Basic compiling and running
---------------------------

To run simulations, the user must first compile his application. Once the application has been defined, the user may compile is by

.. code-block:: bash

   make -s -j 32 DIM=N <application_name>

where *N* may be 2 or 3, and <application_name> is the name of the file that holds the ``main()`` function. This will compile an executable whose name depends on your application name and compiler settings. Please refer to the Chombo manual for explanation of the executable name. You may, of course, rename your application.

Next, applications are run by

.. code-block:: bash

   mpirun -np 32 <application_executable> <input_file>

where <input_file> is your input file. On clusters, this is a little bit different and usually requires passing the above command through a batch system. Note that if you define a parameter multiple times in the input file, the last definition is canon. 

You may also pass input parameters through the command line. For example, running

.. code-block:: bash

   mpirun -np 32 <application_executable> <input_file> driver.max_steps=10

will set the ``driver.max_steps`` parameter to 10. Command-line parameters override definitions in the input_file. 

.. _Chap:ControllingOutput:

Controlling output
------------------

`PlasmaC` comes with controls for adjusting output. Through the :ref:`Chap:driver` class the user may adjust the option ``driver.output_directory`` to specify where output files will be placed. This directory is relative to the location where the application is run. If this directory does not exist, PlasmaC does it's best at creating it. In addition, it will create four more directories

* :file:`output_directory/plt` contains all plot files.
* :file:`output_directory/chk` contains all checkpoint files, which are used for restarting.
* :file:`output_directory/mpi` contains information about individual MPI ranks. 
* :file:`output_directory/geo` contains geometric files that are written by PlasmaC (if you enable ``driver.write_ebis``). 

The files in :file:`output_directory/geo` do *not* represent your geometry in the form of level sets. Instead, the files that are placed here are HDF5 representations of your embedded boundary graph, which can be *read* by PlasmaC if you enable ``driver.read_ebis``. This is a shortcut that allows faster geometry generation when you restart simulations, but geometry generation is typically so fast that it is never used. 

The reason for this structure is that PlasmaC can end up writing thousands of files per simulation and we feel that having a directory structure helps us navigate simulation data.

The driver class :ref:`Chap:driver` is responsible for writing output files at specified intervals, but the user is responsible for specifying what goes into those files. Since not all variables are always of interest, the solver classes themselves have options ``plt_vars`` that specify which output variables will be written to the output file. For example, our convection-diffusion-reaction solver classes have the following output options:

.. code-block:: bash

   cdr_gdnv.plt_vars = phi vel dco src ebflux # Plot variables. Options are 'phi', 'vel', 'dco', 'src', 'ebflux'

where ``phi`` is the state density, ``vel`` is the drift velocity, ``dco`` is the diffusion coefficient, ``src`` is the source term, and ``ebflux`` is the flux at embedded boundaries. Which variables are available for output changes for one class to the next. If you only want to plot the density, then you should put ``cdr_gdnv.plt_vars = phi``. An empty entry like ``cdr_gdnv.plt_vars =`` will lead to run-time errors, so if you do not want a class to provide plot data you may put ``cdr_gdnv.plt_vars = -1``. 


Controlling processor output
----------------------------

By default, Chombo will write a process output file *per MPI process* and this file will be named :file:`pout.n` where ``n`` is the MPI rank. These files are written in the directory where you executed your application, and are *not* related to plot files or checkpoint files. However, PlasmaC prints information to these files as simulations advance (for example by displaying information of the current time step, or convergence rates for multigrid solvers). While it is possible to monitor the evolution of PlasmaC through each MPI, most of these files contain redundant information. To turn off the number of files that will be written, Chombo can read an environment variable ``CH_OUTPUT_INTERVAL``. For example, if you only want the master MPI rank to write :file:`pout.0`, you would do

.. code-block:: bash

   export CH_OUTPUT_INTERVAL=999999999

You can, of course, put the definition in your :file:`.bashrc` file (for Bourne shell). Note that if you run simulations at high concurrencies, you *should* turn off the number of process output files since they impact the performance of the file system. 
   
.. _Chap:RestartingSimulations:

Restarting simulations
----------------------

Restarting simulations is done in exactly the same way as running simulations, although the user must set the ``driver.restart`` parameter. For example,

.. code-block:: bash

   mpirun -np 32 <application_executable> <input_file> driver.restart=10

will restart from step 10. If you set ``driver.restart=0``, you will get a fresh simulation. Specifying anything but an integer is an error. When a simulation is restarted, PlasmaC will look for a checkpoint file with the ``driver.output_names`` variable and the specified restart step. If this file is not found, restarting will not work and ``PlasmaC`` will abort. You must therefore ensure that your executable can locate this file. This also implies that you cannot change the ``driver.output_names`` or ``driver.output_directory`` variables during restarts, unless you also change the name of your checkpoint file and move it to a new directory.

Visualization
-------------

`PlasmaC` output files are written to HDF5 files in the format ``<simulation_name>.step#.dimension.hdf5`` and the files will be written to the directory specified by :ref:`Chap:driver` runtime parameters. Currently, we have only used VisIt for visualizing the plot files.    
