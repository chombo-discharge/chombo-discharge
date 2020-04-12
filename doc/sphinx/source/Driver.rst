.. _Chap:driver:

driver
======

The ``driver`` class is the class that runs `PlasmaC` simulations and is defined in :file:`/src/driver/driver.cpp(H)`.
The constructor for this class is

.. code-block:: c++

   driver(const RefCountedPtr<computational_geometry>& a_compgeom,
          const RefCountedPtr<time_stepper>&           a_timestepper,
	  const RefCountedPtr<amr_mesh>&               a_amr,
	  const RefCountedPtr<cell_tagger>&            a_celltagger = RefCountedPtr<cell_tagger>(NULL),
	  const RefCountedPtr<geo_coarsener>&          a_geocoarsen = RefCountedPtr<geo_coarsener>(NULL));

Observe that the ``driver`` class does not *require* an instance of :ref:`Chap:cell_tagger`.
If users decide to omit a cell tagger, regridding functionality is completely turned off and only the initially generated grids will be used.

The usage of the ``driver`` class is primarily object construction with dependency injection of the geometry, the physics (i.e. ``time_stepper``), the ``amr_mesh`` instance, and possibly a cell tagger.
Usually, only a single routine is used:

.. code-block:: c++

   void setup_and_run();

This routine will set up and run a simulation.
Simulation setup depends on the way a simulation is run.

How fresh simulations are set up
--------------------------------

If a simulation starts from the first time step, the ``driver`` class will perform the following major steps within ``setup_and_run()``.

1. Ask ``computational_geometry`` to generate the cut-cell moments.
2. Collect all the cut-cells and ask ``amr_mesh`` to set up an initial grid where all the cut-cells are refined.
   It is possible to remove some of the cut-cell refinement flags through the auxiliary class ``geo_coarsener``.
3. Ask the ``time_stepper`` to set up all the relevant solvers and fill them with initial data.
4. Perform the number of regrids that the user asks for.

Step 3 will differ significantly depending on the physics that is solved for.

How simulations are restarted
-----------------------------

If a simulation *does not start* the first time step, the ``driver`` class will perform the following major steps within ``setup_and_run()``.

1. Ask ``computational_geometry`` to generate the cut-cell moments.
2. Read a checkpoint file that contains the grids and all the data that have been checkpointed by the solvers. 
3. Ask the ``time_stepper`` to perform a "post-checkpoint" step to initialize any remaining data so that a time step can be taken.
   This functionality has been included because not all data in every solver needs to be checkpointed.
   For example, an electric field solver only needs to write the electric potential to the checkpoint file because the electric field is simply obtained by taking the gradient.
4. Perform the number of initial regrids that the user asks for.

How simulations are run
-----------------------

The algorithm for running a simulation is very simple; the ``driver`` class simply calls ``time_stepper`` for computing a reasonable time step for advancing the equations, and then it asks ``time_stepper`` to actually perform the advance. 
Regrids, plot files, and checkpoint files are written at certain step intervals.
In essence, the algorithm looks like this:

.. code-block:: c++

   driver::run(...){

      while(KeepRunningTheSimulation){
         if(RegridEverything){
	    driver->regrid()
	 }

	 time_stepper->computeTimeStep()
	 time_stepper->advanceAllEquationsOneStep()

         if(WriteAPlotFile || EndOfSimulatoin){
	    driver->writePlotFile();
	 }
	 if(TimeToWriteACheckpointFile || EndOfSimulation){
	    driver->writeCheckpointFile()
	 }

	 KeepRunningTheSimulation = true or false
      }
   }

How regrids are performed
-------------------------

Regrids are called by the ``driver`` class and occur as follows in ``driver::regrid(...)``:

1. Ask ``cell_tagger`` to generate tags for grid refinement and coarsening.
2. The ``time_stepper`` class stores data that is subject to regrids.
   This data various between the type of solvers.
   For grid-based solvers, e.g. CDR solvers, the scalar :math:`\phi` is copied into a scratch space.
   The reason for this backup is that :math:`\phi` will be allocated on the *new* AMR grids, but we must still have access to the previously defined data in order to interpolate to the new grids.
3. If necessary, ``time_stepper`` can deallocate unecessary storage.
   Implementing a deallocation function for ``time_stepper``-derived classes is not a requirement, but can in certain cases be useful, for example when using the Berger-Rigoutsous algorithm at large scale.
4. The ``amr_mesh`` class generates the new grids and defines new AMR operators.
5. The ``time_stepper`` class regrids its solvers and internal data.

In C++ pseudo-code, this looks something like:

.. code-block:: c++

   driver::regrid(){

      // Tag cells
      cell_tagger->tagCellsForRefinement() 

      // Store old data and free up some memory
      time_stepper->storeOldGridData()
      time_stepper->deallocateUnneccesaryData()

      // Generate the new grids
      amr_mesh->regrid() 

      // Regrid physics and all solvers
      time_stepper->regrid()
   }

The full code is defined in ``driver::regrid()`` in file :file:`/src/driver/driver.cpp`. 


Class options
-------------

Various class options are available for adjusting the behavior of the ``driver`` class.

.. literalinclude:: links/driver.options

We now discuss some of these options

* ``verbosity`` describes 
