.. _Chap:TimeStepper:

TimeStepper
***********

``TimeStepper`` is essentially the problem solving class ``chombo-discharge``.
The class is abstract and is subclassed by the programmer in order to solve new problems (or change time discretizations).
Typically, ``TimeStepper`` owns all numerical solvers, has responsibility for setting up solvers, regridding its internal state, and advancing the equations of motion.
Because ``TimeStepper`` and not :ref:`Chap:Driver` owns the solvers as class members, it will also partially coordinate I/O by providing data that :ref:`Chap:Driver` will add to HDF5 files. 

.. tip::

   Here is the `TimeStepper C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classTimeStepper.html>`_

Basic functions
===============

There are numerous functions that must be implemented and coordinated in order to provide a full-fledged ``TimeStepper`` for various problems.
Below, we include the current header file for ``TimeStepper``.

.. literalinclude:: ../../../../Source/Driver/CD_TimeStepper.H
   :language: c++   

Setup routines
==============

Here, we consider the various setup routines in ``TimeStepper``.
The routines are used by :ref:`Chap:Driver` in the simulation setup step, both for fresh simulation setups as well as restarts.

.. _Chap:RegisterRealms:
   
registerRealms
--------------

``chombo-discharge`` permits things to happen on different sets of grids where the the grids themselves cover the same physical region, but where MPI ownership of grids might change between grid sets (see :ref:`Chap:Realm` for details).
To register a :ref:`Chap:Realm`, users will have ``TimeStepper`` register realms in the ``registerRealms()`` routine, as follows:

.. code-block:: c++

   void myTimeStepper::registerRealms(){
      m_amr->registerRealm(Realm::Primal);
      m_amr->registerRealm("particleRealm");
      m_amr->registerRealm("otherParticleRealm");
   }

The above code will ensure that ``chombo-discharge`` generates three :ref:`Chap:Realm`, which can be individually load balanced. 
Since at least one realm is required, :ref:`Chap:Driver` will *always* register the realm ``"Primal"``.
Fundamentally, there is no limitation to the number of realms that can be allocated.

setupSolvers
------------

``setupSolvers`` is used for instantiating the solvers.
This routine is called *prior* to creating grids, so it is not possible to allocate mesh data for the solvers inside this routine.
The rationale for this design is that :ref:`Chap:AmrMesh` must know relatively early which part of the AMR infrastructure that will be instantiated, so the solvers are created before allocating the grids.
It is still quite possible to parse lots of data into the solvers, e.g., setting input variables.

To provide an example, the code snippet below shows the implementation of this routine for the ``TimeStepper`` implementation for the :ref:`Chap:AdvectionDiffusionModel`:

.. _List:AdvectionDiffusionStepperSetupSolvers:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``setupSolvers`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 137-163

.. _Chap:RegisterOperators:

registerOperators
-----------------

Internally, an instantiation of :ref:`Chap:Realm` will provide access to the grids and cut-cell information on that :ref:`Chap:Realm`, as well as any operators that the user has seen fit to *register*.
Various operators are available for, e.g., computing gradients, conservative coarsening, ghost cell interpolation, filling a patch with interpolation data, redistribution, particle-mesh operations, and so on.
Since operators always incur overhead and not all applications require *all* operators, they must be *registered*. 
If a solver needs an operator for, say, piecewise linear ghost cell interpolation, the solver needs to *register* that operator through the ``AmrMesh`` public interface:

.. code-block:: c++

   m_amr->registerOperator(s_eb_pwl_interp, m_realm, m_phase);

Once an operator has been registered, ``Realm`` will automatically instantiate those operators during initialization or regrid operations.
Run-time error messages are issued if an AMR operator is used, but has not been registered.

More commonly, ``chombo-discharge`` solvers will contain a routine that registers the operators that the solver needs.
A valid ``TimeStepper`` implementation *must* register all required operators in the function ``registerOperators()``, but this is usually done by letting the solvers register what they need.
Failure to do so will issue a run-time error.
Solvers will typically allocate a subset of these operators, but for multiphysics code that use both fluid and particles, most of these will be in use.
An example of this is given in :numref:`List:AdvectionDiffusionStepperRegister`, which shows how this routine is implemented for the :ref:`Chap:AdvectionDiffusionModel`:

.. _List:AdvectionDiffusionStepperRegister:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``registerOperators`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 176-186

allocate
--------

``allocate`` is used for allocating particle and mesh data that is required during simulations. 
This step is done *after* the grids have been initialized by :ref:`Chap:AmrMesh` *and* during regrids.
Again using :ref:`Chap:AdvectionDiffusionModel` as an example,

.. _List:AdvectionDiffusionStepperAllocate:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``allocate`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 188-197

In the above snippet, ``TimeStepper`` only calls the solver allocation function (solvers generally know how to allocate their own internal data).
For more complex problems this routine will probably allocate additonal data that only lives within ``TimeStepper`` (and not the solvers). 

initialData
-----------

``initialData`` is called by :ref:`Chap:Driver` setup routines after the ``allocate`` step, and has responsibility of setting up the problem with initial data.
This can occasionally be simple, or for coupled problems it might be highly complex.
For discharge problems, this can involve filling the solvers with initial densities, and solving the Poisson equation for obtaining the electric field.
A simpler example is again given by :ref:`Chap:AdvectionDiffusionModel`:

.. _List:AdvectionDiffusionStepperInitialData:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``initialData`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 199-226

The above code defers the initialization of the density in the advection-diffusion-reaction solver to the actual solver implementation.
Here, one could equally well have fetched the density directly from the solver and set it to something else by simply iterating through the grid cells (or setting it through other means).
In addition, source terms are set to zero, as are boundary fluxes. 

postInitialize
--------------

``postInitialize`` is a special routine that is called *after* :ref:`Chap:Driver` has filled the solvers with initial data *and* :ref:`Chap:Driver` is done with all the initial regrids.
While most data initialization steps can, however, be done in ``initialData``, the function is put there as an open door to the programmer for performing certain post-initialization functions that do not not need be performed in ``initialData`` (which is called once per initial regrid).
For example, the :ref:`Chap:DischargeInceptionModel` uses this function to compute several relevant quantities after the electric has been obtained in ``initialData``.

postCheckpointSetup
-------------------

During simulation restarts, :ref:`Chap:Driver` will open an HDF5 file and have ``TimeStepper`` fill solvers and its own internal data with data from that file.
``postCheckpointSetup`` is a routine which is called immediately after the solvers have performed this step.
Several gas discharge models use this function to compute the electric field from the potential that was saved in the HDF5 file.


I/O routines
============

``TimeStepper`` contains I/O routines primarily serves two purposes:

#. To provide data for HDF5 plot files, used for post-processing analysis.
#. To read and write data from HDF5 checkpoint files, which are used to restart simulations from a specified time step.

In general, plot and checkpoint data do not contain the same data, and ``TimeStepper`` therefore requires that plot and checkpoint files are filled separately.
We discuss these below:

getNumberOfPlotVariables
------------------------

``getNumberOfPlotVariables`` must return the number of components that will be plotted by ``TimeStepper``.
The reason why this routine exists is the :ref:`Chap:Driver` will pre-allocate the necessary memory on each AMR level, and ``TimeStepper`` will then copy solver data into this data holder. 
Specifically, if ``TimeStepper`` will plot a single scalar, it must return a value of one.
If it plots a single vector, it must return a value of ``SpaceDim``.
Below we include the implementation of this routine for the :ref:`Chap:AdvectionDiffusionModel`:

.. _List:AdvectionDiffusionStepperGetNumPlotVars:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``getNumberOfPlotVariables`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 280-290

getPlotVariableNames
--------------------

``getPlotVariableNames`` provides a list of plot variable names for the HDF5 file.
This list must have the same length as the returned value of ``getNumberOfPlotVariables``.
Below we include the implementation of this routine for the :ref:`Chap:AdvectionDiffusionModel`:

.. _List:AdvectionDiffusionStepperGetPlotVarNames:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``getPlotVariableNames`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 292-301

.. tip::

   When writing some vector data :math:`F` to the HDF5-file, one can write the variables as ``x-F``, ``y-F``, and ``z-F``, and VisIt visualization will automatically recognize :math:`F` as a vector field.


writePlotData
-------------

``writePlotData`` will write the plot data to the provided data holder.
The function signature is

.. literalinclude:: ../../../../Source/Driver/CD_TimeStepper.H
   :language: c++
   :lines: 160-171
   :dedent: 2
	   
In this function, ``a_output`` is pre-allocated block of memory that ``TimeStepper`` will write its components to (beginning at ``a_icomp``).
Note that if ``TimeStepper`` writes :math:`N` components, the implementation must increment ``a_icomp`` by :math:`N`.
Usually, solvers will have their own ``writePlotData`` routines which lets ``TimeStepper`` simply call the solver functions.
An example is given below for the :ref:`Chap:AdvectionDiffusionModel`:

.. _List:AdvectionDiffusionStepperWritePlotData:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``writePlotData`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 303-318

writeCheckpointData
-------------------

``writeCheckpointData`` must write necessary data for checkpointing the simulation state.
This data is used when restarting simulations from a checkpoint file. 
Note that checkpoint data is written on a level-by-level basis.
The function signature is

.. literalinclude:: ../../../../Source/Driver/CD_TimeStepper.H
   :language: c++
   :lines: 125-133
   :dedent: 2

Usually, the solvers know themselves what data to put in the checkpoint files and these routines are then pretty simple.
Below, we again include an example for the :ref:`Chap:AdvectionDiffusionModel`:

.. _List:AdvectionDiffusionStepperWriteChkData:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``writeCheckpointData`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 229-238

When implementing the function, it is important to add any data that is required when restarting the simulation.
It is also beneficial to *not* add data that is not required since this will just lead to larger files.
An example of this is the :ref:`Chap:FieldSolver` class, which only checkpoints the potential and not the electric field, as the latter is simply obtained by taking the gradient.

.. tip::
   
   Computational particles can also be added to the checkpoint files.


readCheckpointData
------------------

``readCheckpointData`` is the function that will read data from an HDF5 checkpoint file and populate ``TimeStepper`` with this data.
The data is read on a level-by-level basis, with a function signature

.. literalinclude:: ../../../../Source/Driver/CD_TimeStepper.H
   :language: c++
   :lines: 137-144
   :dedent: 2

Solvers will normally already know what data to read into their data members.
E.g., the example for the :ref:`Chap:AdvectionDiffusionModel` is

.. _List:AdvectionDiffusionStepperReadChkData:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``readCheckpointData`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 242-251

Advance routines
================

computeDt
---------

``computeDt`` is a routine that will compute a trial time step when calling the ``advance`` method.
We have chosen to call this a trial time step because

#. :ref:`Chap:Driver` might choose to use a smaller time step in order to write plot files at specific times.
#. When calling the actual advance method (see below), it is possible to return a different time step than the one computed through ``computeDt``.

The calculation of a time step can be quite involved, depending on the application being imlemented.
Moreover, many ``TimeStepper`` implementations will provide hooks for swapping algorithms, and in this case the time step might be limited differently.
For the :ref:`Chap:AdvectionDiffusionModel` the implementation is as follows:

.. _List:AdvectionDiffusionStepperComputeDt:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``computeDt`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 320-365

In the code above, ``TimeStepper`` supports both fully explicit advection-diffusion advances as well as split-step advances with implicit diffusion.
Depending on how the user chooses to run the code, the time step is therefore computed differently.
At the bottom of the above code, hard limits on the time step are also enforced.

advance
-------

The ``advance`` method has responsibility for advancing physics module one time step, and is called by ``Driver``.
The function signature is

.. literalinclude:: ../../../../Source/Driver/CD_TimeStepper.H
   :language: c++
   :lines: 211-218
   :dedent: 2

As mentioned in the documentation for this method, the function takes a trial time step ``a_dt`` which is the physical time step.
This time step is the one computed by ``computeDt``.
It is, however, quite possible to advance the equations of motion over a time that does not equal ``a_dt``, which is generally the case for adaptive time stepping methods.

The implementation of the ``advance`` method is usually the most time-consuming part of implementation a new ``TimeStepper``, and the implementation of this routine can become substantially complicated.
For simpler problems this routine is relatively straightforward to implement, however, also in a way that involves AMR and cut-cells.

.. _List:AdvectionDiffusionStepperAdvance:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``advance`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 367-376, 380-407, 456-459, 468-469

In the above code we have included the part of the ``advance`` routine that executes Heun's method for advancing the scalar.
This code also includes allocation of temporaries for computing the coefficients, storage for the intermediate state, and enforcement of boundary conditions, all of which include AMR.
At the way out of the routine the solution is coarsened and the ghost cells are updated, and the trial time step (``a_dt``) is returned.

synchronizeSolverTimes
----------------------

``synchronizeSolverTimes`` is called after the ``advance`` method and is used to update the simulation time for all solvers.
Again, this routine exists because there is often a physical time to be tracked by the solvers (e.g., enforcement of time-dependent voltage applications).
This routine simply ensures that ``TimeStepper`` and all solvers that ``TimeStepper`` owns see the same physical time, number of steps, and time step sizes.
The implementation for :ref:`Chap:AdvectionDiffusionModel` is

.. _List:AdvectionDiffusionStepperSync:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``synchronizeSolverTimes`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 471-484

printStepReport
---------------

``printStepReport`` is called after the ``advance`` method, and provides extra information printed to the ``pout.*`` files (see :ref:`Chap:Control`).
This function is called by :ref:`Chap:Driver` after performing a time step, and can be used to print extra information not covered by :ref:`Chap:Driver`, such as how the time step was limited, or other information that is useful for monitoring the behavior of ``TimeStepper``.
For example, the current gas discharge models in ``chombo-discharge`` print the maximum electric field and density at each time step.
Note that ``printStepReport`` has (or should have!) no side-effects that affect the simulation state.

Regrid routines
===============

The regrid routines in ``TimeStepper`` must, in combination, be able to transfer the simulation between old and new grids.
For an explanation to how regridding occurs in ``chombo-discharge``, see :ref:`Chap:DriverRegridding`.
In particular, when regrids occur the old grids are eventually destroyed so it is necessary to cache the old-grid simulation states so that we have something to interpolate from whan transfer the state to the new grids. 

preRegrid
---------

``preRegrid`` should any necessary pre-regrid operations that are necessary in order to call the ``regrid``.
This will virtually always include caching the old-grid simulation state, both for the solvers and also for internal data in ``TimeStepper`` 
Solvers usually know how to do this, and in some cases this function can be deceptively simple, as illustrated by the implementation of this function in :ref:`Chap:AdvectionDiffusionModel`:

.. _List:AdvectionDiffusionStepperPreRegrid:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``preRegrid`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 486-495

Other implementations of ``TimeStepper`` may have substantially more complicated ``preRegrid`` routines.

regrid
------

``regrid`` is the function that performs an actual regrid operation.
At the time when ``regrid`` is called, the old grids are already destroyed and are only available through the cached data.
Solvers are usually implemented with their own regrid routine, and if the only things that need to be regridded are the solvers, the implementation of this routine can be comparatively simple, as illustrated below for the :ref:`Chap:AdvectionDiffusionModel`:

.. _List:AdvectionDiffusionStepperRegrid:
.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp
   :caption: Implementation of the ``regrid`` routine for a simple advection-diffusion problem.
   :language: c++
   :lines: 497-516


For other ``TimeStepper`` implementations this routine can become much more complex.
The :ref:`Chap:ItoKMC`, for example, will regrid not only solver data but also internal mesh and particle data, recompute conductivities, deposit particles, handle superparticles, and prepare the simulation state for the next time step.

postRegrid
----------

The ``postRegrid`` is called after ``regrid`` has completed and can be used to perform any post-regrid specific procedures.
This function is not a pure function, and an implementation of this function is therefore not a requirement.


Load balancing routines
=======================

The default load-balancing method in ``chombo-discharge`` is to distribute grid patches equally among ranks, respecting a space-filling Morton curve on each grid level.
However, during the regrid step, :ref:`Chap:Driver` will check if meshes should be load balanced using different heuristics.
This load balancing can be done separately for each :ref:`Chap:Realm`, and in this case the MPI ranks will have different patch ownership in different grid sets.

If a realm should be load balanced with a different method than the default load balancing scheme, then ``TimeStepper`` can take a ``DisjointBoxLayout`` which originally load balanced using the patch volume, and regenerate the patch-to-rank ownership for the grids.
This functionality is implemented through two routines:

#. ``loadBalanceThisRealm`` which checks if a specific :ref:`Chap:Realm` should be load balanced.
#. ``loadBalanceBoxes`` which load balances the boxes on the specified :ref:`Chap:Realm`.

Note that these functions are not pure functions, and it is perfect fine to use their default implementation, in which case each MPI rank gets approximately the same number of grid patches. 

loadBalanceThisRealm
--------------------

The function signature for this function is

.. literalinclude:: ../../../../Source/Driver/CD_TimeStepper.H
   :language: c++
   :lines: 270-275
   :dedent: 2
	      
This function must return true if the input :ref:`Chap:Realm` (``a_realm``) should be load balanced.

loadBalanceBoxes
----------------

If ``loadBalanceThisRealm`` returns true, the following function is responsible for actually regenerating the grids:

.. literalinclude:: ../../../../Source/Driver/CD_TimeStepper.H
   :language: c++
   :lines: 277-296
   :dedent: 2

This is called if ``loadBalanceThisRealm`` evaluates to true, and in this case the ``TimeStepper`` should compute a new set of rank ownership for the input grid boxes.
Observe that ``loadBalanceBoxes`` occurs for the entire AMR hierarchy, where the outer vector of ``a_procs`` and ``a_boxes`` is the grid level, and the inner vectors describe the ownership of each box.
The default implementation of this function ensures that when we load balance a level, we account for the accumulated load on coarser levels (see :ref:`TimeStepperloadBalanceBoxes`).

.. _TimeStepperloadBalanceBoxes:
.. literalinclude:: ../../../../Source/Driver/CD_TimeStepper.H
   :caption: Default implementation of ``loadBalanceBoxes``.
   :language: c++
   :lines: 270-275
   :dedent: 2

In the above, we use the ``Loads`` class to hold the computational load for each rank, and on each level we compute the load for each patch to be equal to the number of grid cells in the patch.
This is later load balanced in ``LoadBalancing::makeBalance``, which is a routine that ensures that when we assign boxes on some grid level :math:`l`, we account for loads already assigned on coarser grid levels.
