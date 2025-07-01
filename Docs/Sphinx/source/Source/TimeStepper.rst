.. _Chap:TimeStepper:

TimeStepper
***********

``TimeStepper`` represents the EB+AMR equation solving class in ``chombo-discharge``.


It owns the various numerical solvers and is responsible for setting up solvers and advancing the equations of motion.
Because ``TimeStepper`` and not :ref:`Chap:Driver` owns the solvers as class members, it will also coordinate I/O (together with :ref:`Chap:Driver`). 

.. tip::

   `TimeStepper C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classTimeStepper.html>`_

Basic functions
===============

Since it is necessary to implement different solvers for different types of physics, ``TimeStepper`` is an abstract class with the following pure functions:

.. code-block:: c++

   // Setup routines
   virtual void setupSolvers() = 0;
   virtual void allocate() = 0;
   virtual void initialData() = 0;
   virtual void postInitialize() = 0;
   virtual void postCheckpointSetup() = 0;
   virtual void registerRealms() = 0;
   virtual void registerOperators() = 0;
   virtual void parseRuntimeOptions();

   // IO routines
   virtual void writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const = 0;
   virtual void readCheckpointData(HDF5Handle& a_handle, const int a_lvl) = 0;
   virtual void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const = 0;
   virtual int  getNumberOfPlotVariables() const = 0;
   virtual Vector<long int> getCheckpointLoads(const std::string a_realm, const int a_level) const;
   
   // Advance routines
   virtual Real computeDt() = 0;
   virtual Real advance(const Real a_dt) = 0;
   virtual void synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) = 0;
   virtual void printStepReport() = 0;
   
   // Regrid routines
   virtual void preRegrid(const int a_lmin, const int a_oldFinestLevel) = 0;
   virtual void postRegrid() = 0;
   virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) = 0;
   virtual bool needToRegrid();
   virtual bool loadBalanceThisRealm(const std::string a_realm) const;
   virtual void loadBalanceBoxes(Vector<Vector<int> >&            a_procs,
                                 Vector<Vector<Box> >&            a_boxes,
				 const std::string                a_realm,
				 const Vector<DisjointBoxLayout>& a_grids,
				 const int                        a_lmin,
				 const int                        a_finestLevel);

These functions are all used by :ref:`Chap:Driver` class at various stages.

Setup routines
==============

Here, we consider the various setup routines in ``TimeStepper``.
The routines are used by :ref:`Chap:Driver` in the simulation setup step.

.. _Chap:RegisterRealms:
   
registerRealms
--------------

``TimeStepper`` permits solvers to run on different realms (see :ref:`Chap:Realm`) for individual load balancing of various components.
To register a ``Realm``, users will have ``TimeStepper`` register realms in the ``registerRealms()``, as follows:

.. code-block:: c++

   void myTimeStepper::registerRealms(){
      m_amr->registerRealm(Realm::Primal);
      m_amr->registerRealm("particleRealm");
      m_amr->registerRealm("otherParticleRealm");
   }

Since at least one realm is required, :ref:`Chap:Driver` will *always* register the realm ``"Primal"``.
Fundamentally, there is no limitation to the number of realms that can be allocated.

.. _Chap:RegisterOperators:

registerOperators
-----------------

Internally, an instantiation of :ref:`Chap:Realm` contains the grids and the geometric information (e.g. ``EBISLayout``), as well as any operators that the user has seen fit to *register*.
Various operators are available for e.g. gradient stencils, conservative coarsening, ghost cell interpolation, filling a patch with interpolation data, redistribution, and so on.
Since operators always incur overhead and not all applications require *all* operators, they must be *registered*. 
If a solver needs an operator for, say, piecewise linear ghost cell interpolation, the solver needs to *register* that operator through the ``AmrMesh`` public interface:

.. code-block:: c++

   m_amr->registerOperator(s_eb_pwl_interp, m_realm, m_phase);

Once an operator has been registered, ``Realm`` will define those operators during initialization e.g. regrids.
Run-time error messages are issued if an AMR operator is used, but has not been registered.

More commonly, ``chombo-discharge`` solvers will contain a routine that registers the operators that the solver needs.
A valid ``TimeStepper`` implementation *must* register all required operators in the function ``registerOperators()``. 

Currently available operators are:

#. Gradient ``s_eb_gradient``.
#. Irregular cell centroid interpolation, ``s_eb_irreg_interp``.
#. Coarse grid conservative coarsening, ``s_eb_coar_ave``.
#. Piecewise linear interpolation (with slope limiters), ``s_eb_fill_patch``.
#. Linear ghost cell interpolation, ``s_eb_fine_interp``.
#. Flux registers, ``s_eb_flux_reg``.
#. Redistribution registers, ``s_eb_redist``.
#. Non-conservative divergence stencils, ``s_eb_noncons_div``.
#. Multigrid interpolators, ``s_eb_multigrid`` (used for multigrid).     
#. Signed distance function defined on grid, ``s_levelset``.
#. Particle-mesh support, ``s_eb_particle_mesh``.   

Solvers will typically allocate a subset of these operators, but for multiphysics code that use both fluid and particles, most of these will probably be in use.

setupSolvers
------------

``setupSolvers`` is used for setting up solvers.
This step is done *prior* to setting up the grids, so it is not possible to allocate mesh data inside this routine.

allocate
--------

``allocate`` is used for allocating particle and mesh data for the solvers and ``TimeStepper``.
This step is done *after* the grids have been initialized by :ref:`Chap:AmrMesh` *and* during regrids. 

initialData
-----------

``initialData`` is called by :ref:`Chap:Driver` setup routines after the ``allocate`` step.
This routine must fill solvers with initial data.

postInitialize
--------------

``postInitialize`` is called *after* :ref:`Chap:Driver` has filled the solvers with initial data.
Most data initialization steps can, however, be done in ``initialData``.

postCheckpointSetup
-------------------

During simulation restarts, :ref:`Chap:Driver` will open an HDF5 file and have ``TimeStepper`` fill solvers with data from that file.
``postCheckpointSetup`` is a routine which is called immediately after the solvers have been filled with data. 

I/O routines
============

The ``TimeStepper`` I/O routines serves two purposes:

#. To add solver data to HDF5 plot files.
#. To write and read data for checkpoint/restart files.

In general, plot and checkpoint data do not contain the same data.

getNumberOfPlotVariables
------------------------

``getNumberOfPlotVariables`` must return the number of components that will be plotted by ``TimeStepper``.
Note that if ``TimeStepper`` will plot a single scalar, it must return a value of one.
If it plots a single vector, it must return a value of ``SpaceDim``.

The existence of ``getNumberOfPlotVariables`` is due to pre-allocation of memory that will be written to the plot file. 

writePlotData
-------------

``writePlotData`` will write the plot data to the provided data holder.
The signature is

.. code-block:: c++
   
   virtual void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const = 0;

Here, ``a_output`` is pre-allocate block of memory that ``TimeStepper`` will write its components to, and ``a_plotVariable`` are the associated plot variable names.
``a_icomp`` is the starting component in ``a_output`` where we start writing data.

writeCheckpointData
-------------------

``writeCheckpointData`` will write solver data to the provided HDF5 file handle.
This data is used when restarting simulations from a checkpoint file. 
Note that checkpoint data is written on a level-by-level basis.

readCheckpointData
------------------

``readCheckpointData`` will read data from the provided HDF5 file handle and back into the solvers.
Note that the data is read on a level-by-level basis.

Advance routines
================

computeDt
---------

``computeDt`` will compute a time step for :ref:`Chap:Driver` to use when calling the ``advance`` method. 

advance
-------

``advance`` is called by :ref:`Chap:Driver` when advancing the equation of motion one time step.
Note that ``advancpe`` takes a trial time step as input argument and returns the actual time step that was used.
These do not need to be the same. 

synchronizeSolverTimes
----------------------

``synchronizeSolverTimes`` is called after the ``advance`` method and is used to update the simulation time for all solvers. 

printStepReport
---------------

``printStepReport`` called after the ``advance`` method -- it can be left empty but is otherwise used to print some information about the time step that was taken. 

Regrid routines
===============

For an explanation to how regridding occurs in ``chombo-discharge``, see :ref:`Chap:Regridding`.

preRegrid
---------

``preRegrid`` should any necessary pre-regrid operation.
Note that when solvers regrid their data, solution is allocated on new grids and the previously defined data is lost.
For this reason most solvers have the option of putting the old grid data into temporary storage that permits us to interpolate to the new grids. 

regrid
------

``regrid`` will perform the actual regrid operation.

postRegrid
----------

``postRegrid`` is called after ``regrid`` can be used to perform any post-regrid operations.


Load balancing routines
=======================

During the regrid step, :ref:`Chap:Driver` will check if any of the realms should be load balanced.
If a realm should be load balanced then ``TimeStepper`` take a ``DisjointBoxLayout`` which originally load balanced using the patch volume, and generate a new set of grids.

loadBalanceThisRealm
--------------------

Return true if a :ref:`Chap:Realm` should be load balanced and false otherwise.

loadBalanceBoxes
----------------

This is called if ``loadBalanceThisRealm`` evaluates to true, and in this case the ``TimeStepper`` should compute a new set of rank ownership for the input grid boxes. 
