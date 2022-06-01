.. _Chap:TimeStepper:

TimeStepper
============

The ``TimeStepper`` class is the physics class in ``chombo-discharge`` - it has direct responsibility of setting up the solvers and performing time steps.

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
   virtual void computeDt(Real& a_dt, TimeCode& a_timeCode) = 0;
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

These functions are all used in the ``Driver`` class at various stages.
The three functions in the category *setup routines* are, for example using during simulation setup or after reading a checkpoint file for simulation restarts.
The IO routines are there so that users can choose which solvers perform any output, and the advance routines are there such that the user can implement new algorithms for time integration.
Finally, the *regrid routines* are there so that the solvers can back up their data before the old grids are destroyed (``preRegrid()``) and the data is interpolated data onto the new AMR grids ``regrid(...)``.

To see how these functions can be implemented, see :ref:`Chap:Tutorial`. 
