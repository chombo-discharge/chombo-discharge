.. _Chap:time_stepper:

time_stepper
============

The ``time_stepper`` class is the physics class in ``PlasmaC`` - it has direct responsibility of setting up the solvers and performing time steps.

Since it is necessary to implement different solvers for different types of physics, ``time_stepper`` is an abstract class with the following pure functions:

.. code-block:: c++

  // Setup routines
  virtual void setup_solvers() = 0;
  virtual void initial_data() = 0;
  virtual void post_checkpoint_setup() = 0;

  // IO routines
  virtual void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const = 0;
  virtual void read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) = 0;
  virtual int get_num_plot_vars() const = 0;
  virtual void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const = 0;

  // Advance routines
  virtual void compute_dt(Real& a_dt, time_code::which_code& a_timecode) = 0;
  virtual Real advance(const Real a_dt) = 0;
  virtual void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) = 0;
  virtual void print_step_report() = 0;

  // New regrid routines
  virtual void register_operators() = 0; 
  virtual void pre_regrid(const int a_lmin, const int a_old_finest_level) = 0;
  virtual void deallocate() = 0;
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) = 0;


These functions are used in the ``driver`` class at various stages.
The three functions in the category *setup routines* are, for example using during simulation setup or after reading a checkpoint file for simulation restarts.
The IO routines are there so that users can choose which solvers perform any output, and the advance routines are there such that the user can implement new algorithms for time integration.
Finally, the *regrid routines* are there so that the solvers can back up their data before the old grids are destroyed (``cache()``), deallocate unnecessary memory (``deallocate``), and regrid data onto the new AMR grids ``regrid(...)``.

Depending on the physics that is resolved, writing ``time_stepper`` can be a small or a big task.
For example, the code used for advection-diffusion in :file:`/physics/advection_diffusion/` is only a few hundred lines where most of the code is simply calling member functions from :ref:`Chap:CDR`. 
