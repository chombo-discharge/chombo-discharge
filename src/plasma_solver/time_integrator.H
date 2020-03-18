/*!
  @file   time_integrator.H
  @brief  Declaration of an abstract class for integrating the equations of motion
  @author Robert Marskar
  @data   March 2020
*/

#ifndef _TIME_INTEGRATOR_
#define _TIME_INTEGRATOR_

#include "computational_geometry.H"
#include "plasma_kinetics.H"
#include "amr_mesh.H"
#include "mfis.H"
#include "cdr_layout.H"
#include "rte_layout.H"
#include "poisson_solver.H"
#include "sigma_solver.H"

#include <RealVect.H>

/*!
  @brief Abstract class for evolving plasma equations in time
*/
class time_integrator {
public:

  time_integrator(){}

  virtual ~time_integrator(){};

  // Setup routines
  virtual void setup_solvers() = 0;
  virtual void initial_data() = 0;

  // IO routines
  virtual void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const = 0;
  virtual void read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) = 0;
  virtual int get_num_plot_vars() const = 0;
  virtual void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const = 0;

  // Advance routines
  virtual void compute_dt(Real& a_dt, time_code::which_code& a_timecode) = 0;
  virtual Real advance(const Real a_dt) = 0;


  // Regrid routines
  virtual void cache_internals() = 0;
  virtual void deallocate_internals() = 0;
  virtual void cache_states() = 0;
  virtual void deallocate_solver_internals() = 0;
  virtual void regrid_solvers(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) = 0;
  virtual void allocate_internals() = 0;
  virtual void regrid_internals(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) = 0;
  virtual void post_regrid() = 0;

protected:
};

#endif
