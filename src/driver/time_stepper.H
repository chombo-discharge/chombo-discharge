/*!
  @file   time_stepper.H
  @brief  Declaration of an abstract class for integrating the equations of motion
  @author Robert Marskar
  @data   March 2020
*/

#ifndef _TIME_STEPPER_
#define _TIME_STEPPER_

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
class time_stepper {
public:

  time_stepper(){}

  virtual ~time_stepper(){};

  // Setup routines
  virtual void setup_solvers() = 0;
  virtual void initial_data() = 0;
  virtual void post_checkpoint_setup() = 0;
  virtual int get_redistribution_regsize() const;

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
  virtual bool need_to_regrid();
  virtual void cache() = 0;
  virtual void deallocate() = 0;
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) = 0;

  /*!
    @brief Set amr
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Set the physical domain
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom);

protected:

  /*!
    @brief Class verbosity
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_step;
  
  /*!
    @brief TIme
  */
  Real m_time;

  /*!
    @brief Previous time step size
  */
  Real m_dt;

  /*!
    @brief Time code for
  */
  time_code::which_code m_timecode;

  /*!
    @brief AMR
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Physical domain
  */
  RefCountedPtr<computational_geometry> m_compgeom;
};

#endif
