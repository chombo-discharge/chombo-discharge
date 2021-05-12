/*!
  @file   time_stepper.H
  @brief  Declaration of an abstract class for integrating the equations of motion
  @author Robert Marskar
  @date   March 2020
*/

#ifndef _TIME_STEPPER_
#define _TIME_STEPPER_

#include "computational_geometry.H"
#include "CD_AmrMesh.H"
#include "mfis.H"

#include "CD_NamespaceHeader.H"
  
/*!
  @brief Time code for output
*/
enum class time_code {
  advection,
  diffusion,
  advection_diffusion,  
  source,
  relaxation_time,
  restricted,
  hardcap,
  error,

  physics
};

/*!
  @brief Abstract class for performing time stepping
*/
class time_stepper {
public:

  time_stepper(){}

  virtual ~time_stepper(){};

  // Setup routines
  virtual void setup_solvers() = 0;
  virtual void allocate() = 0;
  virtual void initial_data() = 0;
  virtual void post_initialize() = 0;
  virtual void post_checkpoint_setup() = 0;
  virtual int get_redistribution_regsize() const;

  // Registration routines
  virtual void registerRealms() = 0;
  virtual void registerOperators() = 0;

  // Runtime options
  virtual void parse_runtime_options();

  // IO routines
  virtual void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const = 0;
  virtual void read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) = 0;
  virtual void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const = 0;
  virtual int  get_num_plot_vars() const = 0;
  virtual Vector<long int> get_checkpoint_loads(const std::string a_realm, const int a_level) const;

  // Advance routines
  virtual void compute_dt(Real& a_dt, time_code& a_timecode) = 0;
  virtual Real advance(const Real a_dt) = 0;
  virtual void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) = 0;
  virtual void print_step_report() = 0;

  // Regrid routines
  virtual void deallocate() = 0;
  virtual void pre_regrid(const int a_lmin, const int a_old_finest_level) = 0;
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) = 0;

  virtual bool need_to_regrid();
  virtual bool load_balance_realm(const std::string a_realm) const;
  virtual void load_balance_boxes(Vector<Vector<int> >&            a_procs,
				  Vector<Vector<Box> >&            a_boxes,
				  const std::string                a_realm,
				  const Vector<DisjointBoxLayout>& a_grids,
				  const int                        a_lmin,
				  const int                        a_finest_level);
  virtual void post_regrid() = 0;

  // Set routines
  virtual void set_amr(const RefCountedPtr<AmrMesh>& a_amr);
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
  time_code m_timecode;

  /*!
    @brief AMR
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Physical domain
  */
  RefCountedPtr<computational_geometry> m_compgeom;
};
#include "CD_NamespaceFooter.H"

#endif
