/*!
  @file time_stepper.H
  @brief Declaration of a time_stepper, an abstract class for evolving in time across an amr mesh
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _TIME_STEPPER_
#define _TIME_STEPPER_

#include "computational_geometry.H"
#include "plasma_kinetics.H"
#include "physical_domain.H"
#include "amr_mesh.H"
#include "mfis.H"
#include "cdr_layout.H"
#include "rte_layout.H"
#include "poisson_solver.H"
#include "sigma_solver.H"

#include <RealVect.H>

/*!
  @brief Abstract class for evolving in time across an amr mesh
*/
class time_stepper {
public:

  /*!
    @brief Constructor
  */
  time_stepper();

  /*!
    @brief Destructor
  */
  ~time_stepper();

  /*!
    @brief Advance method, advances equations.
    @param[in] a_dt Time step
    @returns Actual time step used (can be different from a_dt for adaptive methods)
  */
  virtual Real advance(const Real a_dt){
#if 0
    MayDay::Warning("time_stepper::advance - this should be PURE");
#endif
  }

  /*!
    @brief Check if RTE solvers are stationary
  */
  virtual bool stationary_rte();

  /*!
    @brief Compute the time step and how it was restricted
  */
  virtual void compute_dt(Real& a_dt, time_code::which_code& a_timecode);

  /*!
    @brief Compute the cell-centered electric field on a specific phase
  */
  virtual void compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase, const MFAMRCellData& a_potential);

  /*!
    @brief Compute the cell-centered electric field on both phases
  */
  virtual void compute_E(MFAMRCellData& a_E, const MFAMRCellData& a_potential);

  /*!
    @brief Compute the current density
  */
  virtual void compute_J(EBAMRCellData& a_J);

  /*!
    @brief Compute photon source terms. 
  */
  virtual void compute_photon_source_terms(Vector<EBAMRCellData*>        a_source,
					   const Vector<EBAMRCellData*>& a_cdr_states,
					   const EBAMRCellData&          a_E,
					   const centering::which_center a_centering);

  /*!
    @brief Compute the cell-centered space charge density by using data inside cdr_layout. 
    @details This assumes that the data contain in cdr_layout is cell-centered
  */
  virtual void compute_rho();

  /*!
    @brief Compute the centroid-centered space charge density
    @details You have the option of specifying the centering for this computation. If the densities are cell-centered, we first 
    compute the cell-centered space charge density and the interpolate to the cell centroid. 
  */
  virtual void compute_rho(MFAMRCellData&                a_rho,
			   const Vector<EBAMRCellData*>  a_densities,
			   const centering::which_center a_centering);

  /*!
    @brief Instantiate solvers
  */
  virtual void instantiate_solvers();

  /*!
    @brief Fill all solvers with initial data
  */
  virtual void initial_data();

  /*!
    @brief Regrid
  */
  virtual void regrid(const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Regrid internal storage. This MUST be overwritten.
  */
  virtual void regrid_internals(){
#if 0
    MayDay::Warning("time_stepper::regrid_internals - this routine should be PURE");
#endif
  }

  /*!
    @brief Regrid individual solvers
  */
  virtual void regrid_solvers(const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Do a sanity check
  */
  virtual void sanity_check();

    /*!
    @brief Set amr
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Set the physical domain
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom);

  /*!
    @brief Set hardcap on the time step
  */
  virtual void set_min_dt(const Real a_min_dt);

    /*!
    @brief Set hardcap on the time step
  */
  virtual void set_max_dt(const Real a_min_dt);

  /*!
    @brief Set CFL condition
  */
  virtual void set_cfl(const Real a_cfl);

  /*!
    @brief Set fudge factor for relaxation time
  */
  virtual void set_relax_time(const Real a_relax_time);

  /*!
    @brief Set the physical domain
  */
  virtual void set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom);

  /*!
    @brief Set plasma kinetics
  */
  virtual void set_plasma_kinetics(const RefCountedPtr<plasma_kinetics>& a_plaskin);

  /*!
    @brief Verbosity
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Set solver verbosity
  */
  virtual void set_solver_verbosity(const int a_verbosity);

  /*!
    @brief Setup the cdr solvers
  */
  virtual void setup_cdr();

  /*!
    @brief Setup the poisson solver
  */
  virtual void setup_poisson();

  /*!
    @brief Setup the rte solvers
  */
  virtual void setup_rte();

  /*!
    @brief Instantiate the sigma solver
  */
  virtual void setup_sigma();

  /*!
    @brief Set potential
  */
  virtual void set_potential(Real (*a_potential)(const Real a_time));

  /*!
    @brief Solve the poisson equation. This calls compute_rho() 
    @details This assumes that the data contain in cdr_layout is cell-centered
  */
  virtual void solve_poisson();

  /*!
    @brief General poisson solver routine. 
    @details You have the option of supplying storage for a right-hand side if you want to. The densities can be either cell-centered or centroid-centered, If they are cell-centered, we will extrapolate the right-hand side to the centroids first. 
  */
  virtual void solve_poisson(MFAMRCellData&                a_potential,
			     MFAMRCellData*                a_rhs,
			     const Vector<EBAMRCellData*>  a_densities,
			     const EBAMRIVData&            a_sigma,
			     const centering::which_center a_centering);

  /*!
    @brief RTE solve using data existing in other solvers
  */
  virtual void solve_rte(const Real a_dt);

  /*!
    @brief General RTE solve routine
  */
  virtual void solve_rte(Vector<EBAMRCellData*>&       a_states,
			 Vector<EBAMRCellData*>&       a_rhs,
			 const Vector<EBAMRCellData*>& a_cdr_states,
			 const EBAMRCellData&          a_E,
			 const Real                    a_dt,
			 const centering::which_center a_centering);

  /*!
    @brief Synchronize solver times
  */
  virtual void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Query about necessary number of ghost cells. This collects that from the individual solvers. 
  */
  virtual int query_ghost();

  /*!
    @brief Compute the relaxation time
  */
  virtual Real compute_relaxation_time();

  /*!
    @brief Restrict dt
  */
  virtual Real restrict_dt(){
#if 0
    MayDay::Warning("time_stepper::restrict_dt - this routine should be PURE");
#endif
  }
  
protected:

  /*!
    @brief Index space
  */
  RefCountedPtr<mfis> m_mfis;

  /*!
    @brief AMR
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Plasma kinetics
  */
  RefCountedPtr<plasma_kinetics> m_plaskin;

  /*!
    @brief Physical domain
  */
  RefCountedPtr<physical_domain> m_physdom;

  /*!
    @brief Physical domain
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief CDR solvers
  */
  RefCountedPtr<cdr_layout> m_cdr;

  /*!
    @brief Radiative transfer solvers
  */
  RefCountedPtr<rte_layout> m_rte;

  /*!
    @brief Poisson solver
  */
  RefCountedPtr<poisson_solver> m_poisson;

  /*!
    @brief Surface charge solver
  */
  RefCountedPtr<sigma_solver> m_sigma;

  /*!
    @brief Set the potential
  */
  Real (*m_potential)(const Real a_time);
  
  /*!
    @brief Minimum allowed time step
  */
  Real m_min_dt;

  /*!
    @brief Maximum allowed time step
  */
  Real m_max_dt;

  /*!
    @brief CFL
  */
  Real m_cfl;

  /*!
    @brief Fudge factor for relaxation time
  */
  Real m_relax_time;

  /*!
    @brief Verbosity
  */
  int m_verbosity;
};

#endif
