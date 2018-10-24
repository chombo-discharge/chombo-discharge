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
  @brief Abstract class for evolving in time across an amr mesh. 
  @details This class has default options relating
  
  Class options
  -------------

  time_stepper.verbosity        = 1     # Class verbosity
  time_stepper.solver_verbosity = 0     # Individual solver verbosities
  time_stepper.fast_rte         = 1     # Solve RTE every this time steps
  time_stepper.fast_poisson     = 1     # Solve Poisson every this time steps
  time_stepper.min_dt           = 0.    # Minimum permitted time step
  time_stepper.min_dt           = 1.    # Maximum permitted time step
  time_stepper.cfl              = 0.5   # CFL number
  time_stepper.relax_time       = 1.0   # Relaxation time constant
  time_stepper.source_growth    = 100.0 # Relaxation time constant

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
  virtual Real advance(const Real a_dt) = 0;

  /*!
    @brief Check if RTE solvers are stationary
  */
  virtual bool stationary_rte();

  /*!
    @brief Cache states
  */
  virtual void cache_states();

  /*!
    @brief Solve the poisson equation. This calls compute_rho() 
    @details This assumes that the data contain in cdr_layout is cell-centered
  */
  virtual bool solve_poisson();

  /*!
    @brief General poisson solver routine. 
    @details You have the option of supplying storage for a right-hand side if you want to. The densities can be either cell-centered or centroid-centered, If they are cell-centered, we will extrapolate the right-hand side to the centroids first. 
  */
  virtual bool solve_poisson(MFAMRCellData&                a_potential,
			     MFAMRCellData&                a_rhs,
			     const Vector<EBAMRCellData*>  a_densities,
			     const EBAMRIVData&            a_sigma,
			     const centering::which_center a_centering);

  /*!
    @brief Compute the time step and how it was restricted
  */
  virtual void compute_dt(Real& a_dt, time_code::which_code& a_timecode);

  /*!
    @brief Compute the cell-centered electric field on both phases
  */
  virtual void compute_E(MFAMRCellData& a_E, const MFAMRCellData& a_potential);

  /*!
    @brief Compute the cell-centered electric field on a specific phase using whatever is available in m_poisson
  */
  virtual void compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase);

  /*!
    @brief Compute the cell-centered electric field on a specific phase
  */
  virtual void compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase, const MFAMRCellData& a_potential);

  /*!
    @brief Compute the face-centered electric field (only normal components)
  */
  virtual void compute_E(EBAMRFluxData& a_E_face, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell);
  
  /*!
    @brief Compute the face-centered electric field (only normal components)
  */
  virtual void compute_E(EBAMRIVData& a_E_eb, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell);

  /*!
    @brief Compute the maximum of the electric field
  */
  virtual void compute_Emax(Real& a_Emax, const phase::which_phase a_phase);

  /*!
    @brief Compute the current density
  */
  virtual void compute_J(EBAMRCellData& a_J);
  
  /*!
    @brief Compute diffusion things using whatever is available in the solvers
  */
  virtual void compute_cdr_diffusion();

  /*!
    @brief Compute diffusion things using whatever is available in the solvers, except for the fields
  */
  virtual void compute_cdr_diffusion(const EBAMRCellData& a_E_cell, const EBAMRIVData& a_E_eb);

  /*!
    @brief Compute diffusion coefficients
    @details The way we do this is that we compute the diffusion coefficients on cell centers, and then we average that
    onto cell faces. 
  */
  virtual void compute_cdr_diffco_face(Vector<EBAMRFluxData*>&       a_diffco_face,
				       const Vector<EBAMRCellData*>& a_cdr_densities,
				       const EBAMRCellData&          a_E,
				       const Real&                   a_time);

  /*!
    @brief Compute diffusion coefficients on the EB. 
  */
  virtual void compute_cdr_diffco_eb(Vector<EBAMRIVData*>&       a_diffco_face,
				     const Vector<EBAMRIVData*>& a_cdr_densities,
				     const EBAMRIVData&          a_E,
				     const Real&                 a_time);

  /*!
    @brief Compute cdr_fluxes
  */
  virtual void compute_cdr_fluxes(Vector<EBAMRIVData*>&       a_fluxes,
				  const Vector<EBAMRIVData*>& a_extrap_cdr_fluxes,
				  const Vector<EBAMRIVData*>& a_extrap_cdr_densities,
				  const Vector<EBAMRIVData*>& a_extrap_cdr_velocities,
				  const Vector<EBAMRIVData*>& a_extrap_cdr_gradients,
				  const Vector<EBAMRIVData*>& a_extrap_rte_fluxes,
				  const EBAMRIVData&          a_field,
				  const Real&                 a_time);

  /*!
    @brief Compute cdr fluxes on domain faces
    @note Everything that comes in here should have been computed on the domain faces
  */
  virtual void compute_cdr_domain_fluxes(Vector<EBAMRIFData*>&       a_fluxes,
					 const Vector<EBAMRIFData*>& a_extrap_cdr_fluxes,
					 const Vector<EBAMRIFData*>& a_extrap_cdr_densities,
					 const Vector<EBAMRIFData*>& a_extrap_cdr_velocities,
					 const Vector<EBAMRIFData*>& a_extrap_cdr_gradients,
					 const Vector<EBAMRIFData*>& a_extrap_rte_fluxes,
					 const EBAMRIFData&          a_field,
					 const Real&                 a_time);

  /*!
    @brief Compute the gradient at the boundary (we also take the dot product w.r.t. to the normal vector)
    @note This version computes the gradients from phi
  */
  virtual void compute_gradients_at_eb(Vector<EBAMRIVData*>&         a_grad,
				       const phase::which_phase&     a_phase,
				       const Vector<EBAMRCellData*>& a_phi);

  /*!
    @brief Compute gradients at domain faces
  */
  virtual void compute_gradients_at_domain_faces(Vector<EBAMRIFData*>&         a_grad,
						 const phase::which_phase&     a_phase,
						 const Vector<EBAMRCellData*>& a_phi);

  /*!
    @brief Compute CDR source terms
  */
  virtual void compute_cdr_sources();


  /*!
    @brief Compute cdr source terms
    @note This version computes the cdr_gradients and then calls the full compute_cdr_sources
  */
  virtual void compute_cdr_sources(Vector<EBAMRCellData*>&        a_sources,
				   const Vector<EBAMRCellData*>&  a_cdr_densities,
				   const Vector<EBAMRCellData*>&  a_rte_densities,
				   const EBAMRCellData&           a_E,
				   const Real&                    a_time,
				   const centering::which_center  a_centering);

  /*!
    @brief Compute cdr source terms
  */
  virtual void compute_cdr_sources(Vector<EBAMRCellData*>&        a_sources,
				   const Vector<EBAMRCellData*>&  a_cdr_densities,
				   const Vector<EBAMRCellData*>&  a_cdr_gradients,
				   const Vector<EBAMRCellData*>&  a_rte_densities,
				   const EBAMRCellData&           a_E,
				   const Real&                    a_time,
				   const centering::which_center  a_centering);

  /*!
    @brief Compute cdr velocities
  */
  virtual void compute_cdr_velocities();

  /*!
    @brief Compute the cell-centered cdr velocities
  */
  virtual void compute_cdr_velocities(Vector<EBAMRCellData*>&       a_velocities,
				      const Vector<EBAMRCellData*>& a_cdr_densities,
				      const EBAMRCellData&          a_E,
				      const Real&                   a_time);


  /*!
    @brief Compute photon source terms. 
  */
  virtual void compute_rte_sources(Vector<EBAMRCellData*>        a_source,
				   const Vector<EBAMRCellData*>& a_cdr_states,
				   const EBAMRCellData&          a_E,
				   const Real&                   a_time,
				   const centering::which_center a_centering);

  /*!
    @brief Compute the cell-centered fluxes and extrapolate these to the EB centroid
  */
  virtual void compute_extrapolated_fluxes(Vector<EBAMRIVData*>&        a_fluxes,
					   const Vector<EBAMRCellData*> a_densities,
					   const Vector<EBAMRCellData*> a_velocities,
					   const phase::which_phase     a_phase);

  /*!
    @brief Compute flux
  */
  virtual void compute_flux(EBAMRCellData&       a_flux,
			    const EBAMRCellData& a_density,
			    const EBAMRCellData& a_velocity);

  /*!
    @brief Compute charge flux
  */
  virtual void compute_charge_flux(EBAMRIVData& a_flux, Vector<EBAMRIVData*>& a_cdr_flxues);

  /*!
    @brief Compute RTE source terms
  */
  virtual void compute_rte_sources();



  /*!
    @brief Compute the cell-centered space charge density by using data inside cdr_layout. 
    @details This assumes that the data contain in cdr_layout is cell-centered
  */
  virtual void compute_rho();

  /*!
    @brief Compute cell-centered space charge density on a specific phase using whatever is available in cdr_layout
  */
  virtual void compute_rho(EBAMRCellData& a_rho, const phase::which_phase a_phase);

  /*!
    @brief Compute the centroid-centered space charge density
    @details You have the option of specifying the centering for this computation. If the densities are cell-centered, we first 
    compute the cell-centered space charge density and the interpolate to the cell centroid. 
  */
  virtual void compute_rho(MFAMRCellData&                 a_rho,
			   const Vector<EBAMRCellData*>&  a_densities,
			   const centering::which_center  a_centering);

  /*!
    @brief Deallocate internals
  */
  virtual void deallocate_internals(){}
  
  /*!
    @brief Deallocate internal solver storages
  */
  virtual void deallocate_solver_internals();

  /*!
    @brief Extrapolate data to EB centroid
  */
  virtual void extrapolate_to_eb(EBAMRIVData& a_extrap, const phase::which_phase a_phase, const EBAMRCellData& a_data);

  /*!
    @brief Extrapolate data to EB centroid
  */
  virtual void extrapolate_to_eb(Vector<EBAMRIVData*>&         a_extrap,
				 const phase::which_phase      a_phase,
				 const Vector<EBAMRCellData*>& a_data);

  /*!
    @brief Extrapolate to domain faces
  */
  virtual void extrapolate_to_domain_faces(EBAMRIFData& a_extrap, const phase::which_phase a_phase, const EBAMRCellData& a_data);

  /*!
    @brief Extrapolate to domain faces
  */
  virtual void extrapolate_to_domain_faces(Vector<EBAMRIFData*>&         a_extrap,
					   const phase::which_phase      a_phase,
					   const Vector<EBAMRCellData*>& a_data);

  /*!
    @brief Get maximum value of the cdr solver stuff
  */
  virtual void get_cdr_max(Real& a_cdr_max, std::string& a_solver_name);

  /*!
    @brief Instantiate solvers
  */
  virtual void instantiate_solvers();

  /*!
    @brief Fill all solvers with initial data
  */
  virtual void initial_data();

  /*!
    @brief Fill initial cdr data
  */
  virtual void initial_cdr_data();

  /*!
    @brief Fill initial rte data
  */
  virtual void initial_rte_data();

  /*!
    @brief Fill initial sigma data
  */
  virtual void initial_sigma_data();

  /*!
    @brief Extrapolated flux to EB and project it onto the normal
  */
  virtual void project_flux(EBAMRIVData& a_projected_flux, const EBAMRIVData& a_flux);

  /*!
    @brief Project flux on the boundary normal
  */
  virtual void project_domain(EBAMRIFData& a_projected_flux, const EBAMRIFData& a_flux);

  /*!
    @brief Regrid
  */
  virtual void regrid(const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Regrid internal storage. This MUST be overwritten.
  */
  virtual void regrid_internals() = 0;

  /*!
    @brief Regrid individual solvers
  */
  virtual void regrid_solvers(const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Set data defined over dielectric cells to zero
  */
  virtual void reset_dielectric_cells(EBAMRIVData& a_data);

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
    @brief Set fast rte
  */
  virtual void set_fast_rte(const int a_fast_rte);

  /*!
    @brief Set fast rte
  */
  virtual void set_fast_poisson(const int a_fast_poisson);

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
    @brief Set valid levels for relaxation time
  */
  virtual void set_relax_level(const Real a_relax_level);

  /*!
    @brief Set source growth
  */
  virtual void set_source_growth(const Real a_src_growth);

  /*!
    @brief Set source growth
  */
  virtual void set_source_growth_tolerance(const Real a_src_tolerance);

  /*!
    @brief Set source growth
  */
  virtual void set_source_growth_elec_only(const bool a_src_elec_only);

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
    @brief Set stationary/nonstationary RTE
  */
  virtual void set_stationary_rte(const bool a_stationary);

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
    @brief Solver dumps data files
  */
  virtual void solver_dump();

  /*!
    @brief Set potential
  */
  virtual void set_potential(Real (*a_potential)(const Real a_time));

  /*!
    @brief Set wall function for Poisson equation boundary conditions
  */
  virtual void set_poisson_wall_func(const int a_dir, const Side::LoHiSide a_side, Real (*a_func)(const RealVect a_pos));

  /*!
    @brief RTE advance using data existing in other solvers
  */
  virtual void solve_rte(const Real a_dt);

  /*!
    @brief General RTE advance routine
    @details This computes the source terms at the specified time and then advances the solution
  */
  virtual void solve_rte(Vector<EBAMRCellData*>&       a_states,
			 Vector<EBAMRCellData*>&       a_rhs,
			 const Vector<EBAMRCellData*>& a_cdr_states,
			 const EBAMRCellData&          a_E,
			 const Real                    a_time,
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
    @brief Compute the ohmic current on electrodes
  */
  virtual Real compute_electrode_current();

  /*!
    @brief Compute the ohmic current on dielectrics
  */
  virtual Real compute_dielectric_current();

  /*!
    @brief Compute induced current in external circuit due to ohmic conduction
  */
  virtual Real compute_ohmic_induction_current();
			     
  /*!
    @brief Compute the relaxation time
  */
  virtual Real compute_relaxation_time();

  /*!
    @brief Restrict dt
  */
  virtual Real restrict_dt() = 0;

  /*!
    @brief Return time
  */
  virtual Real get_time();

  /*!
    @brief Get dt
  */
  virtual Real get_dt();

  /*!
    @brief Get cdr solvers
  */
  virtual RefCountedPtr<cdr_layout>& get_cdr();

  /*!
    @brief Get cdr solvers
  */
  virtual RefCountedPtr<poisson_solver>& get_poisson();
  
  /*!
    @brief Get cdr solvers
  */
  virtual RefCountedPtr<rte_layout>& get_rte();

  /*!
    @brief Get sigma solver
  */
  virtual RefCountedPtr<sigma_solver>& get_sigma();
  
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
    @brief Function pointer for Poisson BC on domain walls
  */
  Real (*m_wall_func_x_lo)(const RealVect a_pos);

  /*!
    @brief Function pointer for Poisson BC on domain walls
  */
  Real (*m_wall_func_x_hi)(const RealVect a_pos);

  /*!
    @brief Function pointer for Poisson BC on domain walls
  */
  Real (*m_wall_func_y_lo)(const RealVect a_pos);

  /*!
    @brief Function pointer for Poisson BC on domain walls
  */
  Real (*m_wall_func_y_hi)(const RealVect a_pos);

#if CH_SPACEDIM==3
  /*!
    @brief Function pointer for Poisson BC on domain walls
  */
  Real (*m_wall_func_z_lo)(const RealVect a_pos);

  /*!
    @brief Function pointer for Poisson BC on domain walls
  */
  Real (*m_wall_func_z_hi)(const RealVect a_pos);
#endif
  
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
    @brief Source growth
  */
  Real m_src_growth;

  /*!
    @brief TIme
  */
  Real m_time;

  /*!
    @brief Previous time step size
  */
  Real m_dt;

  /*!
    @brief Tolerance for source term growth
  */
  Real m_src_tolerance;

  /*!
    @brief Finest level on which to impose relaxation time constraints
  */
  int m_relax_level;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Verbosity for solvers
  */
  int m_solver_verbosity;

  /*!
    @brief Fast RTE solver
  */
  int m_fast_rte;

  /*!
    @brief Fast RTE solver
  */
  int m_fast_poisson;

  /*!
    @brief Current step
  */
  int m_step;

  /*!
    @brief Do source term time step restriction only for electrons
    @details If this one is true, compute_dt passes information down the toolchain to only compute the source term
    restriction time for the first species defined in the plasma kinetics
  */
  bool m_src_elec_only;
};

#endif
