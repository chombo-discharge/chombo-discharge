/*!
  @file cdr_solver.H
  @brief Declaration of an abstract class for evolving cdr equations
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _CDR_SOLVER_
#define _CDR_SOLVER_

#include "type_definitions.H"
#include "computational_geometry.H"
#include "amr_mesh.H"
#include "physical_domain.H"
#include "species.H"

/*!
  @brief Abstract class for cdr-equations. This is (almost) a whole solver. 
*/
class cdr_solver {
public:

  /*!
    @brief Default constructor. Stay away.
  */
  cdr_solver();

  /*!
    @brief Constructor
  */
  virtual ~cdr_solver();
  
  /*!
    @brief Cache internal state
  */
  virtual void cache_state();

    /*!
    @brief Compute div(nv - D*grad(n))
  */
  virtual void compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt) = 0;

  /*!
    @brief Compute advective derivative
  */
  virtual void compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt, const bool a_redist) = 0;

  /*!
    @brief Compute diffusion term
  */
  virtual void compute_divD(EBAMRCellData& a_diffusive_term, const EBAMRCellData& a_state) = 0;

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocate_internals();

  /*!
    @brief Sanity check
  */
  virtual void sanity_check();

  /*!
    @brief Set species
  */
  virtual void set_species(const RefCountedPtr<species> a_species);
  
  /*!
    @brief Set computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom);

  /*!
    @brief Set the physical domain
  */
  virtual void set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom);

  /*!
    @brief Set the amr object
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng set_computational_geometry
  */
  virtual void set_phase(phase::which_phase a_phase = phase::gas);

  /*!
    @brief Set verbosity.
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Get current time
  */
  virtual void set_time(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Set velocity
  */
  virtual void set_velocity(const EBAMRCellData& a_velo);

  /*!
    @brief Set constant velocity
  */
  virtual void set_velocity(const RealVect a_velo);

  /*!
    @brief Set diffusion coefficients
  */
  virtual void set_diffco(const EBAMRFluxData& a_diffco, const EBAMRIVData& a_diffco_eb);

  /*!
    @brief Set constant
  */
  virtual void set_diffco(const Real a_diffco);

  /*!
    @brief Set source terms
  */
  virtual void set_source(const EBAMRCellData& a_source);

  /*!
    @brief Set constant source terms
  */
  virtual void set_source(const Real a_source);

  /*!
    @brief Set the EB flux
  */
  virtual void set_ebflux(const EBAMRIVData& a_ebflux);

  /*!
    @brief Set the EB flux
  */
  virtual void set_ebflux(const Real a_ebflux);

  /*!
    @brief Fill internal state with initial data
  */
  virtual void initial_data();

  /*!
    @brief Allocate internal storage
  */
  virtual void allocate_internals();
  
  /*!
    @brief Advance method
  */
  virtual void advance(const Real& a_dt);

  /*!
    @brief Advance state. Default is second order Runge-Kutta.
  */
  virtual void advance(EBAMRCellData& a_state, const Real& a_dt);

  /*!
    @brief Advective advance on internal data holder
  */
  virtual void advance_advect(const Real a_dt);

  /*!
    @brief Pure advective advance. I.e. advance with D = 0 and S = 0
  */
  virtual void advance_advect(EBAMRCellData& a_state, const Real a_dt) = 0;
  
  /*!
    @brief Diffusion advance on internal data holder
  */
  virtual void advance_diffusion(const Real a_dt);

  /*!
    @brief Pure diffusion advance. I.e. advance with div(nv) = 0 and S = 0
  */
  virtual void advance_diffusion(EBAMRCellData& a_state, const Real a_dt) = 0;

  /*!
    @brief Explicit RK2 advance
  */
  virtual void advance_rk2(EBAMRCellData& a_state, const Real a_dt, const Real a_alpha);

  /*!
    @brief Compute right-hand side using the internal state.
  */
  virtual void compute_rhs(EBAMRCellData& a_rhs, const Real& a_dt); 

  /*!
    @brief Compute the right-hand side using supplied solution. This is the basic function that is required for using the method
    of lines (MOL)
  */
  virtual void compute_rhs(EBAMRCellData& a_rhs, const EBAMRCellData& a_state, const Real& a_dt);

#ifdef CH_USE_HDF5
  /*!
    @brief Write plot file
  */
  virtual void write_plot_file();
#endif

  /*!
    @brief Use mass-weighted redistribution
  */
  virtual void set_mass_redist(const bool a_mass_redist);

  /*!
    @brief Regrid internals
  */
  virtual void regrid(const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Get solver name
  */
  virtual std::string get_name();

  /*!
    @brief Get the ghost cell requirement
  */
  virtual int query_ghost() const = 0;

  /*!
    @brief Compute the largest possible advection time step (for explicit methods)
  */
  virtual Real compute_cfl_dt();

  /*!
    @brief Compute the largest possible diffusive time step (for explicit methods)
  */
  virtual Real compute_diffusive_dt();

  /*!
    @brief Compute the largest possible source time step (for explicit methods
  */
  virtual Real compute_source_dt(const Real a_max, const Real a_tolerance);
  
  /*!
    @brief Compute the total mass in m_state
  */
  virtual Real compute_mass();

  /*!
    @brief Compute the total charge in m_state
  */
  virtual Real compute_charge();

  /*!
    @brief Diffusive or not
  */
  virtual bool is_diffusive();

  /*!
    @brief Mobile or not
  */
  virtual bool is_mobile();
  
  /*!
    @brief Get the internal state
  */
  virtual EBAMRCellData& get_state();

  /*!
    @brief Get the source term
  */
  virtual EBAMRCellData& get_source();

  /*!
    @brief Get the cell-centered velocity
  */
  virtual EBAMRCellData& get_velo_cell();

  /*!
    @brief Get the face-centered velocities
  */
  virtual EBAMRFluxData& get_velo_face();

  /*!
    @brief Get the eb-centered velocities
  */
  virtual EBAMRIVData& get_velo_eb();

  /*!
    @brief Get the face-centered diffusion coefficient
  */
  virtual EBAMRFluxData& get_diffco_face();

  /*!
    @brief Get the EB-centered diffusion coefficient
  */
  virtual EBAMRIVData& get_diffco_eb();

  /*!
    @brief Get the eb flux data holder
  */
  virtual EBAMRIVData& get_ebflux();

protected:
  
  /*!
    @brief Phase
  */
  RefCountedPtr<EBIndexSpace> m_ebis;

  /*!
    @brief Species
  */
  RefCountedPtr<species> m_species;
  
  /*!
    @brief Computational geometry
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief Physical domain
  */
  RefCountedPtr<physical_domain> m_physdom;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Stencils for interpolating flux to centroids
  */
  Vector<RefCountedPtr<LayoutData<BaseIFFAB<FaceStencil> > > > m_interp_stencils[SpaceDim];

  /*!
    @brief Stencils for computing the non-conservative divergence
  */
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > > m_stencils_nc;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Name
  */
  std::string m_name;

  /*!
    @brief State data (i.e. the ion density)
  */
  EBAMRCellData m_state;

  /*!
    @brief Source term 
  */
  EBAMRCellData m_source;

  /*!
    @brief Cell-centered velocities
  */
  EBAMRCellData m_velo_cell;

  /*!
    @brief Scratch storage. 
  */
  EBAMRCellData m_scratch;

  /*!
    @brief Cached internal state. DO NOT TOUCH!
  */
  EBAMRCellData m_cache;

  /*!
    @brief Face-centered velocities (only normal components)
  */
  EBAMRFluxData m_velo_face;

  /*!
    @brief EB-centered velocities
  */
  EBAMRIVData m_velo_eb;

  /*!
    @brief Flux through the boundary
  */
  EBAMRIVData m_ebflux;

  /*!
    @brief Diffusion coefficients (can
  */
  EBAMRFluxData m_diffco;

  /*!
    @brief Diffusion coefficients (can
  */
  EBAMRIVData m_diffco_eb;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_step;

  /*!
    @brief Current time
  */
  Real m_time;

  /*!
    @brief Time step
  */
  Real m_dt;

  /*!
    @brief Also solve for diffusion
  */
  bool m_diffusive;

  /*!
    @brief Also solve for advection/convection
  */
  bool m_mobile;

  /*!
    @brief Use mass-weighted redistribution
  */
  bool m_mass_redist;

  /*!
    @brief Average cell-centered velocities to face centers
  */
  virtual void average_velo_to_faces(EBAMRFluxData& a_velo_face, const EBAMRCellData& a_velo_cell);


  /*!
    @brief Compute the conservative divergence using face-centered states and velocities
  */
  virtual void conservative_divergence(EBAMRCellData&       a_cons_div,
				       const EBAMRFluxData& a_face_vel,
				       const EBAMRFluxData& a_face_state);
  
  /*!
    @brief Set up flux
  */
  virtual void compute_flux(EBAMRFluxData& a_flux, const EBAMRFluxData& a_face_state, const EBAMRFluxData& a_face_vel);

  /*!
    @brief Compute conservative divergence from fluxes. 
  */
  virtual void conservative_divergence(EBAMRCellData& a_cons_div, const EBAMRFluxData& a_flux);

  /*!
    @brief Conservative divergence over regular cells
  */
  virtual void consdiv_regular(LevelData<EBCellFAB>& a_divJ, const LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Compute flux interpolant
  */
  virtual void setup_flux_interpolant(LevelData<BaseIFFAB<Real> >   a_interpolant[SpaceDim],
				      const LevelData<EBFluxFAB>&   a_flux,
				      const int                     a_lvl);

  /*!
    @brief Interpolate flux to centroids
  */
  virtual void interpolate_flux_to_centroids(LevelData<BaseIFFAB<Real> >       a_flux[SpaceDim],
					     const int                         a_lvl);

  /*!
    @brief Compute div(F) on irregular cells
  */
  virtual void compute_divF_irreg(LevelData<EBCellFAB>&              a_divF,
				  const LevelData<BaseIFFAB<Real> >  a_flux[SpaceDim],
				  const LevelData<BaseIVFAB<Real> >& a_ebflux,
				  const int                          a_lvl);
				  
  /*!
    @brief Compute the non-conservative divergence
  */
  virtual void nonconservative_divergence(EBAMRIVData& a_div_nc, const EBAMRCellData& a_divF, const EBAMRFluxData& a_face_state);


  /*!
    @brief Use the non-conservative divergence to make the conservative divergence hold the hybrid divergence
  */
  virtual void hybrid_divergence(EBAMRCellData&       a_divF,
				 EBAMRIVData&         a_mass_diff,
				 const EBAMRIVData&   a_noncons_div);

  /*!
    @brief Increment flux registers. 
  */
  virtual void increment_flux_register(const EBAMRFluxData& a_face_state, const EBAMRFluxData& a_velo_face);

  /*!
    @brief Increment flux registers
  */
  virtual void increment_flux_register(const EBAMRFluxData& a_flux);

  /*!
    @brief Increment the coarse-fine redist with mass difference (only for EB x CF)
  */
  virtual void coarse_fine_increment(const EBAMRIVData& m_mass_diff);

  /*!
    @brief Do mass redistribution
  */
  virtual void hyperbolic_redistribution(EBAMRCellData&       a_del_vel_rho,
					 const EBAMRIVData&   a_mass_diff,
					 const EBAMRCellData& a_redist_weights);

  /*!
    @brief Increment redistribution flux registers
  */
  virtual void increment_redist_flux();

  /*!
    @brief Increment redistribution objects with mass to be redistributed
  */
  virtual void increment_redist(const EBAMRIVData& a_mass_diff);

  /*!
    @brief Do coarse-fine redistribution
  */
  virtual void coarse_fine_redistribution(EBAMRCellData& a_state);

  /*!
    @brief Reflux state
  */
  virtual void reflux(EBAMRCellData& a_state);

  /*!
    @brief Set ebis
  */
  virtual void set_ebis(const RefCountedPtr<EBIndexSpace>& a_ebis);

  /*!
    @brief Define stencils
  */
  virtual void define_interp_stencils();

  /*!
    @brief Define stencils for computing the non-conservative divergence
  */
  virtual void define_divFnc_stencils();
};

#endif
