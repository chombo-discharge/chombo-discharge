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
#include "ebconductivityopfactory.H"

#include <AMRMultiGrid.H>
#include <AMRTGA.H>
#include <BiCGStabSolver.H>
#include <EBBackwardEuler.H>
#include <EBSimpleSolver.H>

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
    @brief Set solver parameters. These are normally good for this class. 
  */
  virtual void set_gmg_solver_parameters(relax::which_relax a_relax_type = relax::gsrb_fast,
					 amrmg::which_mg a_gmg_type      = amrmg::vcycle,
					 const int a_verbosity           = -1,
					 const int a_pre_smooth          = 16,
					 const int a_post_smooth         = 16,
					 const int a_bot_smooth          = 16,
					 const int a_max_iter            = 16,
					 const int a_min_iter            = 5,
					 const Real a_eps                = 1.E-8,
					 const Real a_hang               = 0.2);

  /*!
    @brief Set the bottom solver.  
    @param[in] a_which_solver Solver type. 0 = smoother, 1 = BiCGStab
  */
  virtual void set_bottom_solver(const int a_whichsolver);

  /*!
    @brief Set the depth for bottom drop
  */
  virtual void set_bottom_drop(const int a_bottom_drop);

  /*!
    @brief Set number of smoothing for the bottom solver (simplified solver only)
  */
  virtual void set_botsolver_smooth(const int a_smooth);

  /*!
    @brief Use tga for temporal integration
  */
  virtual void set_tga(const bool a_use_tga);

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
    @brief Advance state. Current is explicit Euler.
    @details Make the default RK2 for non-diffusive media and semi-implicit TGA for diffusive media. 
  */
  virtual void advance(EBAMRCellData& a_state, const Real& a_dt);

  /*!
    @brief Compute right-hand side using the internal state.
  */
  virtual void compute_rhs(EBAMRCellData& a_rhs, const Real& a_dt); 

  /*!
    @brief Compute the right-hand side using supplied solution
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
    @brief Compute the total mass in m_state
  */
  virtual Real compute_mass();

  /*!
    @brief Diffusive or not
  */
  virtual bool is_diffusive();
  
  /*!
    @brief Get the internal state
  */
  virtual EBAMRCellData& get_state();

  /*!
    @brief Get the source term
  */
  virtual EBAMRCellData& get_source();

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
    @brief Geometric multigrid solver
  */
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > m_gmg_solver;

  /*!
    @brief TGA solver
  */
  RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > > m_tgasolver;

  /*!
    @brief Euler solver
  */
  RefCountedPtr<EBBackwardEuler> m_eulersolver;

  /*!
    @brief Operator factory
  */
  RefCountedPtr<ebconductivityopfactory> m_opfact;

  /*!
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<EBCellFAB> > m_bicgstab;

  /*!
    @brief Simple solver
  */
  EBSimpleSolver m_simple_solver;

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
    @brief Relaxation type for gmg
  */
  relax::which_relax m_gmg_relax_type;

  /*!
    @brief GMG multigrid type
  */
  amrmg::which_mg m_gmg_type;

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
    @brief Verbosity for geometric multigrid
  */
  int m_gmg_verbosity;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_gmg_pre_smooth;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_gmg_post_smooth;

  /*!
    @brief Number of smoothing before bottom solver
  */
  int m_gmg_bot_smooth;

  /*!
    @brief Maximum number of iterations
  */
  int m_gmg_max_iter;

  /*!
    @brief Minimum number of iterations
  */
  int m_gmg_min_iter;

  /*!
    @brief Bottom solver type
  */
  int m_bottomsolver;

  /*!
    @brief Number of smoothing for bottom solver
  */
  int m_numsmooth;

  /*!
    @brief Set bottom drop depth
  */
  int m_bottom_drop;

  /*!
    @brief Use tga
  */
  bool m_use_tga;

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
    @brief 
  */
  Real m_gmg_eps;
  
  /*!
    @brief 
  */
  Real m_gmg_hang;


  /*!
    @brief Also solve for diffusion
  */
  bool m_diffusive;

  /*!
    @brief Use mass-weighted redistribution
  */
  bool m_mass_redist;

  /*!
    @brief Explicit RK2 advance
  */
  virtual void advance_rk2(EBAMRCellData& a_state, const Real a_dt, const Real a_alpha);

  /*!
    @brief TGA advancement algorithm
  */
  virtual void advance_tga(EBAMRCellData& a_state, const Real a_dt);

  /*!
    @brief Compute div(nv - D*grad(n))
  */
  virtual void compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt);

  /*!
    @brief Compute advective derivative
  */
  virtual void compute_advection_term(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt);

  /*!
    @brief Compute diffusion term
  */
  virtual void compute_diffusion_term(EBAMRCellData& a_diffusive_term, const EBAMRCellData& a_state);

  /*!
    @brief Average cell-centered velocities to face centers
  */
  virtual void average_velo_to_faces(EBAMRFluxData& a_velo_face, const EBAMRCellData& a_velo_cell);

  /*!
    @brief Extrapolate state to faces. 
  */
  virtual void extrapolate_to_faces(EBAMRFluxData& a_face_state, const EBAMRCellData& a_state, const Real a_extrap_dt) = 0;

  /*!
    @brief Compute the conservative divergence using face-centered states and velocities
  */
  virtual void conservative_divergence(EBAMRCellData&       a_cons_div,
				       const EBAMRFluxData& a_face_vel,
				       const EBAMRFluxData& a_face_state);

  /*!
    @brief Compute the advective derivative. This computes the divF = sum_i(vel*phi) using cell centered data.
  */
  virtual void advective_derivative(LevelData<EBCellFAB>&       a_divF,
				    const LevelData<EBFluxFAB>& a_vel,
				    const LevelData<EBFluxFAB>& a_phi,
				    const int                   a_lvl);

  virtual void compute_flux_interpolant(LevelData<BaseIFFAB<Real> >   a_interpolant[SpaceDim],
					const LevelData<EBFluxFAB>&   a_rho,
					const LevelData<EBFluxFAB>&   a_vel,
					const int                     a_lvl);

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

  /*!
    @brief Set up geometric multigrid
  */
  virtual void setup_gmg();

  /*!
    @brief Setup the operator factory
  */
  virtual void setup_operator_factory();

  /*!
    @brief Setup multigrid
  */
  virtual void setup_multigrid();

  /*!
    @brief Setup the TGA solver
  */
  virtual void setup_tga();

    /*!
    @brief Setup the Euler solver
  */
  virtual void setup_euler();
};

#endif
