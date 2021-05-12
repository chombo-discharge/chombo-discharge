/*!
  @file cdr_solver.H
  @brief Declaration of an abstract class for evolving cdr equations
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _CDR_SOLVER_
#define _CDR_SOLVER_

#include "computational_geometry.H"
#include "CD_AmrMesh.H"
#include "cdr_species.H"

#include <random>
#include <time.h>
#include <chrono>

#include "CD_NamespaceHeader.H"

/*!
  @brief Boundary condition types for the CDR solvers
*/
namespace cdr_bc {
  enum which_bc {
    external = 0,
    wall     = 1,
    outflow  = 2,
    extrap   = 3,
  };
};

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
    @brief Parse class options from the command line
  */
  virtual void parseOptions() = 0;

  /*!
    @brief Parse runtime options
  */
  virtual void parseRuntimeOptions() = 0;

  /*!
    @brief Implicit diffusion TGA advance without source terms
  */
  virtual void advance_tga(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt) = 0;

  /*!
    @brief Implicit diffusion TGA advance with source terms
  */
  virtual void advance_tga(EBAMRCellData&       a_new_state,
			   const EBAMRCellData& a_old_state,
			   const EBAMRCellData& a_source,
			   const Real           a_dt) = 0;

  /*!
    @brief Implicit diffusion Euler advance without source terms
  */
  virtual void advance_euler(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt) = 0;

  /*!
    @brief Implicit diffusion Euler advance with source terms
  */
  virtual void advance_euler(EBAMRCellData&       a_new_state,
			     const EBAMRCellData& a_old_state,
			     const EBAMRCellData& a_source,
			     const Real           a_dt) = 0;

  /*!
    @brief Compute div(nv - D*grad(n)) explicitly
    @note a_state is non-const because ghost cells will be re-filled
  */
  virtual void compute_divJ(EBAMRCellData& a_divJ,
			    EBAMRCellData& a_state,
			    const Real     a_extrap_dt,
			    const bool     a_ebflux = true) = 0;

  /*!
    @brief Compute div(v*phi) explicitly
    @note a_state is non-const because ghost cells will be re-filled
  */
  virtual void compute_divF(EBAMRCellData&  a_divF,
			    EBAMRCellData&  a_state,
			    const Real      a_extrap_dt,
			    const bool      a_ebflux = true) = 0;

  /*!
    @brief Compute div(D*grad(phi)) explicitly
    @note a_state is non-const because ghost cells will be re-filled
  */
  virtual void compute_divD(EBAMRCellData& a_diffusive_term,
			    EBAMRCellData& a_state,
			    const bool     a_ebflux = false) = 0;

  /*!
    @brief Compute div(G) where G is a general face-centered flux on face centers and EB centers
    @param[in]    a_divG div(G). This includes redistribution
    @param[inout] a_G Must contain face-centered fluxes on input. Contains face-centroid fluxes on output. 
    @param[in]    a_ebG Flux on the EB centroids
    @note If you also want charge injection through the domain sides, you must put those in a_G before entering this routine
  */
  virtual void compute_divG(EBAMRCellData& a_divG, EBAMRFluxData& a_G, const EBAMRIVData& a_ebG);

  /*!
    @brief Random source for diffusion
  */
  virtual void GWN_diffusion_source(EBAMRCellData& a_ransource, const EBAMRCellData& a_cell_states);

  /*!
    @brief Smoothed face average
  */
  virtual void smooth_heaviside_faces(EBAMRFluxData& a_face_states, const EBAMRCellData& a_cell_states);

  /*!
    @brief Gaussian noise field
  */
  virtual void fill_GWN(EBAMRFluxData& a_noise, const Real a_mean);

  /*!
    @brief Average velocities to faces
  */
  virtual void average_velo_to_faces();

  /*!
    @brief Pre-regrid
  */
  virtual void pre_regrid(const int a_lbase, const int a_old_finest_level);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocate_internals();

  /*!
    @brief Register operators for AMR operations
  */
  virtual void registerOperators();

  /*!
    @brief Sanity check
  */
  virtual void sanityCheck();

  /*!
    @brief Set the realm
  */
  virtual void set_realm(const std::string a_realm);

  /*!
    @brief Set species
  */
  virtual void set_species(const RefCountedPtr<cdr_species> a_species);
  
  /*!
    @brief Set computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom);

  /*!
    @brief Set the amr object
  */
  virtual void set_amr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set domain bc type
  */
  virtual void set_domain_bc(const cdr_bc::which_bc a_bctype);

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
    @brief Set the domain flux
  */
  virtual void set_domain_flux(const Real a_domain_flux);

  /*!
    @brief Fill internal state with initial data
  */
  virtual void initial_data();

  /*!
    @brief Fill initial data from a distribution function
  */
  virtual void initial_data_distribution();

  /*!
    @brief Fill initial data from particles
  */
  virtual void initial_data_particles();

  /*!
    @brief Inject the EB flux. 
    @details This routine is equivalent to an update phi^(k+1) = phi^k - dt*div(G)
    where G is non-zero only on the EB
  */
  virtual void inject_ebflux(EBAMRCellData& a_phi, const EBAMRIVData& a_ebflux, const Real a_dt);

  /*!
    @brief Allocate internal storage
  */
  virtual void allocate_internals();

  /*!
    @brief Write plot file
  */
  virtual void write_plot_file();

  /*!
    @brief Write output data to a_output
  */
  virtual void write_plot_data(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write data to output. Convenience function. 
  */
  virtual void write_data(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp);

  /*!
    @brief Write checkpoint data into handle
  */
  virtual void write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void read_checkpoint_level(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Regrid internals
  */
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Set which variables to output
  */
  virtual void set_plot_variables();

  /*!
    @brief Reset redistribution weights
  */
  virtual void reset_redist_weights(const EBAMRCellData& a_state);

  /*!
    @brief Get the realm
  */
  virtual const std::string get_realm() const;

  /*!
    @brief Get solver name
  */
  virtual std::string get_name();

  /*!
    @brief Get output plot names
  */
  virtual Vector<std::string> get_plotvar_names() const;

  /*!
    @brief Get the ghost cell requirement
  */
  virtual int query_ghost() const = 0;

  /*!
    @brief Get number of output fields
  */
  virtual int get_num_plotvars() const;

  /*!
    @brief Compute the largest possible diffusive time step (for explicit methods)
  */
  virtual Real compute_advection_dt();

  /*!
    @brief Compute the largest possible diffusive time step (for explicit methods)
  */
  virtual Real compute_diffusion_dt();

  /*!
    @brief Compute the largest possible diffusive time step (for explicit methods)
  */
  virtual Real compute_advection_diffusion_dt();

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
    @brief Check if we can extrapolate
  */
  virtual bool extrap_source() const;
  
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

  /*!
    @brief Get the domain flux data holder
  */
  virtual EBAMRIFData& get_domainflux();


protected:
  
  /*!
    @brief Phase
  */
  RefCountedPtr<EBIndexSpace> m_ebis;

  /*!
    @brief Species
  */
  RefCountedPtr<cdr_species> m_species;
  
  /*!
    @brief Computational geometry
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Stencils for interpolating flux to centroids
  */
  Vector<RefCountedPtr<LayoutData<BaseIFFAB<FaceStencil> > > > m_interp_stencils[SpaceDim];

  /*!
    @brief Stencils for computing the non-conservative divergence
  */
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > > m_stencils_nc;

  /*!
    @brief Storage for flux interpolant
  */
  Vector<RefCountedPtr<LevelData<BaseIFFAB<Real> > > > m_interpolant[SpaceDim];

  /*!
    @brief Cells on which to compute the flux interpolant
  */
  Vector<RefCountedPtr<LayoutData<IntVectSet> > > m_interp_sets[SpaceDim];

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Domain bc condition
  */
  cdr_bc::which_bc m_dombc;

  /*!
    @brief Name
  */
  std::string m_name;

  /*!
    @brief Class name. 
  */
  std::string m_class_name;

  /*!
    @brief Realm
  */
  std::string m_realm;

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
    @brief Holder for face centered states
  */
  EBAMRFluxData m_face_states;

  /*!
    @brief Scratch storage. 
  */
  EBAMRCellData m_scratch;

  /*!
    @brief Scratch storage
  */
  EBAMRFluxData m_scratchFluxOne;

  /*!
    @brief Scratch storage
  */
  EBAMRFluxData m_scratchFluxTwo;

  /*!
    @brief Scratch storage
  */
  EBAMRIVData m_divG_nc;

  /*!
    @brief Scratch storage
  */
  EBAMRIVData m_mass_diff;

  /*!
    @brief Zero flux
  */
  EBAMRIVData m_eb_zero;

  /*!
    @brief Cached state vector for regrid. DO NOT TOUCH!
  */
  EBAMRCellData m_cache_state;

  /*!
    @brief Cached source term for regrids. DO NOT TOUCH!
  */
  EBAMRCellData m_cache_source;

  /*!
    @brief Face-centered velocities (only normal components)
  */
  EBAMRFluxData m_velo_face;

  /*!
    @brief EB-centered velocities
  */
  EBAMRIVData m_velo_eb;

  /*!
    @brief Flux through the embedded boundary
  */
  EBAMRIVData m_ebflux;

  /*!
    @brief Domain flux
  */
  EBAMRIFData m_domainflux;

  /*!
    @brief Dummy data
  */
  EBAMRCellData m_aco;

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
    @brief Flag for mass-weighted redistribution
  */
  bool m_redist_mass_weighted;

  /*!
    @brief Flag for blending the hybrid divergence
  */
  bool m_blend_conservation;

  /*!
    @brief Also solve for diffusion
  */
  bool m_diffusive;

  /*!
    @brief Also solve for advection/convection
  */
  bool m_mobile;
  
  /*!
    @brief Turn on/off source terms when time-extrapolating
  */
  bool m_extrap_source;

  /*!
    @brief Output state
  */
  bool m_plot_phi;

  /*!
    @brief Output velocities
  */
  bool m_plot_vel;

  /*!
    @brief Output diffusion coefficients
  */
  bool m_plot_dco;

  /*!
    @brief Output EB fluxes
  */
  bool m_plot_ebf;

  /*!
    @brief Output source term
  */
  bool m_plot_src;

  /*!
    @brief Plot numbers or densities
  */
  bool m_plot_numbers;

  /*!
    @brief RNG engine
  */
  std::mt19937_64 m_rng;

  /*!
    @brief RNG seed
  */
  int m_seed;

  /*!
    @brief Average cell-centered velocities to face centers
  */
  virtual void average_velo_to_faces(EBAMRFluxData& a_velo_face, const EBAMRCellData& a_velo_cell);

  /*!
    @brief Set up face-centered advection flux
  */
  virtual void compute_flux(EBAMRFluxData&       a_flux,
			    const EBAMRFluxData& a_face_state,
			    const EBAMRFluxData& a_face_vel,
			    const EBAMRIFData&   a_domain_flux);

  /*!
    @brief Set up face-centered advection flux
  */
  virtual void compute_flux(LevelData<EBFluxFAB>&              a_flux,
			    const LevelData<EBFluxFAB>&        a_face_state,
			    const LevelData<EBFluxFAB>&        a_face_vel,
			    const LevelData<DomainFluxIFFAB>&  a_domain_flux,
			    const int                          a_lvl);

  /*!
    @brief Compute the face-centered diffusion flux
  */
  virtual void compute_diffusion_flux(EBAMRFluxData& a_flux, const EBAMRCellData& a_state);

  /*!
    @brief Compute the face-centered difufsion flux
  */
  virtual void compute_diffusion_flux(LevelData<EBFluxFAB>& a_flux, const LevelData<EBCellFAB>& a_state, const int a_lvl);

  /*!
    @brief Compute conservative divergence from fluxes. 
  */
  virtual void conservative_divergence(EBAMRCellData& a_cons_div, EBAMRFluxData& a_flux, const EBAMRIVData& a_ebflux);

  /*!
    @brief Compute the conservative divergence with the assumption that only the eb-flux is non-zero
  */
  virtual void conservative_divergence_eb(EBAMRCellData& a_consdiv, const EBAMRIVData& a_ebflux);


  /*!
    @brief Compute the non-conservative divergence
  */
  virtual void nonconservative_divergence(EBAMRIVData& a_div_nc, const EBAMRCellData& a_divG);

  /*!
    @brief Use the non-conservative divergence to make the conservative divergence hold the hybrid divergence
  */
  virtual void hybrid_divergence(EBAMRCellData&       a_divF,
				 EBAMRIVData&         a_mass_diff,
				 const EBAMRIVData&   a_NonConservativeDivergenceStencil);

  
  /*!
    @brief Conservative divergence over regular cells
  */
  virtual void consdiv_regular(LevelData<EBCellFAB>& a_divJ, const LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Make the hybrid divergence. On the way in, a_divF_H must hold the conservative divergence
  */
  virtual void hybrid_divergence(LevelData<EBCellFAB>&              a_divF_H,
				 LevelData<BaseIVFAB<Real> >&       a_mass_diff,
				 const LevelData<BaseIVFAB<Real> >& a_divF_nc,
				 const int                          a_lvl);

  /*!
    @brief Compute flux interpolant
  */
  virtual void setup_flux_interpolant(const LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Interpolate flux to centroids
  */
  virtual void interpolate_flux_to_centroids(LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Compute div(F) on irregular cells
  */
  virtual void compute_divG_irreg(LevelData<EBCellFAB>&              a_divG,
				  const LevelData<BaseIVFAB<Real> >& a_ebflux,
				  const int                          a_lvl);

  /*!
    @brief Increment flux registers. 
  */
  virtual void increment_flux_register(const EBAMRFluxData& a_face_state, const EBAMRFluxData& a_velo_face);

  /*!
    @brief Increment flux registers
  */
  virtual void increment_flux_register(const EBAMRFluxData& a_flux);

  /*!
    @brief Reset the flux register
  */
  virtual void reset_flux_register();

  /*!
    @brief Increment the coarse-fine redist with mass difference (only for EB x CF)
  */
  virtual void coarse_fine_increment(const EBAMRIVData& m_mass_diff);

  /*!
    @brief Do mass redistribution
  */
  virtual void hyperbolic_redistribution(EBAMRCellData& a_del_vel_rho, const EBAMRIVData&   a_mass_diff);

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
    @brief Define interpolant
  */
  virtual void define_interpolant();

  /*!
    @brief Parses domain BC options
  */
  virtual void parseDomain_bc();    // Parses domain BC options

  /*!
    @brief Parses plot variables
  */
  virtual void parse_plot_vars();    // Parses plot variables

  /*!
    @brief Parse random number generator
  */
  virtual void parse_rng_seed();

  /*!
    @brief Parse plot mode
  */
  virtual void parse_plotmode();

  /*!
    @brief Parse the conservation
  */
  virtual void parse_conservation();

  /*!
    ®brief Parse source term stuff
  */
  virtual void parse_extrap_source();
};
#include "CD_NamespaceFooter.H"

#endif
