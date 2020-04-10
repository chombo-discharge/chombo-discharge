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
#include "cdr_species.H"

#include <random>
#include <time.h>
#include <chrono>

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
  virtual void parse_options() = 0;

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
  */
  virtual void compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt) = 0;

  /*!
    @brief Compute div(v*phi) explicitly
  */
  virtual void compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt) = 0;

  /*!
    @brief Compute div(D*grad(phi)) explicitly
  */
  virtual void compute_divD(EBAMRCellData& a_diffusive_term, const EBAMRCellData& a_state) = 0;

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
    @brief Make a_phi non-negative by redistributing mass in the vicinity of the cut cells
  */
  virtual void make_non_negative(EBAMRCellData& a_phi);

  /*!
    @brief Average velocities to faces
  */
  virtual void average_velo_to_faces();

  /*!
    @brief Cache internal state
  */
  virtual void cache_state();

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
  virtual void set_species(const RefCountedPtr<cdr_species> a_species);
  
  /*!
    @brief Set computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom);

  /*!
    @brief Set the amr object
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

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
    @brief Use mass-weighted redistribution
  */
  virtual void set_mass_redist(const bool a_mass_redist);

  /*!
    @brief Regrid internals
  */
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Set which variables to output
  */
  virtual void set_plot_variables();

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
    @brief Compute the non-conservative divergence
  */
  virtual void nonconservative_divergence(EBAMRIVData& a_div_nc, const EBAMRCellData& a_divG);

  /*!
    @brief Use the non-conservative divergence to make the conservative divergence hold the hybrid divergence
  */
  virtual void hybrid_divergence(EBAMRCellData&       a_divF,
				 EBAMRIVData&         a_mass_diff,
				 const EBAMRIVData&   a_noncons_div);

  
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
    @brief Mass injection from regular cells into cut cells in order to avoid negative densities. Preserves mass. 
  */
  virtual void concentration_redistribution(EBAMRCellData&       a_phi,
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
    @brief Increment redistribution register to ensure non-negative densities in cut cells
  */
  virtual void increment_concentration_redist(const EBAMRIVData& a_mass_diff);

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
    @brief Define interpolant
  */
  virtual void define_interpolant();

  /*!
    @brief Smoothed face average
  */
  virtual void smooth_heaviside_faces(EBAMRFluxData& a_face_states, const EBAMRCellData& a_cell_states);

  /*!
    @brief Gaussian noise field
  */
  virtual void fill_GWN(EBAMRFluxData& a_noise, const Real a_mean);

  /*!
    @brief Options for parsing some class options
  */
  virtual void parse_mass_redist();  // Parses mass redistribution

  /*!
    @brief Parses domain BC options
  */
  virtual void parse_domain_bc();    // Parses domain BC options

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
};

#endif
