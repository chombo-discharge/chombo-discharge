/*!
  @file   mc_photo.H
  @author Robert Marskar
  @date   Jan. 2018
  @todo   Particle deposition should use the new EBMeshInterp and not MeshInterp
*/

#ifndef _MC_PHOTO_
#define _MC_PHOTO_

#include "rte_solver.H"
#include "EBParticleInterp.H"
#include "particle_container.H"

#include <Particle.H>
#include <ParticleData.H>
#include <ParticleValidRegion.H>
#include <BinItem.H>

#include <random>

/*!
  @brief Radiative tranfer equation solver using stationary Monte-Carlo
*/
class mc_photo : public rte_solver {
public:

  /*!
    @brief Constructor
  */
  mc_photo();

  /*!
    @brief Destructor
  */
  virtual ~mc_photo();

  /*!
    @brief Advance RTE and deposit photons on a_state
  */
  virtual bool advance(const Real a_dt, EBAMRCellData& a_state, const EBAMRCellData& a_source, const bool a_zerophi = false);

  /*!
    Parse options
  */
  virtual void parse_options();
  
  /*!
    @brief Allocate internal storage
  */
  virtual void allocate_internals();

  /*!
    @brief Stuff to do before a regrid
  */
  virtual void pre_regrid(const int a_base, const int a_old_finest_level);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocate_internals();

  /*!
    @brief Regrid function for this class
  */
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Register operators
  */
  virtual void register_operators();

  /*!
    @brief Compute the boundary flux. For Eddington, the boundary flux is = c*phi/2
  */
  virtual void compute_boundary_flux(EBAMRIVData& a_ebflux, const EBAMRCellData& a_state);

  /*!
    @brief Compute the domain flux. For Eddington, the domain flux is = c*phi/2
  */
  virtual void compute_domain_flux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_state);

  /*!
    @brief Compute the flux
    @details For diffusive models, the flux will be something like grad(a_state). For higher-order models, the flux 
    will be contained in a_state (somehow);
  */
  virtual void compute_flux(EBAMRCellData& a_flux, const EBAMRCellData& a_state);

  /*!
    @brief Get isotropic part. 
  */
  virtual void compute_density(EBAMRCellData& a_isotropic, const EBAMRCellData& a_state);

  /*!
    @brief Clear data holder
  */
  virtual void clear();

  /*!
    @brief Clear data holder
  */
  virtual void clear(particle_container<photon>& a_photon);

  /*!
    @brief Clear data holder
  */
  virtual void clear(AMRParticles<photon>& a_photons);

  /*!
    @brief Write plot file
  */
  virtual void write_plot_file();

  /*!
    @brief Write plot data
  */
  virtual void write_plot_data(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write checkpoint data into handle
  */
  virtual void write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void read_checkpoint_level(HDF5Handle& a_handle, const int a_level);  

  /*!
    @brief Return required number of ghost cells
  */
  virtual int query_ghost() const;

  /*!
    @brief Count number of photons in particle list
  */
  virtual int count_photons(const AMRParticles<photon>& a_photons) const;

  /*!
    @brief Count number of outcast photons in particle list
  */
  virtual int count_outcast(const AMRParticles<photon>& a_photons) const;

  /*!
    @brief Get photons
  */
  virtual particle_container<photon>& get_photons();

  /*!
    @brief Get EB photons
  */
  virtual particle_container<photon>& get_eb_photons();

  /*!
    @brief Get domain photons
  */
  virtual particle_container<photon>& get_domain_photons();
  
protected:

  enum photon_generation {
    deterministic = 0,
    stochastic = 1,
  };

  
  enum source_type {
    number   = 0,
    per_vols = 1, 
    per_vol  = 2,
    per_s    = 3,
  };

  /*!
    @brief If true, the NUMBER of of photons will be deposited in each cell
  */
  bool m_deposit_numbers;

  /*!
    @brief Switch for plotting numbers or densities
  */
  bool m_plot_numbers;

  /*!
    @brief Number of photons per pseudophoton
  */
  int  m_max_photons;

  /*!
    @brief RNG seed
  */
  int  m_seed;

  /*!
    @brief Threshold for swapping Poisson and exponential distributions
  */
  int  m_poiss_exp_swap;

  /*!
    @brief PVR bufer
  */
  int  m_pvr_buffer;

  /*!
    @brief Number of photons
  */
  long long m_num_photons;

  /*!
    @brief Photon generation type
  */
  photon_generation m_photogen;

  /*!
    @brief Source type
  */
  source_type m_src_type;

  /*!
    @brief Bisection step length for boundary intersections
  */
  Real m_bisect_step;

  /*!
    @brief RNG engine
  */
  std::mt19937_64* m_rng; // Engine

  /*!
    @brief Distribution on [0,1]
  */
  std::uniform_real_distribution<Real>* m_udist01; 

  /*!
    @brief Distribution on [-1,1]
  */
  std::uniform_real_distribution<Real>* m_udist11; 

  /*!
    @brief Deposition type
  */
  DepositionType::Which m_deposition;

  /*!
    @brief Plot deposition type
  */
  DepositionType::Which m_plot_deposition;

  /*!
    @brief Coarse data for interpolation of deposition clouds
  */
  EBAMRCellData m_scratch;

  /*!
    @brief Scratch storage for holding the non-conservative deposition
  */
  EBAMRIVData m_depositionNC;

  /*!
    @brief Scratch storage for holding the mass difference when using hybrid deposition
  */
  EBAMRIVData m_massDiff;

  /*!
    @brief All particles
  */
  particle_container<photon> m_photons;

  /*!
    @brief This is a particle container which acts as scratch space when you want to add particles through the EB
  */
  particle_container<photon> m_eb_photons;

  /*!
    @brief This is a particle container which acts as scratch space when you want to add particles through the EB
  */
  particle_container<photon> m_domain_photons;

  /*!
    @brief Domain boundary conditions
  */
  Vector<wallbc::which_bc> m_domainbc;

  /*!
    @brief Generate photons
  */
  void generate_photons(particle_container<photon>& a_particles, const EBAMRCellData& a_source, const Real a_dt);

  /*!
    @brief Move photons and absorb them. 
  */
  void move_and_absorb_photons(particle_container<photon>& a_absorbed, particle_container<photon>& a_original, const Real a_dt);

  /*!
    @brief Remap photons
  */
  void remap();

  /*!
    @brief Remap photons
  */
  void remap(particle_container<photon>& a_photons);

  /*!
    @brief Deposit photons
  */
  void deposit_photons();

  /*!
    @brief Deposit photons
  */
  void deposit_photons(EBAMRCellData&                    a_state,
		       const particle_container<photon>& a_particles,
		       const DepositionType::Which&      a_deposition);
  
  /*!
    @brief Deposition photons on mesh
  */
  void deposit_photons(EBAMRCellData& a_state, const AMRParticles<photon>& a_particles, const DepositionType::Which& a_deposition);

  /*!
    @brief Draw photons
  */
  int draw_photons(const Real a_source, const Real a_volume, const Real a_dt);

  /*!
    @brief Mapping function for domain boundary conditions
  */
  int domainbc_map(const int a_dir, const Side::LoHiSide a_side);

  /*!
    @brief Random Poisson trial
  */
  int random_poisson(const Real a_mean);

  /*!
    @brief Random exponential trial
  */
  Real random_exponential(const Real a_mean);

  /*!
    @brief Random direction
  */
  RealVect random_direction();

#if CH_SPACEDIM==2
  /*!
    @brief Random direction in 2D
  */
  RealVect random_direction2D();
#elif CH_SPACEDIM==3
  /*!
    @brief Random direction in 2D
  */
  RealVect random_direction3D();
#endif

    /*!
    @brief This computes the "conservative" deposition, multiplied by kappa
  */
  void deposit_kappaConservative(EBAMRCellData&              a_state,
				 const AMRParticles<photon>& a_particles,
				 const DepositionType::Which a_deposition);

  /*!
    @brief Make the "non-conservative" kappa deposition
  */
  void deposit_nonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC);

  /*!
    @brief Make the hybrid deposition. Also compute the mass difference
  */
  void deposit_hybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_mass_diff, const EBAMRIVData& a_depositionNC);

  /*!
    @brief Increment redistribution objects with mass to be redistributed
  */
  void increment_redist(const EBAMRIVData& a_mass_diff);

  /*!
    @brief Do redistribution on all levels, no coarse-fine stuff
  */
  void level_redistribution(EBAMRCellData& a_state);

  /*!
    @brief Increment the coarse-fine redist with mass difference (only for EB x CF)
  */
  void coarse_fine_increment(const EBAMRIVData& m_mass_diff);

  /*!
    @brief Do coarse-fine redistribution
  */
  void coarse_fine_redistribution(EBAMRCellData& a_state);

  /*!
    @brief Parse RNG options
  */
  void parse_rng();

  /*!
    @brief Parse pseudophotons
  */
  void parse_pseudophotons();

  /*!
    @brief Parse photogeneration type
  */
  void parse_photogen();

  /*!
    @brief Parse source type
  */
  void parse_source_type();

  /*!
    @brief Parse deposition
  */
  void parse_deposition();

  /*!
    @brief Parse bisection step
  */
  void parse_bisect_step();

  /*!
    @brief Parse domain BCs
  */
  void parse_domain_bc();

  /*!
    @brief Parse PVR buffer
  */
  void parse_pvr_buffer();

  /*!
    @brief Parse plot variables
  */
  void parse_plot_vars();
};

#endif
