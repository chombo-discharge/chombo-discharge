/*!
  @file mc_photo.H
  @brief Declaration of a first order Eddington solver
  @author Robert Marskar
  @date Jan. 2018
*/

#ifndef _MC_PHOTO_
#define _MC_PHOTO_

#include "rte_solver.H"
#include "joint_photon.H"

#include <Particle.H>
#include <ParticleData.H>
#include <ParticleValidRegion.H>
#include <BinItem.H>
#include <MeshInterp.H>

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
    @brief Advance RTE onto state a_state
  */
  virtual bool advance(const Real a_dt, EBAMRCellData& a_state, const EBAMRCellData& a_source, const bool a_zerophi = false);
  
  /*!
    @brief Allocate internal storage
  */
  virtual void allocate_internals();
  
  /*!
    @brief Cache state
  */
  virtual void cache_state();

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocate_internals();

  /*!
    @brief Regrid function for this class
  */
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);

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
    @brief Set grid deposition type
  */
  virtual void set_deposition_type();

  /*!
    @brief Set photon generation type
  */
  virtual void set_photon_generation();

  /*!
    @brief Declare what the source contains. Since the physics interface is arbitrary, this might be an integer number, the number per volume, a volumetric rate, or just a rate. In this case some normalization is required. 
  */
  virtual void set_source_type();

  /*!
    @brief Set pseudophotons
  */
  virtual void set_pseudophotons();

  /*!
    @brief Set bisection step legnth
  */
  virtual void set_bisect_step();

  /*!
    @brief Switch for selecting absorption length randomly
    @details Used e.g. if the photon frequency is "random"
  */
  virtual void set_random_kappa();

  /*!
    @brief Clear data holder
  */
  virtual void clear(EBAMRPhotons& a_photons);
  
  /*!
    @brief 
  */
  virtual void set_domain_bc();
  
  /*!
    @brief Write plot file
  */
  virtual void write_plot_file();

  /*!
    @brief Write checkpoint data into handle
  */
  virtual void write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void read_checkpoint_level(HDF5Handle& a_handle, const int a_level);  

  /*!
    @brief Set PVR buffer
  */
  virtual int set_pvr_buffer();

  /*!
    @brief Return required number of ghost cells
  */
  virtual int query_ghost() const;

  /*!
    @brief Count number of photons in particle list
  */
  virtual int count_photons(const EBAMRPhotons& a_photons) const;

  /*!
    @brief Count number of outcast photons in particle list
  */
  virtual int count_outcast(const EBAMRPhotons& a_photons) const;
  
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
    
  bool m_random_kappa;    // Use random kappa or not
  bool m_deposit_numbers; // Deposit numbers of photons in each cell

  int  m_max_photons;     // Number of photons per pseudophoton
  int  m_seed;            // RNG seed
  int  m_poiss_exp_swap;  // Threshold for swapping Poisson and exponential distributions
  int  m_pvr_buffer;      // PVR buffer

  long long m_num_photons;

  photon_generation m_photogen;
  source_type m_src_type;

  Real m_bisect_step;     // Bisection step length
  
  std::mt19937_64* m_rng; // Engine
  std::uniform_real_distribution<Real>* m_udist01; // Uniform real distribuion on [0,1]
  std::uniform_real_distribution<Real>* m_udist11; // Uniform real distribuion on [-1,1]

  InterpType m_deposition; // Mass deposition type

  EBAMRPhotons m_photons;    // Particles
  EBAMRPhotons m_photocache; // Photon cache
  EBAMRPVR m_pvr;            // Valid region
  EBAMRCellData m_crse;      // Photoionization profile from the coarse level photons


  Vector<wallbc::which_bc> m_domainbc;
  void set_rng();
  
  void generate_photons(EBAMRPhotons& a_particles, const EBAMRCellData& a_source, const Real a_dt);

  void move_and_absorb_photons(EBAMRPhotons& a_absorbed, EBAMRPhotons& a_original, const Real a_dt);
  void remap_photons(EBAMRPhotons& a_photons);
  void aggregate_photons(const EBAMRPhotons& a_photons);
  void deposit_photons(EBAMRCellData& a_state, const EBAMRPhotons& a_particles);

  // For mapping photons into bins
  void binmapPhotons(std::map<IntVect, joint_photon, CompIntVect>& a_mip,
                     const List<photon>&                           a_photons,
                     const RealVect&                               a_meshSpacing,
                     const RealVect&                               a_origin);

  int draw_photons(const Real a_source, const Real a_volume, const Real a_dt);
  int domainbc_map(const int a_dir, const Side::LoHiSide a_side);
  int random_poisson(const Real& a_mean);
  Real random_exponential(const Real a_mean);
  RealVect random_direction();
#if CH_SPACEDIM==2
  RealVect random_direction2D();
#elif CH_SPACEDIM==3
  RealVect random_direction3D();
#endif
};

#endif
