/*!
  @file mc_photo.H
  @brief Declaration of a first order Eddington solver
  @author Robert Marskar
  @date Jan. 2018
*/

#ifndef _MC_PHOTO_
#define _MC_PHOTO_

#include "rte_solver.H"

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
  virtual void regrid(const int a_old_finest_level, const int a_new_finest_level);

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
    @brief Turn on/off subcycling in time
  */
  virtual void set_subcycle();

  /*!
    @brief Set maximum number of photons per pseudophoton
  */
  virtual void set_pseudophotons();
    
#ifdef CH_USE_HDF5
  /*!
    @brief Write plot file
  */
  virtual void write_plot_file();
#endif

  /*!
    @brief Return required number of ghost cells
  */
  virtual int query_ghost() const;

protected:

  bool m_subcycle;        // Subcycling or not
  int  m_max_photons;     // Number of photons per pseudophoton
  int  m_seed;            // RNG seed
  std::mt19937_64* m_rng; // Engine
  std::uniform_real_distribution<Real>* m_udist01; // Uniform real distribuion on [0,1]
  std::uniform_real_distribution<Real>* m_udist11; // Uniform real distribuion on [-1,1]

  bool draw_random;

  EBAMRPhotons m_photons; // Particles
  EBAMRPVR m_pvr;         // Valid region

  void set_rng();
  
  void generate_photons(EBAMRPhotons& a_particles, const EBAMRCellData& a_source, const Real a_dt);
  void move_and_absorb_photons(EBAMRPhotons& a_absorbed, EBAMRPhotons& a_original, const Real a_dt);
  void deposit_photons(EBAMRCellData& a_state, const EBAMRPhotons& a_particles);

  RealVect random_direction();

#if CH_SPACEDIM==2
  RealVect random_direction2D();
#elif CH_SPACEDIM==3
  RealVect random_direction3D();
#endif
};

#endif
