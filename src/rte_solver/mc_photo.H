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

  /*!
    @brief Needs setup
  */
  bool m_needs_setup;

  EBAMRParticles m_particles;
  EBAMRPVR m_pvr;

  void generate_photons(EBAMRParticles& a_particles, const EBAMRCellData& a_source, const Real a_dt);
  void move_photons(EBAMRParticles& a_particles, const Real a_dt);
  void deposit_photons(EBAMRCellData& a_state, const EBAMRParticles& a_particles);
};

#endif
