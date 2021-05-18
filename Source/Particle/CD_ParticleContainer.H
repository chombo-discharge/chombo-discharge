/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ParticleContainer.H
  @brief  Declaration of a class for holding particles on an AMR hierarchy
  @author Robert Marskar
*/

#ifndef CD_ParticleContainer_H
#define CD_ParticleContainer_H

// Chombo includes
#include <Particle.H>
#include <ParticleData.H>
#include <ParticleValidRegion.H>
#include <BinItem.H>
#include <BinFab.H>
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>
#include "BaseEBCellFAB.H"

// Our includes
#include <CD_NamespaceHeader.H>
  
template <class P>
using AMRParticles     = Vector<RefCountedPtr<ParticleData<P> > >;

template <class P>
using AMRCellParticles = Vector<RefCountedPtr<LayoutData<BinFab<P> > > >;

using AMRPVR           = Vector<RefCountedPtr<ParticleValidRegion> >;


/*!
  @brief Templated class for holding particles on an AMR hierarchy with particle remapping
  template <class P>
*/
template <class P>
class ParticleContainer {
public:

  /*!
    @brief No need for weak construction
  */
  ParticleContainer();

  /*!
    @brief Full constructor
    @note It is extremely important that a_grids uses only square blocks for the grids. 
  */
  ParticleContainer(const Vector<DisjointBoxLayout>& a_grids,
		    const Vector<ProblemDomain>&     a_domains,
		    const Vector<Real>&              a_dx,
		    const Vector<int>&               a_ref_rat,
		    const RealVect&                  a_probLo,
		    const int                        a_blocking_factor,
		    const int                        a_finestLevel,
		    const int                        a_pvr_buffer,
		    const std::string                a_realm);

  /*!
    @brief Destructor
  */
  ~ParticleContainer();

  /*!
    @brief Define the container. This will do a clear-out of all particles. 
  */
  void define(const Vector<DisjointBoxLayout>& a_grids,
	      const Vector<ProblemDomain>&     a_domains,
	      const Vector<Real>&              a_dx,
	      const Vector<int>&               a_ref_rat,
	      const RealVect&                  a_probLo,
	      const int                        a_blocking_factor,
	      const int                        a_finestLevel,
	      const int                        a_pvr_buffer,
	      const std::string                a_realm);

  /*!
    @brief Regrid function. a_base is the coarsest grid level which did NOT change
  */
  void regrid(const Vector<DisjointBoxLayout>& a_grids,
	      const Vector<ProblemDomain>&     a_domains,
	      const Vector<Real>&              a_dx,
	      const Vector<int>&               a_ref_rat,
	      const int                        a_base,
	      const int                        a_newFinestLevel);

  /*!
    @brief Cache particles in a pre-regrid form. a_base is the coarsest grid level which will NOT change
  */
  void preRegrid(const int a_base);

  /*!
    @brief Copy particles to halo
  */
  void copyMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool> > > >& a_halo) const;

  /*!
    @brief Copy non-halo particles to non-halo storage
  */
  void copyNonMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool> > > >& a_halo) const;

  /*!
    @brief Copy particles to halo. a_level is the coarse level. 
  */
  void copyMaskParticles(const int a_level, const LevelData<BaseFab<bool> >& a_halo) const;

  /*!
    @brief Copy non-particles to non-halo data storage. a_level is the coarse level. 
  */
  void copyNonMaskParticles(const int a_level, const LevelData<BaseFab<bool> >& a_halo) const;

  /*!
    @brief Frees up halo particles
  */
  void clearMaskParticles() const;

  /*!
    @brief Frees up halo particles
  */
  void clearNonMaskParticles() const;

  /*!
    @brief Clear all particles
  */
  void clearParticles();

  /*!
    @brief Get finest level
  */
  int getFinestLevel() const;

  /*!
    @brief Get the Realm
  */
  const std::string getRealm() const;

  /*!
    @brief Get all particles on all levels
  */
  AMRParticles<P>& getParticles();

  /*!
    @brief Get all particles on a level
  */
  const AMRParticles<P>& getParticles() const;

  /*!
    @brief Get the "ghosted" coarse particles. I.e. particles that lie on the coarse side of the CF interface. 
  */
  AMRParticles<P>& getMaskParticles();

  /*!
    @brief Get the "ghosted" coarse particles. I.e. particles that lie on the coarse side of the CF interface. 
  */
  const AMRParticles<P>& getMaskParticles() const;

  /*!
    @brief Get the "ghosted" coarse particles. I.e. particles that lie on the coarse side of the CF interface. 
  */
  AMRParticles<P>& getNonMaskParticles();

  /*!
    @brief Get the "ghosted" coarse particles. I.e. particles that lie on the coarse side of the CF interface. 
  */
  const AMRParticles<P>& getNonMaskParticles() const;

  /*!
    @brief Get cache particles
  */
  Vector<List<P> >& getCacheParticles() ;

  /*!
    @brief Get cache particles
  */
  const Vector<List<P> >& getCacheParticles() const;

  /*!
    @brief Get PVR
  */
  const AMRPVR& getPVR() const;

  /*!
    @brief Get grids
  */
  const Vector<DisjointBoxLayout>& getGrids() const;

  /*!
    @brief Get data on a level
  */
  ParticleData<P>& operator[](const int a_level);

  /*!
    @brief Get data on a level
  */
  const ParticleData<P>& operator[](const int a_level) const;

  /*!
    @brief Get cell particles
  */
  AMRCellParticles<P>& getCellParticles();

  /*!
    @brief Get cell particles, const version. 
  */
  const AMRCellParticles<P>& getCellParticles() const;
  
  /*!
    @brief Get cell particles
  */
  LayoutData<BinFab<P> >& getCellParticles(const int a_level);

  /*!
    @brief Get cell particles, const version. 
  */
  const LayoutData<BinFab<P> >& getCellParticles(const int a_level) const;

  /*!
    @brief Get cell particles
  */
  BinFab<P>& getCellParticles(const int a_level, const DataIndex a_dit);

  /*!
    @brief Get cell particles, const version. 
  */
  const BinFab<P>& getCellParticles(const int a_level, const DataIndex a_dit) const;

  /*!
    @brief Get particles and sort them by cell
  */
  void getCellParticles(BinFab<P>& a_cellParticles, const int a_lvl, const DataIndex a_dit) const;

  /*!
    @brief Get particles and sort them by cell
  */
  void getCellParticlesDestructive(BinFab<P>& a_cellParticles, const int a_lvl, const DataIndex a_dit);

  /*!
    @brief Sort particles by cell
  */
  void sortParticlesByCell();

  /*!
    @brief Sort particles by cell
  */
  void sortParticlesByPatch();

  /*!
    @brief Add particles to container
  */
  void addParticles(const List<P>& a_particles);

  /*!
    @brief Add particles to container destructively
  */
  void addParticlesDestructive(List<P>& a_particles);

  /*!
    @brief Add particles to a certain level and patch
  */
  void addParticles(const BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit);

  /*!
    @brief Add particles to a certain level and patch
  */
  void addParticlesDestructive(BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit);

  /*!
    @brief Add particles from other container to this one
  */
  void addParticles(const ParticleContainer<P>& a_otherContainer);

  /*!
    @brief Add particles from other container to this one. This destroys the particles in the other container. 
  */
  void addParticlesDestructive(ParticleContainer<P>& a_otherContainer);

  /*!
    @brief Remap over the entire AMR hierarchy
  */
  void remap();

  /*!
    @brief Only do level remaps
  */
  void levelRemap();

  /*!
    @brief Do a level remap
  */
  void levelRemap(const int a_level);

  /*!
    @brief Move the particles to a different container. 
    @note  The grids MUST be the same for both containers
  */
  void cacheParticles(ParticleContainer<P>& a_otherContainer);

  /*!
    @brief Discard particles that are under a PVR level
  */
  void discardInvalidParticles();

  /*!
    @brief Get local number of particles
  */
  size_t getNumberOfValidParticesLocal() const;

  /*!
    @brief Get global number of particles
  */
  size_t getNumberOfValidParticesGlobal() const;

  /*!
    @brief Get local number of particles
  */
  size_t getNumberOfOutcastParticesLocal() const;

  /*!
    @brief Get global number of particles
  */
  size_t getNumberOfOutcastParticesGlobal() const;

  /*!
    @brief Get the number particles in the halo cells. 
  */
  size_t getNumberOfMaskParticesLocal() const;

  /*!
    @brief Get the number particles in the halo cells. 
  */
  size_t getNumberOfMaskParticesGlobal() const;

  /*!
    @brief Get the total weight of all particles. Local version. 
  */
  Real getWeightValidLocal() const;

  /*!
    @brief Get the total weight of all particles. Global version
  */
  Real getWeightValidGlobal() const;

protected:

  /*!
    @brief Realm on which the ParticleContainer is defined
  */
  std::string m_realm;
  
  /*!
    @brief Blocking factor, aka grid size. 
  */
  int m_blockingFactor;  // Blocking facto

  /*!
    @brief PVR buffer - to be deprecated...
  */
  int m_buffer;

  /*!
    @brief Finest grid level
  */
  int m_finestLevel;

  /*!
    @brief Lower left corner of the problem domain (aka the "origin")
  */
  RealVect m_probLo;

  /*!
    @brief Check if particle container is defined
  */
  bool m_is_defined;

  /*!
    @brief Check if particle container is "cell sorted"
  */
  bool m_cell_sorted;

  /*!
    @brief AMR grids
  */
  Vector<DisjointBoxLayout> m_grids;

  /*!
    @brief Problem domains
  */
  Vector<ProblemDomain> m_domains;

  /*!
    @brief Resolutions on each grid level
  */
  Vector<RealVect> m_dx;

  /*!
    @brief Refinement ratios. Entry at index 'l' is the refinement between level 'l' and level 'l+1'
  */
  Vector<int> m_ref_rat;

  /*!
    @brief Cache particles. Use during pre-regrid operations, and for transferring particles between coarse and fine before
    deposition steps. 
  */
  Vector<List<P> > m_cacheParticles;

  /*!
    @brief The actual particles that this ParticleContainer represents
  */
  AMRParticles<P> m_particles;

  /*!
    @brief Special data holder that holds copies of particles on the coarse level that are within a specified range of the
    refinement boundary. Used for deposition of coarse-level particle clouds across refinement boundaries. 
  */
  mutable AMRParticles<P> m_halo_particles;

  /*!
    @brief Special data holder that holds copies of particles on the coarse level that are within a specified range of the
    refinement boundary. Used for deposition of coarse-level particle clouds across refinement boundaries. 
  */
  mutable AMRParticles<P> m_non_halo_particles;

  /*!
    @brief Valid region for each AMR level
  */
  AMRPVR m_pvr;

  /*!
    @brief Cell particles. 
  */
  AMRCellParticles<P> m_cell_particles;

  /*!
    @brief Setup function for the particle valid region
  */
  void setupPVR(const int a_base, const int a_finestLevel);

  /*!
    @brief Setup function for the particle data (m_particles and m_halo_particles)
  */
  void setupParticleData(const int a_base, const int a_finestLevel);

  /*!
    @brief Special remapping function in case particles have moved across more than one refinement level. Used
    for getting those particles and putting them in the correct data holder. 
  */
  void remapLostParticles();
};

#include <CD_NamespaceFooter.H>

#include <CD_ParticleContainerImplem.H>

#endif
