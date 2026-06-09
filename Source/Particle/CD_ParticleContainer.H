/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleContainer.H
  @brief  Declaration of a class for holding particles on an AMR hierarchy
  @author Robert Marskar
*/

#ifndef CD_PARTICLECONTAINER_H
#define CD_PARTICLECONTAINER_H

// Chombo includes
#include <Particle.H>
#include <ParticleData.H>
#include <ParticleValidRegion.H>
#include <BinItem.H>
#include <BinFab.H>
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>
#include <BaseEBCellFAB.H>

// Our includes
#include <CD_OpenMP.H>
#include <CD_LevelTiles.H>
#include <CD_ParticleMap.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Alias for ParticleData on each AMR level.
*/
template <class P>
using AMRParticles = Vector<RefCountedPtr<ParticleData<P>>>;

/**
  @brief Alias for cell-sorted particles on each AMR level.
*/
template <class P>
using AMRCellParticles = Vector<RefCountedPtr<LayoutData<BinFab<P>>>>;

/**
  @brief Templated class for holding particles on an AMR hierarchy with particle remapping. 
  @details The template parameter P should derive from Chombo's BinItem. 
  @note ParticleContainers are, because of Chombo, restricted to square-block grids. This means
  that the blocking factor and maximum box size must be the same. 
*/
template <class P>
class ParticleContainer
{
public:
  /** @brief alias type alias. */
  using ValidMask = RefCountedPtr<LevelData<BaseFab<bool>>>;

  /**
    @brief Default constructor. Leaves object in undefined state. 
  */
  ParticleContainer();

  /**
    @brief Full constructor
    @param[in] a_grids          AMR grids
    @param[in] a_domains        AMR domains
    @param[in] a_dx             Grid resolutions
    @param[in] a_refRat         Refinement ratios
    @param[in] a_validMask      Valid cells (false for cells covered by a finer grid)
    @param[in] a_levelTiles     Tiled AMR grids
    @param[in] a_probLo         Lower-left corner of the physical domain
    @param[in] a_blockingFactor Grid blocking factor
    @param[in] a_finestLevel    Finest grid level
    @param[in] a_realm          Realm where these particles live. 
  */
  ParticleContainer(const Vector<DisjointBoxLayout>&         a_grids,
                    const Vector<ProblemDomain>&             a_domains,
                    const Vector<Real>&                      a_dx,
                    const Vector<int>&                       a_refRat,
                    const Vector<ValidMask>&                 a_validMask,
                    const Vector<RefCountedPtr<LevelTiles>>& a_levelTiles,
                    const RealVect&                          a_probLo,
                    int                                      a_blockingFactor,
                    int                                      a_finestLevel,
                    std::string                              a_realm);

  /**
    @brief Destructor ( does nothing)
  */
  virtual ~ParticleContainer();

  ParticleContainer(const ParticleContainer&) = delete;
  ParticleContainer&
  operator=(const ParticleContainer&)    = delete;
  ParticleContainer(ParticleContainer&&) = delete;
  ParticleContainer&
  operator=(ParticleContainer&&) = delete;

  /**
    @brief Define the container. This will do a clear-out of all particles. 
    @param[in] a_grids          AMR grids
    @param[in] a_domains        AMR domains
    @param[in] a_dx             Grid resolutions
    @param[in] a_refRat         Refinement ratios
    @param[in] a_validMask      Valid cells (false for cells covered by a finer grid)
    @param[in] a_levelTiles     Tiled AMR grids
    @param[in] a_probLo         Lower-left corner of the physical domain
    @param[in] a_blockingFactor Grid blocking factor
    @param[in] a_finestLevel    Finest grid level
    @param[in] a_realm          Realm where these particles live. 
  */
  void
  define(const Vector<DisjointBoxLayout>&         a_grids,
         const Vector<ProblemDomain>&             a_domains,
         const Vector<Real>&                      a_dx,
         const Vector<int>&                       a_refRat,
         const Vector<ValidMask>&                 a_validMask,
         const Vector<RefCountedPtr<LevelTiles>>& a_levelTiles,
         const RealVect&                          a_probLo,
         int                                      a_blockingFactor,
         int                                      a_finestLevel,
         std::string                              a_realm);

  /**
    @brief Regrid function. a_base is the coarsest grid level which did not change
    @param[in] a_grids          AMR grids
    @param[in] a_domains        AMR domains
    @param[in] a_dx             Grid resolutions
    @param[in] a_refRat         Refinement ratios
    @param[in] a_validMask      Valid cells
    @param[in] a_levelTiles     Tiled AMR grids
    @param[in] a_base           Coarsest grid level that did not change. 
    @param[in] a_newFinestLevel New finest grid level
  */
  void
  regrid(const Vector<DisjointBoxLayout>&         a_grids,
         const Vector<ProblemDomain>&             a_domains,
         const Vector<Real>&                      a_dx,
         const Vector<int>&                       a_refRat,
         const Vector<ValidMask>&                 a_validMask,
         const Vector<RefCountedPtr<LevelTiles>>& a_levelTiles,
         int                                      a_base,
         int                                      a_newFinestLevel);

  /**
    @brief Cache particles before calling regrid. 
    @param[in] a_base Coarsest grid level which will not change. 
  */
  void
  preRegrid(int a_lmin);

  /**
    @brief Compute new particle IDs.
  */
  void
  resetParticleIDs() noexcept;

  /**
    @brief Copy particles to mask particle data holder.
    @param[in] a_mask Mask
    @note If the mask is nullptr on any of the levels, those levels will be ignored. 
  */
  void
  copyMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_mask) const;

  /**
    @brief Copy particles to masked particle data holder. 
    @param[in] a_level Grid level
    @param[in] a_mask  Mask
  */
  void
  copyMaskParticles(int a_level, const LevelData<BaseFab<bool>>& a_mask) const;

  /**
    @brief Copy particles to the mask particle data holder.
    @param[in] a_mask Mask
    @note If the mask is nullptr on any of the levels, those levels will be ignored. 
  */
  void
  transferMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_mask);

  /**
    @brief Copy particles to masked particle data holder. 
    @param[in] a_level Grid level
    @param[in] a_mask  Mask
  */
  void
  transferMaskParticles(int a_level, const LevelData<BaseFab<bool>>& a_mask);

  /**
    @brief Sort particles according to the < operator in the particle.
  */
  void
  sortParticles() noexcept;

  /**
    @brief Clear all particles
  */
  void
  clearParticles();

  /**
    @brief Clear the buffer particles. 
  */
  void
  clearBufferParticles() const;

  /**
    @brief Clear the "mask" particles. 
  */
  void
  clearMaskParticles() const;

  /**
    @brief Clear outcast particles
  */
  void
  clearOutcast() noexcept;

  /**
    @brief Clear particles on input data holder
    @param[in,out] a_particles
  */
  void
  clear(AMRParticles<P>& a_particles) const;

  /**
    @brief Is cell-sorted or not
    @return True if organized by cell, false otherwise
  */
  bool
  isOrganizedByCell() const;

  /**
    @brief Get finest AMR level
    @return Finest level
  */
  int
  getFinestLevel() const;

  /**
    @brief Get the realm where this ParticleContainer lives. 
    @return Realm
  */
  const std::string
  getRealm() const;

  /**
    @brief Get all particles on all levels
    @return m_particles
  */
  AMRParticles<P>&
  getParticles();

  /**
    @brief Get all particles on a level. Const version.
    @return m_particles
  */
  const AMRParticles<P>&
  getParticles() const;

  /**
    @brief Get buffer particles on all levels
    @return m_bufferParticles
  */
  AMRParticles<P>&
  getBufferParticles();

  /**
    @brief Get buffer particles on all levels. Const version. 
    @return m_bufferParticles
  */
  const AMRParticles<P>&
  getBufferParticles() const;

  /**
    @brief Get the mask particles. 
    @return m_maskParticles
  */
  AMRParticles<P>&
  getMaskParticles();

  /**
    @brief Get the mask particles. Const version. 
    @return m_maskParticles
  */
  const AMRParticles<P>&
  getMaskParticles() const;

  /**
    @brief Get the AMR grids
    @return Grids
  */
  const Vector<DisjointBoxLayout>&
  getGrids() const;

  /**
    @brief Get lower-left corner
    @return m_probLo
  */
  const RealVect
  getProbLo() const noexcept;

  /**
    @brief Get grid resolutions
    @return m_dx
  */
  const Vector<RealVect>
  getDx() const noexcept;

  /**
    @brief Get particle data on a level
    @param[in] a_level AMR level
    @return Return value
  */
  ParticleData<P>&
  operator[](int a_lvl);

  /**
    @brief Get data on a level. Const version.
    @param[in] a_lvl AMR level
    @return Return value
  */
  const ParticleData<P>&
  operator[](int a_level) const;

  /**
    @brief Get all cell particles. 
    @note This will issue an error if the particles are not sorted by cell. 
    @return Cell particles
  */
  AMRCellParticles<P>&
  getCellParticles();

  /**
    @brief Get cell particles, const version. 
    @note This will issue an error if the particles are not sorted by cell. 
    @return Cell particles
  */
  const AMRCellParticles<P>&
  getCellParticles() const;

  /**
    @brief Get cell particles
    @param[in] a_level Grid level
    @note This will issue an error if the particles are not sorted by cell. 
    @return Cell particles
  */
  LayoutData<BinFab<P>>&
  getCellParticles(int a_level);

  /**
    @brief Get cell particles, const version.
    @param[in] a_level Grid level
    @note This will issue an error if the particles are not sorted by cell.
    @return Cell particles
  */
  const LayoutData<BinFab<P>>&
  getCellParticles(int a_level) const;

  /**
    @brief Get cell particles on a specific level and grid inex.
    @param[in] a_level Grid level
    @param[in] a_dit   Grid index
    @note This will issue an error if the particles are not sorted by cell.
    @return Cell particles
  */
  BinFab<P>&
  getCellParticles(int a_level, DataIndex a_dit);

  /**
    @brief Get cell particles, const version.
    @param[in] a_level Grid level
    @param[in] a_dit   Grid index
    @note This will issue an error if the particles are not sorted by cell.
    @return Cell particles
  */
  const BinFab<P>&
  getCellParticles(int a_level, DataIndex a_dit) const;

  /**
    @brief Fill a cell-sorted particle data holder with all the particles in the grid patch. 
    @param[out] a_cellParticles Particles in the input grid patch, sorted by cell. 
    @param[in]  a_dit   Grid index
    @note This will define a_cellParticles over the DisjointBoxLayout grid box and fill it with the particles in the patch. 
    @param[in] a_lvl Lvl
  */
  void
  getCellParticles(BinFab<P>& a_cellParticles, int a_lvl, DataIndex a_dit) const;

  /**
    @brief Fill a cell-sorted particle data holder with all the particles in the grid patch. The original particles are destroyed. 
    @param[out] a_cellParticles Particles in the input grid patch, sorted by cell. 
    @param[in]  a_dit   Grid index
    @note This will define a_cellParticles over the DisjointBoxLayout grid box and fill it with the particles in the patch. 
    @param[in] a_lvl Lvl
  */
  void
  getCellParticlesDestructive(BinFab<P>& a_cellParticles, int a_lvl, DataIndex a_dit);

  /**
    @brief Sort particles by cell
    @details This will fill m_cellSortedParticles and destroy the patch-sorted particles. 
  */
  void
  organizeParticlesByCell();

  /**
    @brief Sort particles by cell
    @details This will fill m_particles from m_cellSortedParticles and destroy the cell-sorted particles. 
  */
  void
  organizeParticlesByPatch();

  /**
    @brief Add particles to container
    @param[in] a_particles particles to add to this container. 
  */
  void
  addParticles(const List<P>& a_particles);

  /**
    @brief Add particles to container. The input particles are destroyed by this routine. 
    @param[in] a_particles particles to add to this container. . 
  */
  void
  addParticlesDestructive(List<P>& a_particles);

  /**
    @brief Add particles to a certain level and patch. The input particles are cell-sorted. 
    @param[in] a_particles Input particles (cell-sorted)
    @param[in] a_lvl       Grid level
    @param[in] a_dit       Grid index
  */
  void
  addParticles(const BinFab<P>& a_particles, int a_lvl, DataIndex a_dit);

  /**
    @brief Add particles to a certain level and patch. The input particles are cell-sorted (and destroyed on output). 
    @param[in,out] a_particles Input particles (cell-sorted). Destroyed on exit. 
    @param[in]    a_lvl       Grid level
    @param[in]    a_dit       Grid index
  */
  void
  addParticlesDestructive(BinFab<P>& a_particles, int a_lvl, DataIndex a_dit);

  /**
    @brief Add particles from other container to this one
    @param[in] a_otherContainer Other particle container. 
  */
  void
  addParticles(const ParticleContainer<P>& a_otherContainer);

  /**
    @brief Add particles from other container to this one. This destroys the particles in the other container. 
    @param[in,out] a_otherContainer Other particle container (destroyed on output). 
  */
  void
  addParticlesDestructive(ParticleContainer<P>& a_otherContainer);

  /**
    @brief Run a sanity check and make sure all particles are in their correctly assigned box
  */
  void
  sanityCheck() const noexcept;

  /**
    @brief Remap over the entire AMR hierarchy
  */
  void
  remap();

  /**
    @brief Move particles into this container. 
    @param[in,out] a_otherContainer Other particle container. 
    @note  The grids MUST be the same for both containers. That is, the containers MUST be defined over the same realm. 
  */
  void
  transferParticles(ParticleContainer<P>& a_otherContainer);

  /**
    @brief Move particles into this container. 
    @param[in,out] a_otherContainer Other particle container. 
    @note  The grids MUST be the same for both containers. That is, the containers MUST be defined over the same realm. 
  */
  void
  transferParticles(AMRParticles<P>& a_otherContainer);

  /**
    @brief Evict particles if they move out of the valid region.
    @param[in,out] a_evictedParticles List that the particles are moved into.
    @param[in,out] a_particles        Particles.
    @param[in] a_coarseLevel Coarse level
  */
  void
  admitValidParticles(List<P>& a_evictedParticles, ParticleData<P>& a_particles, int a_coarseLevel);

  /**
    @brief Set the particle member to the input value
    @details This is a jack-of-all-trades kind of routine for setting a particle field. The template parameter indicates the field to be set, it must be a
    pointer-to-member-function that returns the field to be set. For example, if the particle class P has a field mass which one wants to set, and this field
    is returned by a function Real& P::mass(), one can set the quantity for all particles by calling ParticleContainer<P>::setValue<&P::mass>(Real)
    @param[in] a_value Value 
  */
  template <Real& (P::*particleScalarField)()>
  inline void
  setValue(Real a_value);

  /**
    @brief Set the particle member to the input value
    @details This is a jack-of-all-trades kind of routine for setting a particle field. The template parameter indicates the field to be set, it must be a
    pointer-to-member-function that returns the field to be set. For example, if the particle class P has a field vector which one wants to set, and this field
    is returned by a function RealVect& P::velocity(), one can set the quantity for all particles by calling ParticleContainer<P>::setValue<&P::velocity>(RealVect);
    @param[in] a_value Value 
  */
  template <RealVect& (P::*particleVectorField)()>
  inline void
  setValue(RealVect a_value);

  /**
    @brief Get local number of particles
    @return Number of valid particles local
  */
  unsigned long long
  getNumberOfValidParticlesLocal() const;

  /**
    @brief Get global number of particles
    @return Number of valid particles global
  */
  unsigned long long
  getNumberOfValidParticlesGlobal() const;

  /**
    @brief Get local number of particles
    @return Number of outcast particles local
  */
  unsigned long long
  getNumberOfOutcastParticlesLocal() const;

  /**
    @brief Get global number of particles
    @return Number of outcast particles global
  */
  unsigned long long
  getNumberOfOutcastParticlesGlobal() const;

  /**
    @brief Get the number particles in the halo cells. 
    @return Number of mask particles local
  */
  unsigned long long
  getNumberOfMaskParticlesLocal() const;

  /**
    @brief Get the number particles in the halo cells. 
    @return Number of mask particles global
  */
  unsigned long long
  getNumberOfMaskParticlesGlobal() const;

protected:
  // Reductions for performing catenation and join operations in a thread-safe manner.
#pragma omp declare reduction(catenate : List<P> : ThreadSafeCatenation<P>(omp_out, omp_in))
#pragma omp declare reduction(join : List<P> : ThreadSafeJoin<P>(omp_out, omp_in))

  /**
    @brief Realm on which the ParticleContainer is defined
  */
  std::string m_realm;

  /**
    @brief Blocking factor, aka grid size. 
  */
  int m_blockingFactor;

  /**
    @brief Finest grid level
  */
  int m_finestLevel;

  /**
    @brief Lower left corner of the physical domain. 
  */
  RealVect m_probLo;

  /**
    @brief Check if particle container is defined.
  */
  bool m_isDefined;

  /**
    @brief Check if particle container is "cell sorted". 
  */
  bool m_isOrganizedByCell;

  /**
    @brief Profile or not
  */
  bool m_profile;

  /**
    @brief Debug or not
  */
  bool m_debug;

  /**
    @brief Verbose or not
  */
  bool m_verbose;

  /**
    @brief Tiled AMR space
  */
  Vector<RefCountedPtr<LevelTiles>> m_levelTiles;

  /**
    @brief AMR grids
  */
  Vector<DisjointBoxLayout> m_grids;

  /**
    @brief Old AMR grids. Used during regrids
  */
  Vector<DisjointBoxLayout> m_oldGrids;

  /**
    @brief Problem domains
  */
  Vector<ProblemDomain> m_domains;

  /**
    @brief Map of valid cells; evaluates to false if a cell is covered by a finer grid. 
  */
  Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> m_validRegion;

  /**
    @brief Grown layouts
  */
  Vector<BoxLayout> m_grownGrids;

  /**
    @brief Resolutions on each grid level
  */
  Vector<RealVect> m_dx;

  /**
    @brief Refinement ratios. Entry at index 'l' is the refinement between level 'l' and level 'l+1'
  */
  Vector<int> m_refRat;

  /**
    @brief The actual particles that this ParticleContainer represents
  */
  mutable AMRParticles<P> m_particles;

  /**
    @brief Special data holder that holds particles on a grown grid layout. 
    @details This particle data holder uses grids that exceed the DisjointBoxLayout on each level (boxes are grown by m_refRat). This
    means that this particle data holder can store particles that are out the DisjointBoxLayout by m_refRat cells. This is very useful
    when one wants to deposit particles on the fine level, but where the particles live on the coarse level. In such cases one will typically
    copy some (or all) particles from the coarse level to this data holder and deposit from this data holder onto the fine level. 
  */
  mutable AMRParticles<P> m_bufferParticles;

  /**
    @brief Special data holder that holds copies of particles on the coarse level that are within a specified range of the
    refinement boundary. Used for deposition of coarse-level particle clouds across refinement boundaries. 
  */
  mutable AMRParticles<P> m_maskParticles;

  /**
    @brief Particles stored on the old grid. 
    @details Populated during the preRegrid call.
  */
  AMRParticles<P> m_cacheParticles;

  /**
    @brief Cell particles. 
  */
  AMRCellParticles<P> m_cellSortedParticles;

  /**
    @brief Set up grown grids
    @param[in] a_base        Base level
    @param[in] a_finestLevel Finest AMR level
  */
  void
  setupGrownGrids(int a_base, int a_finestLevel);

  /**
    @brief Setup function for the particle data (m_particles and m_maskParticles)
    @param[in] a_base        Base level
    @param[in] a_finestLevel Finest AMR level
  */
  void
  setupParticleData(int a_base, int a_finestLevel);

  /**
    @brief Gather the particles onto a single list
    @param[in,out] a_list List containing all the particles in a_particles
    @param[in,out] a_particles Particles to transfer to list
    @note This is supposed to get inside an OpenMP parallel region.
  */
  inline void
  transferParticlesToSingleList(List<P>& a_list, AMRParticles<P>& a_particles) const noexcept;

  /**
    @brief Copy the input particles onto a single list
    @param[in,out] a_list List containing all the particles in a_particles
    @param[in,out] a_particles Particles to transfer to list
    @note This is supposed to get inside an OpenMP parallel region.
  */
  inline void
  copyParticlesToSingleList(List<P>& a_list, const AMRParticles<P>& a_particles) const noexcept;

  /**
    @brief Iterate through the unmapped particles and map them to proper level, grid indices, and grid ranks
    @param[in,out] a_mappedParticles Particles mapped according to rank ownership (vector index) and grid level and index (std::pair)
    @param[in,out] a_unmappedParticles Particles that are not mapped (yet)
  */
  inline void
  mapParticlesToAMRGrid(std::vector<ParticleMap<List<P>>>& a_mappedParticles,
                        List<P>&                           a_unmappedParticles) const noexcept;

  /**
    @brief Catenate the particles. This is usually called within OpenMP parallel regions.
    @param[in,out] a_globalParticles Global particle map
    @param[in,out] a_localParticles Local particle map
   */
  inline void
  catenateParticleMaps(std::vector<ParticleMap<List<P>>>& a_globalParticles,
                       std::vector<ParticleMap<List<P>>>& a_localParticles) const noexcept;

  /**
    @brief Gather particles locally
    @param[in,out] a_mappedParticles Particles that have been mapped to this rank. 
    @param[in,out] a_particleData Where to assign the mapped particles
  */
  inline void
  assignLocalParticles(ParticleMap<List<P>>& a_mappedParticles, AMRParticles<P>& a_particleData) const noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_ParticleContainerImplem.H>

#endif
