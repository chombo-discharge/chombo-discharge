/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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
#include <BaseEBCellFAB.H>

// Our includes
#include <CD_OpenMP.H>
#include <CD_LevelTiles.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Alias for ParticleData on each AMR level.
*/
template <class P>
using AMRParticles = Vector<RefCountedPtr<ParticleData<P>>>;

/*!
  @brief Alias for cell-sorted particles on each AMR level.
*/
template <class P>
using AMRCellParticles = Vector<RefCountedPtr<LayoutData<BinFab<P>>>>;

/*!
  @brief Templated class for holding particles on an AMR hierarchy with particle remapping. 
  @details The template parameter P should derive from Chombo's BinItem. 
  @note ParticleContainers are, because of Chombo, restricted to square-block grids. This means
  that the blocking factor and maximum box size must be the same. 
*/
template <class P>
class ParticleContainer
{
public:
  using ValidMask = RefCountedPtr<LevelData<BaseFab<bool>>>;

  /*!
    @brief Default constructor. Leaves object in undefined state. 
  */
  ParticleContainer();

  /*!
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
                    const int                                a_blockingFactor,
                    const int                                a_finestLevel,
                    const std::string                        a_realm);

  /*!
    @brief Destructor ( does nothing)
  */
  virtual ~ParticleContainer();

  /*!
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
         const int                                a_blockingFactor,
         const int                                a_finestLevel,
         const std::string                        a_realm);

  /*!
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
         const int                                a_base,
         const int                                a_newFinestLevel);

  /*!
    @brief Cache particles before calling regrid. 
    @param[in] a_base Coarsest grid level which will not change. 
  */
  void
  preRegrid(const int a_base);

  /*!
    @brief Copy particles to mask particle data holder.
    @param[in] a_mask Mask
    @note If the mask is nullptr on any of the levels, those levels will be ignored. 
  */
  void
  copyMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_mask) const;

  /*!
    @brief Copy particles to masked particle data holder. 
    @param[in] a_level Grid level
    @param[in] a_mask  Mask
  */
  void
  copyMaskParticles(const int a_level, const LevelData<BaseFab<bool>>& a_mask) const;

  /*!
    @brief Copy particles to the mask particle data holder.
    @param[in] a_mask Mask
    @note If the mask is nullptr on any of the levels, those levels will be ignored. 
  */
  void
  transferMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_mask);

  /*!
    @brief Copy particles to masked particle data holder. 
    @param[in] a_level Grid level
    @param[in] a_mask  Mask
  */
  void
  transferMaskParticles(const int a_level, const LevelData<BaseFab<bool>>& a_mask);

  /*!
    @brief Sort particles according to the < operator in the particle.
  */
  void
  sortParticles() noexcept;

  /*!
    @brief Clear all particles
  */
  void
  clearParticles();

  /*!
    @brief Clear the buffer particles. 
  */
  void
  clearBufferParticles() const;

  /*!
    @brief Clear the "mask" particles. 
  */
  void
  clearMaskParticles() const;

  /*!
    @brief Clear outcast particles
  */
  void
  clearOutcast() noexcept;

  /*!
    @brief Clear particles on input data holder
    @param[inout] a_particles
  */
  void
  clear(AMRParticles<P>& a_particles) const;

  /*!
    @brief Is cell-sorted or not
  */
  bool
  isOrganizedByCell() const;

  /*!
    @brief Get finest AMR level
  */
  int
  getFinestLevel() const;

  /*!
    @brief Get the realm where this ParticleContainer lives. 
  */
  const std::string
  getRealm() const;

  /*!
    @brief Get all particles on all levels
    @return m_particles
  */
  AMRParticles<P>&
  getParticles();

  /*!
    @brief Get all particles on a level. Const version.
    @return m_particles
  */
  const AMRParticles<P>&
  getParticles() const;

  /*!
    @brief Get buffer particles on all levels
    @return m_bufferParticles
  */
  AMRParticles<P>&
  getBufferParticles();

  /*!
    @brief Get buffer particles on all levels. Const version. 
    @return m_bufferParticles
  */
  const AMRParticles<P>&
  getBufferParticles() const;

  /*!
    @brief Get the mask particles. 
    @return m_maskParticles
  */
  AMRParticles<P>&
  getMaskParticles();

  /*!
    @brief Get the mask particles. Const version. 
    @return m_maskParticles
  */
  const AMRParticles<P>&
  getMaskParticles() const;

  /*!
    @brief Get the AMR grids
  */
  const Vector<DisjointBoxLayout>&
  getGrids() const;

  /*!
    @brief Get lower-left corner
    @return m_probLo
  */
  const RealVect
  getProbLo() const noexcept;

  /*!
    @brief Get grid resolutions
    @return m_dx
  */
  const Vector<RealVect>
  getDx() const noexcept;

  /*!
    @brief Get particle data on a level
    @param[in] a_level AMR level
  */
  ParticleData<P>&
  operator[](const int a_level);

  /*!
    @brief Get data on a level. Const version. 
    @param[in] a_level AMR level
  */
  const ParticleData<P>&
  operator[](const int a_level) const;

  /*!
    @brief Get all cell particles. 
    @note This will issue an error if the particles are not sorted by cell. 
  */
  AMRCellParticles<P>&
  getCellParticles();

  /*!
    @brief Get cell particles, const version. 
    @note This will issue an error if the particles are not sorted by cell. 
  */
  const AMRCellParticles<P>&
  getCellParticles() const;

  /*!
    @brief Get cell particles
    @param[in] a_level Grid level
    @note This will issue an error if the particles are not sorted by cell. 
  */
  LayoutData<BinFab<P>>&
  getCellParticles(const int a_level);

  /*!
    @brief Get cell particles, const version. 
    @param[in] a_level Grid level
    @note This will issue an error if the particles are not sorted by cell. 
  */
  const LayoutData<BinFab<P>>&
  getCellParticles(const int a_level) const;

  /*!
    @brief Get cell particles on a specific level and grid inex. 
    @param[in] a_level Grid level
    @param[in] a_dit   Grid index
    @note This will issue an error if the particles are not sorted by cell. 
  */
  BinFab<P>&
  getCellParticles(const int a_level, const DataIndex a_dit);

  /*!
    @brief Get cell particles, const version. 
    @param[in] a_level Grid level
    @param[in] a_dit   Grid index
    @note This will issue an error if the particles are not sorted by cell. 
  */
  const BinFab<P>&
  getCellParticles(const int a_level, const DataIndex a_dit) const;

  /*!
    @brief Fill a cell-sorted particle data holder with all the particles in the grid patch. 
    @param[out] a_cellParticles Particles in the input grid patch, sorted by cell. 
    @param[in]  a_level Grid level
    @param[in]  a_dit   Grid index
    @note This will define a_cellParticles over the DisjointBoxLayout grid box and fill it with the particles in the patch. 
  */
  void
  getCellParticles(BinFab<P>& a_cellParticles, const int a_lvl, const DataIndex a_dit) const;

  /*!
    @brief Fill a cell-sorted particle data holder with all the particles in the grid patch. The original particles are destroyed. 
    @param[out] a_cellParticles Particles in the input grid patch, sorted by cell. 
    @param[in]  a_level Grid level
    @param[in]  a_dit   Grid index
    @note This will define a_cellParticles over the DisjointBoxLayout grid box and fill it with the particles in the patch. 
  */
  void
  getCellParticlesDestructive(BinFab<P>& a_cellParticles, const int a_lvl, const DataIndex a_dit);

  /*!
    @brief Sort particles by cell
    @details This will fill m_cellSortedParticles and destroy the patch-sorted particles. 
  */
  void
  organizeParticlesByCell();

  /*!
    @brief Sort particles by cell
    @details This will fill m_particles from m_cellSortedParticles and destroy the cell-sorted particles. 
  */
  void
  organizeParticlesByPatch();

  /*!
    @brief Add particles to container
    @param[in] a_particles particles to add to this container. 
  */
  void
  addParticles(const List<P>& a_particles);

  /*!
    @brief Add particles to container. The input particles are destroyed by this routine. 
    @param[in] a_particles particles to add to this container. . 
  */
  void
  addParticlesDestructive(List<P>& a_particles);

  /*!
    @brief Add particles to a certain level and patch. The input particles are cell-sorted. 
    @param[in] a_particles Input particles (cell-sorted)
    @param[in] a_lvl       Grid level
    @param[in] a_dit       Grid index
  */
  void
  addParticles(const BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit);

  /*!
    @brief Add particles to a certain level and patch. The input particles are cell-sorted (and destroyed on output). 
    @param[inout] a_particles Input particles (cell-sorted). Destroyed on exit. 
    @param[in]    a_lvl       Grid level
    @param[in]    a_dit       Grid index
  */
  void
  addParticlesDestructive(BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit);

  /*!
    @brief Add particles from other container to this one
    @param[in] a_otherContainer Other particle container. 
  */
  void
  addParticles(const ParticleContainer<P>& a_otherContainer);

  /*!
    @brief Add particles from other container to this one. This destroys the particles in the other container. 
    @param[inout] a_otherContainer Other particle container (destroyed on output). 
  */
  void
  addParticlesDestructive(ParticleContainer<P>& a_otherContainer);

  /*!
    @brief Remap over the entire AMR hierarchy
  */
  void
  remap();

  /*!
    @brief Do level remaps. This ignores particles moving between levels and moving out of valid regions.
  */
  void
  levelRemap();

  /*!
    @brief Do a level remap. This ignores particles moving between levels and moving out of valid regions. 
    @param[in] a_level Grid level
  */
  void
  levelRemap(const int a_level);

  /*!
    @brief Move particles into this container. 
    @param[inout] a_otherContainer Other particle container. 
    @note  The grids MUST be the same for both containers. That is, the containers MUST be defined over the same realm. 
  */
  void
  transferParticles(ParticleContainer<P>& a_otherContainer);

  /*!
    @brief Move particles into this container. 
    @param[inout] a_otherContainer Other particle container. 
    @note  The grids MUST be the same for both containers. That is, the containers MUST be defined over the same realm. 
  */
  void
  transferParticles(AMRParticles<P>& a_otherContainer);

  /*!
    @brief Evict particles if they move out of the valid region.
    @param[inout] a_evictedParticles List that the particles are moved into.
    @param[inout] a_particles        Particles.
    @param[in]    a_level            Grid level.
  */
  void
  evictInvalidParticles(List<P>& a_evictedParticles, ParticleData<P>& a_particles, const int a_coarseLevel);

  /*!
    @brief Evict particles if they move out of the valid region.
    @param[inout] a_evictedParticles List that the particles are moved into.
    @param[inout] a_particles        Particles.
    @param[in]    a_level            Grid level.
  */
  void
  admitValidParticles(List<P>& a_evictedParticles, ParticleData<P>& a_particles, const int a_coarseLevel);

  /*!
    @brief Set the particle member to the input value
    @details This is a jack-of-all-trades kind of routine for setting a particle field. The template parameter indicates the field to be set, it must be a
    pointer-to-member-function that returns the field to be set. For example, if the particle class P has a field mass which one wants to set, and this field
    is returned by a function Real& P::mass(), one can set the quantity for all particles by calling ParticleContainer<P>::setValue<&P::mass>(Real)
    @param[in] a_value Value 
  */
  template <Real& (P::*particleScalarField)()>
  inline void
  setValue(const Real a_value);

  /*!
    @brief Set the particle member to the input value
    @details This is a jack-of-all-trades kind of routine for setting a particle field. The template parameter indicates the field to be set, it must be a
    pointer-to-member-function that returns the field to be set. For example, if the particle class P has a field vector which one wants to set, and this field
    is returned by a function RealVect& P::velocity(), one can set the quantity for all particles by calling ParticleContainer<P>::setValue<&P::velocity>(RealVect);
    @param[in] a_value Value 
  */
  template <RealVect& (P::*particleVectorField)()>
  inline void
  setValue(const RealVect a_value);

  /*!
    @brief Get local number of particles
  */
  unsigned long long
  getNumberOfValidParticlesLocal() const;

  /*!
    @brief Get global number of particles
  */
  unsigned long long
  getNumberOfValidParticlesGlobal() const;

  /*!
    @brief Get local number of particles
  */
  unsigned long long
  getNumberOfOutcastParticlesLocal() const;

  /*!
    @brief Get global number of particles
  */
  unsigned long long
  getNumberOfOutcastParticlesGlobal() const;

  /*!
    @brief Get the number particles in the halo cells. 
  */
  unsigned long long
  getNumberOfMaskParticlesLocal() const;

  /*!
    @brief Get the number particles in the halo cells. 
  */
  unsigned long long
  getNumberOfMaskParticlesGlobal() const;

protected:
  // Reductions for performing catenation and join operations in a thread-safe manner.
#pragma omp declare reduction(catenate : List <P> : ThreadSafeCatenation <P>(omp_out, omp_in))
#pragma omp declare reduction(join : List <P> : ThreadSafeJoin <P>(omp_out, omp_in))

  /*!
    @brief Realm on which the ParticleContainer is defined
  */
  std::string m_realm;

  /*!
    @brief Blocking factor, aka grid size. 
  */
  int m_blockingFactor;

  /*!
    @brief Finest grid level
  */
  int m_finestLevel;

  /*!
    @brief Lower left corner of the physical domain. 
  */
  RealVect m_probLo;

  /*!
    @brief Check if particle container is defined.
  */
  bool m_isDefined;

  /*!
    @brief Check if particle container is "cell sorted". 
  */
  bool m_isOrganizedByCell;

  /*!
    @brief Profile or not
  */
  bool m_profile;

  /*!
    @brief Debug or not
  */
  bool m_debug;

  /*!
    @brief Verbose or not
  */
  bool m_verbose;

  /*!
    @brief Tiled AMR space
  */
  Vector<RefCountedPtr<LevelTiles>> m_levelTiles;

  /*!
    @brief AMR grids
  */
  Vector<DisjointBoxLayout> m_grids;

  /*!
    @brief Old AMR grids. Used during regrids
  */
  Vector<DisjointBoxLayout> m_oldGrids;

  /*!
    @brief Problem domains
  */
  Vector<ProblemDomain> m_domains;

  /*!
    @brief Map of valid cells; evaluates to false if a cell is covered by a finer grid. 
  */
  Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> m_validRegion;

  /*!
    @brief Grown layouts
  */
  Vector<BoxLayout> m_grownGrids;

  /*!
    @brief Resolutions on each grid level
  */
  Vector<RealVect> m_dx;

  /*!
    @brief Refinement ratios. Entry at index 'l' is the refinement between level 'l' and level 'l+1'
  */
  Vector<int> m_refRat;

  /*!
    @brief The actual particles that this ParticleContainer represents
  */
  mutable AMRParticles<P> m_particles;

  /*!
    @brief Special data holder that holds particles on a grown grid layout. 
    @details This particle data holder uses grids that exceed the DisjointBoxLayout on each level (boxes are grown by m_refRat). This
    means that this particle data holder can store particles that are out the DisjointBoxLayout by m_refRat cells. This is very useful
    when one wants to deposit particles on the fine level, but where the particles live on the coarse level. In such cases one will typically
    copy some (or all) particles from the coarse level to this data holder and deposit from this data holder onto the fine level. 
  */
  mutable AMRParticles<P> m_bufferParticles;

  /*!
    @brief Special data holder that holds copies of particles on the coarse level that are within a specified range of the
    refinement boundary. Used for deposition of coarse-level particle clouds across refinement boundaries. 
  */
  mutable AMRParticles<P> m_maskParticles;

  /*!
    @brief Cache particles. Use during pre-regrid operations, and for transferring particles between coarse and fine before
    deposition steps. 
  */
  Vector<List<P>> m_cacheParticles;

  /*!
    @brief Particles stored on the old grid. 
    @details Populated during the preRegrid call.
  */
  AMRParticles<P> m_regridParticles;

  /*!
    @brief Cell particles. 
  */
  AMRCellParticles<P> m_cellSortedParticles;

  /*!
    @brief Set up grown grids
    @param[in] a_base        Base level
    @param[in] a_finestLevel Finest AMR level
  */
  void
  setupGrownGrids(const int a_base, const int a_finestLevel);

  /*!
    @brief Setup function for the particle data (m_particles and m_maskParticles)
    @param[in] a_base        Base level
    @param[in] a_finestLevel Finest AMR level
  */
  void
  setupParticleData(const int a_base, const int a_finestLevel);

  /*!
    @brief Remap the outcast particles in the particle data
    @param[inout] Particles where outcast particles are gathered.
  */
  void
  gatherOutcast(ParticleData<P>& a_particles) noexcept;

  /*!
    @brief Remap the outcast particles in the particle data
    @param[inout] Particles where outcast particles are gathered.
  */
  void
  remapOutcast(ParticleData<P>& a_particles) noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_ParticleContainerImplem.H>

#endif
