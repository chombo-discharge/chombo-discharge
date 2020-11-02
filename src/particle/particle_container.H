/*!
  @file particle_container.H
  @brief Declaration of a class for holding particles on an AMR hierarchy
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _PARTICLE_CONTAINER_
#define _PARTICLE_CONTAINER_

#include <Particle.H>
#include <ParticleData.H>
#include <ParticleValidRegion.H>
#include <BinItem.H>
#include <BinFab.H>
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>
#include "BaseEBCellFAB.H"

template <class P>
using AMRParticles     = Vector<RefCountedPtr<ParticleData<P> > >;

template <class P>
using AMRCellParticles = Vector<RefCountedPtr<LayoutData<BinFab<P> > > >;

using AMRPVR           = Vector<RefCountedPtr<ParticleValidRegion> >;
//using EBAMRPPC         = Vector<RefCountedPtr<BaseEBCellFAB<int> > >;



/*!
  @brief Templated class for holding particles on an AMR hierarchy with particle remapping
template <class P>
*/
template <class P>
class particle_container {
public:

  /*!
    @brief No need for weak construction
  */
  particle_container();

  /*!
    @brief Full constructor
    @note It is extremely important that a_grids uses only square blocks for the grids. 
  */
  particle_container(const Vector<DisjointBoxLayout>& a_grids,
		     const Vector<ProblemDomain>&     a_domains,
		     const Vector<Real>&              a_dx,
		     const Vector<int>&               a_ref_rat,
		     const RealVect&                  a_prob_lo,
		     const int                        a_blocking_factor,
		     const int                        a_finest_level,
		     const int                        a_pvr_buffer,
		     const int                        a_halo,
		     const std::string                a_realm);

  /*!
    @brief Destructor
  */
  ~particle_container();

  /*!
    @brief Define the container. This will do a clear-out of all particles. 
  */
  void define(const Vector<DisjointBoxLayout>& a_grids,
	      const Vector<ProblemDomain>&     a_domains,
	      const Vector<Real>&              a_dx,
	      const Vector<int>&               a_ref_rat,
	      const RealVect&                  a_prob_lo,
	      const int                        a_blocking_factor,
	      const int                        a_finest_level,
	      const int                        a_pvr_buffer,
	      const int                        a_halo_buffer,
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
  void pre_regrid(const int a_base);

  /*!
    @brief Copy particles to halo
  */
  void copy_halo_particles() const;

  /*!
    @brief Copy non-halo particles to non-halo storage
  */
  void copy_non_halo_particles() const;

  /*!
    @brief Frees up halo particles
  */
  void clear_halo_particles() const;

    /*!
    @brief Frees up halo particles
  */
  void clear_non_halo_particles() const;

  /*!
    @brief Clear all particles
  */
  void clear_particles();

  /*!
    @brief Get the realm
  */
  const std::string get_realm() const;

  /*!
    @brief Get all particles on all levels
  */
  AMRParticles<P>& get_particles();

  /*!
    @brief Get all particles on a level
  */
  const AMRParticles<P>& get_particles() const;

  /*!
    @brief Get the "ghosted" coarse particles. I.e. particles that lie on the coarse side of the CF interface. 
  */
  AMRParticles<P>& get_halo_particles();

  /*!
    @brief Get the "ghosted" coarse particles. I.e. particles that lie on the coarse side of the CF interface. 
  */
  const AMRParticles<P>& get_halo_particles() const;

    /*!
    @brief Get the "ghosted" coarse particles. I.e. particles that lie on the coarse side of the CF interface. 
  */
  AMRParticles<P>& get_non_halo_particles();

  /*!
    @brief Get the "ghosted" coarse particles. I.e. particles that lie on the coarse side of the CF interface. 
  */
  const AMRParticles<P>& get_non_halo_particles() const;

  /*!
    @brief Get PVR
  */
  const AMRPVR& get_pvr() const;

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
  AMRCellParticles<P>& get_cell_particles();

  /*!
    @brief Get cell particles, const version. 
  */
  const AMRCellParticles<P>& get_cell_particles() const;
  
  /*!
    @brief Get cell particles
  */
  LayoutData<BinFab<P> >& get_cell_particles(const int a_level);

  /*!
    @brief Get cell particles, const version. 
  */
  const LayoutData<BinFab<P> >& get_cell_particles(const int a_level) const;

  /*!
    @brief Get cell particles
  */
  BinFab<P>& get_cell_particles(const int a_level, const DataIndex a_dit);

  /*!
    @brief Get cell particles, const version. 
  */
  const BinFab<P>& get_cell_particles(const int a_level, const DataIndex a_dit) const;

  /*!
    @brief Get particles and sort them by cell
  */
  void get_cell_particles(BinFab<P>& a_cellParticles, const int a_lvl, const DataIndex a_dit) const;

  /*!
    @brief Get particles and sort them by cell
  */
  void get_cell_particles_destructive(BinFab<P>& a_cellParticles, const int a_lvl, const DataIndex a_dit);

  /*!
    @brief Sort particles by cell
  */
  void sort_particles_by_cell();

  /*!
    @brief Sort particles by cell
  */
  void sort_particles_by_patch();

  /*!
    @brief Add particles to container
  */
  void add_particles(const List<P>& a_particles);

  /*!
    @brief Add particles to container destructively
  */
  void add_particles_destructive(List<P>& a_particles);

  /*!
    @brief Add particles to a certain level and patch
  */
  void add_particles(const BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit);

  /*!
    @brief Add particles to a certain level and patch
  */
  void add_particles_destructive(BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit);

  /*!
    @brief Add particles from other container to this one
  */
  void add_particles(const particle_container<P>& a_otherContainer);

  /*!
    @brief Add particles from other container to this one. This destroys the particles in the other container. 
  */
  void add_particles_destructive(particle_container<P>& a_otherContainer);

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
  void cache_particles(particle_container<P>& a_otherContainer);

  /*!
    @brief Discard particles that are under a PVR level
  */
  void discard_invalid_particles();

  /*!
    @brief Get local number of particles
  */
  size_t get_num_valid_local() const;

  /*!
    @brief Get global number of particles
  */
  size_t get_num_valid_global() const;

  /*!
    @brief Get local number of particles
  */
  size_t get_num_outcast_local() const;

  /*!
    @brief Get global number of particles
  */
  size_t get_num_outcast_global() const;

  /*!
    @brief Get the number particles in the halo cells. 
  */
  size_t get_num_halo_local() const;

    /*!
    @brief Get the number particles in the halo cells. 
  */
  size_t get_num_halo_global() const;

  /*!
    @brief Get the total weight of all particles. Local version. 
  */
  Real get_weight_valid_local() const;

  /*!
    @brief Get the total weight of all particles. Global version
  */
  Real get_weight_valid_global() const;

protected:

  /*!
    @brief Realm on which the particle_container is defined
  */
  std::string m_realm;
  
  /*!
    @brief Blocking factor, aka grid size. 
  */
  int m_blocking_factor;  // Blocking facto

  /*!
    @brief PVR buffer - to be deprecated...
  */
  int m_buffer;

  /*!
    @brief Finest grid level
  */
  int m_finest_level;

  /*!
    @brief Coarse side buffer. Coarse grid particles that are within this number of cells from the refinement boundary
           can be placed on the fine level. 
  */
  int m_halo_buffer;

  /*!
    @brief Lower left corner of the problem domain (aka the "origin")
  */
  RealVect m_prob_lo;

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
  Vector<List<P> > m_cache_particles;

  /*!
    @brief Mask on the coarse level. Cells that are within m_coarse_buffer steps away from the boundary are "true". The rest
           are false.
  */
  Vector<RefCountedPtr<LevelData<BaseFab<bool> > > > m_halo; 

  /*!
    @brief The actual particles that this particle_container represents
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
  void setup_pvr(const int a_base, const int a_finest_level);

  /*!
    @brief Setup function for the particle data (m_particles and m_halo_particles)
  */
  void setup_particle_data(const int a_base, const int a_finest_level);

  /*!
    @brief Special remapping function in case particles have moved across more than one refinement level. Used
           for getting those particles and putting them in the correct data holder. 
  */
  void remap_lost_particles();

  /*!
    @brief Define function for the coarse-level mask (m_halo_mask)
  */
  void define_halo_masks();

  /*!
    @brief Define function for the coarse-level mask (m_halo_mask). Level version. 
  */
  void define_halo_mask(LevelData<BaseFab<bool> >& a_coarMask,
		      const ProblemDomain&       a_domainCoar,
		      const ProblemDomain&       a_domainFine,
		      const DisjointBoxLayout&   a_gridsCoar,
		      const DisjointBoxLayout&   a_gridsFine,
		      const int                  a_buffer,
		      const int                  a_ref_rat);

  /*!
    @brief Gather and broadcast function when setting up the coar mask.
  */
  void gather_and_broadcast(IntVectSet& a_set);

  /*!
    @brief Copy particles to halo. a_level is the coarse level. 
  */
  void copy_halo_particles(const int a_level) const;

    /*!
    @brief Copy non-particles to non-halo data storage. a_level is the coarse level. 
  */
  void copy_non_halo_particles(const int a_level) const;
};

#include "particle_containerI.H"

#endif
