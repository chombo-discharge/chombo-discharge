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

template <class P>
using AMRParticles     = Vector<RefCountedPtr<ParticleData<P> > >;

template <class P>
using AMRCellParticles = Vector<RefCountedPtr<LayoutData<BinFab<P> > > >;

using AMRPVR           = Vector<RefCountedPtr<ParticleValidRegion> >;
using EBAMRPPC         = Vector<RefCountedPtr<BaseEBCellFAB<int> > >;

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
		     const int                        a_pvr_buffer);

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
	      const int                        a_pvr_buffer);

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
    @brief Get all particles on all levels
  */
  AMRParticles<P>& get_particles();

  /*!
    @brief Get all particles on a level
  */
  const AMRParticles<P>& get_particles() const;

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
    @brief Add particles from other container to this one.
  */
  void add_particles_destructive(const particle_container<P>& a_otherContainer);

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

protected:

  Vector<DisjointBoxLayout> m_grids;
  Vector<ProblemDomain> m_domains;
  Vector<RealVect> m_dx;
  Vector<int> m_ref_rat;
  RealVect m_prob_lo;
  int m_blocking_factor;
  int m_buffer;
  int m_finest_level;

  bool m_is_defined;
  bool m_cell_sorted;

  AMRParticles<P> m_particles;
  AMRPVR m_pvr;
  AMRCellParticles<P> m_cell_particles;

  void setup_pvr(const int a_base, const int a_finest_level);
  void setup_particle_data(const int a_base, const int a_finest_level);
  void remap_lost_particles();
};

#include "particle_containerI.H"

#endif
