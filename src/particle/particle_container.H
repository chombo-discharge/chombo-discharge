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
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>

template <class P>
using AMRParticles   = Vector<RefCountedPtr<ParticleData<P> > >;
using EBAMRPVR       = Vector<RefCountedPtr<ParticleValidRegion> >;
using EBAMRParticles = Vector<RefCountedPtr<ParticleData<Particle> > >;

/*!
  @brief Templated class for holding particles on an AMR hierarchy with particle remapping
template <class P>
*/
template <class P>
class particle_container {
public:

  /*!
    @brief Full constructor
  */
  particle_container(const Vector<DisjointBoxLayout>& a_grids,
		     const Vector<ProblemDomain>&     a_domains,
		     const Vector<Real>&              a_dx,
		     const Vector<int>&               a_ref_rat,
		     const RealVect&                  a_origin);

  /*!
    @brief Destructor
  */
  ~particle_container();

  /*!
    @brief Get all particles on all levels
  */
  AMRParticles<P>& get_particles();

  /*!
    @brief Get data on a level
  */
  ParticleData<P>& operator[](const int a_level);

  /*!
    @brief Add particles to container
  */
  void add_particles(const List<P>& a_particles);

  /*!
    @brief Add particles to container destructively
  */
  void add_particles_destructive(List<P>& a_particles);

protected:

  Vector<DisjointBoxLayout> m_grids;
  Vector<ProblemDomain> m_domains;
  Vector<RealVect> m_dx;
  Vector<int> m_ref_rat;
  RealVect m_origin;

  AMRParticles<P> m_particles;
  EBAMRPVR m_pvr;
  
  /*!
    @brief No need for weak construction
  */
  particle_container();
};

#include "particle_containerI.H"

#endif
