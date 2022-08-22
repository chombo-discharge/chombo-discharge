/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SuperParticles.H
  @brief  Namespace containing various particle merging utilities.
  @author Robert Marskar
*/

#ifndef CD_SuperParticlePartitioners_H
#define CD_SuperParticlePartitioners_H

#include <CD_SuperParticles.H>
#include <CD_NamespaceHeader.H>

namespace SuperParticles {
  /*!
    @brief Template function for splitting a particle into two new particles
  */
  template <class P>
  using SplitParticle = std::function<std::pair<P, P>(P& a_particles)>;
  
  /*!
    @brief Partitioner that partitions a set of particles into two subsets having approximately the same total mass. 
    @details P is the particle type and P::*weight is a pointer to a member function providing const-reference access
    to the particle weight. Likewise, P::*position should provide const-reference access to the particle position. 
  */
  template <class P,
	    const Real& (P::*weight)() const,
	    const RealVect& (P::*position)() const
	    const SplitParticle<P>>
  ParticlePartitioner<P> PartitionEqualMass =
    [](ParticleList<P>& a_particles) -> std::pair<ParticleList<P>, ParticleList<P>>
    {
     // 1. Figure out along which direction we should split. We split along the "longest axis"
     RealVect loCorner = std::numeric_limits<Real>::max();
     RealVect hiCorner = -std::numeric_limits<Real>::max();
     Real totalWeight = 0.0;
     
     for (const auto& particle : a_particles) {
       const Real& w = (particle.*weight)();
       const RealVect& x = (particle.*position)();
       
       for (int dir = 0; dir < SpaceDim; dir++) {
	 loCorner[dir] = std::min(loCorner[dir], x[dir]);
	 hiCorner[dir] = std::max(hiCorner[dir], x[dir]);	 
       }

       totalWeight += w;
     }

     const int splitDir = (hiCorner-loCorner).maxDir(true);

     // 2. Sort particles along splitDir
     auto sortCrit =
     [d=splitDir](const P& p1, const P&p2) -> bool
     {
       return p1[d] < p2[d];
     };
     std::sort(std::begin(a_particles), std::end(a_particles), sortCrit);

     // 3. Figure out which one is the "median particle". We do this by assuming we have three "sets" of particles; a left
     //    set, a right set, and the median particle itself (which could potentially be split!).

    };
}

#include <CD_NamespaceFooter.H>

#endif
