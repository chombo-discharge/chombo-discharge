/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ParticleOps.H
  @brief  Declaration of a static class containing some common useful particle routines that would otherwise be replicated across classes. 
  @author Robert Marskar
  @todo   In ebIntersectionRaycast we should check if the path is tangential to the EB since that can consume a lot of resources. 
*/

#ifndef CD_ParticleOps_H
#define CD_ParticleOps_H

// Chombo includes
#include <RefCountedPtr.H>
#include <BaseIF.H>
#include <RealVect.H>

// Our includes
#include <CD_ParticleContainer.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief A shallow static class for doing various kinds of particle-related operations.
*/
class ParticleOps
{
public:
  /*!
    @brief Compute the intersection point between a particle path and a domain side.
    @details If the particle crossed the domain side, the output argument a_s determines the intersection point by x(s) = x0 + s*(x1-x0) where
    x0 is the starting position (a_oldPos), x1 is the end position (a_newPos), and x(s) is the intersection point with the boundary. 
    @param[in]  a_oldPos Particle starting position
    @param[in]  a_newPos Particle end position
    @param[in]  a_probLo Lower-left corner in physical coordinates
    @param[in]  a_probHi Upper-right corner in physical coordinates
    @param[out] a_s      Relative length along the path 
    @return Returns true if the particle crossed the domain boundary.
  */
  static inline bool
  domainIntersection(const RealVect& a_oldPos,
                     const RealVect& a_newPos,
                     const RealVect& a_probLo,
                     const RealVect& a_probHi,
                     Real&           a_s);

  /*!
    @brief Compute the intersection point between a particle path and an implicit function using a bisection algorithm. 
    @details This routine will bisect the particle path into intervals of length a_bisectStep. Each interval is examined for a crossing with the implicit function
    and if the particle crossed the into the EB, the output argument a_s determines the intersection point by x(s) = x0 + s*(x1-x0) where x0 is the starting
    position (a_oldPos), x1 is the end position (a_newPos), and x(s) is the intersection point with the geometry. 
    @param[in]  a_impFunc    Implicit function. 
    @param[in]  a_oldPos     Particle starting position
    @param[in]  a_newPos     Particle end position
    @param[in]  a_bisectStep Bisection step length
    @param[out] a_s          Relative length along the path 
    @return Returns true if the particle crossed into the EB. 
  */
  static inline bool
  ebIntersectionBisect(const RefCountedPtr<BaseIF>& a_impFunc,
                       const RealVect&              a_oldPos,
                       const RealVect&              a_newPos,
                       const Real&                  a_bisectStep,
                       Real&                        a_s);

  /*!
    @brief Compute the intersection point between a particle path and an implicit function using a ray-casting algorithm.
    @details This routine will assume that a_impcFunc is an approximation to the signed distance function and compute the intersection point using ray-marching. Starting
    on the starting position we compute the distance to the EB and move the particle the compute distance along the particle trajectory. If the particle comes too close
    to the EB (within a_tolerance) we consider it to be absorbed. In that case the output argument a_s determines the intersection point by x(s) = x0 + s*(x1-x0) where 
    x0 is the starting position (a_oldPos), x1 is the end position (a_newPos), and x(s) is the intersection point with the geometry. 
    @param[in]  a_impFunc    Implicit function. 
    @param[in]  a_oldPos     Particle starting position
    @param[in]  a_newPos     Particle end position
    @param[in]  a_tolerace   Tolerance for intersectinon method. If the distance to the EB is less than this, the particle is absorbed. 
    @param[out] a_s          Relative length along the path 
    @return Returns true if the particle crossed into the EB.
  */
  static inline bool
  ebIntersectionRaycast(const RefCountedPtr<BaseIF>& a_impFunc,
                        const RealVect&              a_oldPos,
                        const RealVect&              a_newPos,
                        const Real&                  a_tolerance,
                        Real&                        a_s);

  /*!
    @brief Compute particle weights. 
    @param[out] a_weight               Particle weights for "whole" computational particles. 
    @param[out] a_num                  Number of "whole" computational particles
    @param[out] a_remainder            Weight remainder
    @param[in]  a_numPhysicalParticles Number of physical particles
    @param[in]  a_ppc                  Desired number of computational particles
    @details This takes a specified number of physical particles and "desired" number of computational particles and computes
    the number of "whole" computational particles and a remainder. So, if we have 10 physical particles but want 3 computational particles
    we have a_weight = 3, a_num=3, and a_remainder = 1
  */
  static inline void
  computeParticleWeights(unsigned long long&      a_weight,
                         unsigned long long&      a_num,
                         unsigned long long&      a_remainder,
                         const unsigned long long a_numPhysicalParticles,
                         const int                a_ppc);

  /*!
    @brief Partition particles so that all MPI rank draw a cumulative number of particles equal to a_numParticles
    @param[in] a_numParticles Total number of particles to be drawn. 
    @return Returns the number of particles drawn by each rank. 
  */
  template <typename T>
  static inline T
  partitionParticles(const T a_numParticles);

  /*!
    @brief Draw a specified number of random particles distributed according to a predefined distribution. 
    @details The implementation will partition a_numParticles among the ranks so the total number of particles drawn is a_numParticles.
    P is the particle type (must derived from BinItem) and T is an integer type. 
    @param[out] a_particles    Output particles (different for each rank)
    @param[in]  a_numParticles Number of particles
    @param[in]  a_distribution Particle distribution. Must return a random position when calling the function. 
  */
  template <typename P, typename T>
  static inline void
  drawRandomParticles(List<P>& a_particles, const T a_numParticles, const std::function<RealVect()>& a_distribution);

  /*!
    @brief Draw particles in a sphere
    @details The implementation will partition a_numParticles among the ranks so the total number of particles drawn is a_numParticles.
    P is the particle type (must have a RealVect& position() function) and T is an integer type. 
    @param[out] a_particles    Output particles (different for each rank)
    @param[in]  a_numParticles Number of particles
    @param[in]  a_center       Sphere center
    @param[in]  a_radius       Sphere radius
  */
  template <typename P, typename T>
  static inline void
  drawSphereParticles(List<P>&       a_particles,
                      const T        a_numParticles,
                      const RealVect a_center,
                      const Real     a_radius) noexcept;

  /*!
    @brief Draw particles in a box
    @details The implementation will partition a_numParticles among the ranks so the total number of particles drawn is a_numParticles.
    P is the particle type (must have a RealVect& position() function) and T is an integer type. 
    @param[out] a_particles    Output particles (different for each rank)
    @param[in]  a_numParticles Number of particles
    @param[in]  a_loCorner     Low corner of box
    @param[in]  a_hiCorner     High corner of box
  */
  template <typename P, typename T>
  static inline void
  drawBoxParticles(List<P>&       a_particles,
                   const T        a_numParticles,
                   const RealVect a_loCorner,
                   const RealVect a_hiCorner) noexcept;

  /*!
    @brief Draw particles from a Gaussian distribution
    @details The implementation will partition a_numParticles among the ranks so the total number of particles drawn is a_numParticles.
    P is the particle type (must have a RealVect& position() function) and T is an integer type. 
    @param[out] a_particles    Output particles (different for each rank)
    @param[in]  a_numParticles Number of particles
    @param[in]  a_center       Gaussian blob Center
    @param[in]  a_radius       Gaussian blob radius
  */
  template <typename P, typename T>
  static inline void
  drawGaussianParticles(List<P>&       a_particles,
                        const T        a_numParticles,
                        const RealVect a_center,
                        const Real     a_radius) noexcept;

  /*!
    @brief Copy all the particles from the a_src to a_dst
    @param[out] a_dst Copy of original particles. 
    @param[out] a_src Source particles
  */
  template <typename P>
  static inline void
  copy(ParticleContainer<P>& a_dst, const ParticleContainer<P>& a_src) noexcept;

  /*!
    @brief Copy all the particles from the a_src to a_dst. This destroys the source particles
    @param[out] a_dst Copy of original particles. 
    @param[out] a_src Source particles
  */
  template <typename P>
  static inline void
  copyDestructive(ParticleContainer<P>& a_dst, ParticleContainer<P>& a_src) noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_ParticleOpsImplem.H>

#endif
