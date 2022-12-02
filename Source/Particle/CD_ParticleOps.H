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
    @brief Get the number of physical particles per cell
    @param[out] a_ppc Physical particles per cell
    @param[in]  a_src Particle source container. 
    @details The second template parameter is a pointer to a member function in P that provides the number of physical particles for the computational
    particles. If P has a function const Real& P::weight() const, we can call this as

    ParticleOps::getPhysicalParticlesPerCell<P, &P::weight>(...)
  */
  template <typename P, const Real& (P::*weight)() const>
  static inline void
  getPhysicalParticlesPerCell(EBAMRCellData& a_ppc, const ParticleContainer<P>& a_src) noexcept;

  /*!
    @brief Get the number of computational particles per cell
    @param[out] a_ppc Physical particles per cell
    @param[in]  a_src Particle source container. 
  */
  template <typename P>
  static inline void
  getComputationalParticlesPerCell(EBAMRCellData& a_ppc, const ParticleContainer<P>& a_src) noexcept;

  /*!
    @brief Get the grid cell where the particle lives
    @param[in] a_particlePosition Physical particle position
    @param[in] a_probLo           Lower-left physical corner of simulation domain 
    @param[in] a_dx               Grid resolution.
  */
  static inline IntVect
  getParticleGridCell(const RealVect& a_particlePosition, const RealVect& a_probLo, const RealVect& a_dx) noexcept;

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
