/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleOps.H
  @brief  Declaration of a static class containing some common useful particle routines that would otherwise be replicated across classes. 
  @author Robert Marskar
  @todo   In ebIntersectionRaycast we should check if the path is tangential to the EB since that can consume a lot of resources. 
*/

#ifndef CD_PARTICLEOPS_H
#define CD_PARTICLEOPS_H

// Chombo includes
#include <RefCountedPtr.H>
#include <BaseIF.H>
#include <RealVect.H>

// Our includes
#include <CD_EBAMRData.H>
#include <CD_ParticleMap.H>
#include <CD_ParticleContainer.H>
#include <CD_NamespaceHeader.H>

/**
  @brief A shallow static class for doing various kinds of particle-related operations.
*/
class ParticleOps
{
public:
  /**
    @brief Get the cell index corresponding to the particle position
    @param[in] a_particlePosition Particle position
    @param[in] a_probLo Lower-left corner of computational domain
    @param[in] a_dx Grid resolution
    @return Particle cell index
  */
  static inline IntVect
  getParticleCellIndex(const RealVect& a_particlePosition, const RealVect& a_probLo, const Real& a_dx) noexcept;

  /**
    @brief Get the cell index corresponding to the particle position
    @param[in] a_particlePosition Particle position
    @param[in] a_probLo Lower-left corner of computational domain
    @param[in] a_dx Grid resolution
    @return Particle cell index
  */
  static inline IntVect
  getParticleCellIndex(const RealVect& a_particlePosition, const RealVect& a_probLo, const RealVect& a_dx) noexcept;

  /**
    @brief Get the number of physical particles per cell (SoA overload).
    @details Sums the container-owned weight column into a_ppc, deposited NGP by grid cell.
    @param[out] a_ppc Physical particles per cell.
    @param[in]  a_src SoA particle source container.
    @tparam P      Payload type.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  static inline void
  getPhysicalParticlesPerCell(EBAMRCellData& a_ppc, const ParticleContainer<P, Traits>& a_src) noexcept;

  /**
    @brief Get the number of computational particles per cell (SoA overload).
    @param[out] a_ppc Computational particles per cell.
    @param[in]  a_src SoA particle source container.
    @tparam P      Payload type.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  static inline void
  getComputationalParticlesPerCell(EBAMRCellData& a_ppc, const ParticleContainer<P, Traits>& a_src) noexcept;

  /**
    @brief Get the grid cell where the particle lives
    @param[in] a_particlePosition Physical particle position
    @param[in] a_probLo           Lower-left physical corner of simulation domain 
    @param[in] a_dx               Grid resolution.
    @return Particle grid cell
  */
  static inline IntVect
  getParticleGridCell(const RealVect& a_particlePosition, const RealVect& a_probLo, const RealVect& a_dx) noexcept;

  /**
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

  /**
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

  /**
    @brief Compute the intersection point between a particle path and an implicit function using a ray-casting algorithm.
    @details This routine will assume that a_impcFunc is an approximation to the signed distance function and compute the intersection point using ray-marching. Starting
    on the starting position we compute the distance to the EB and move the particle the compute distance along the particle trajectory. If the particle comes too close
    to the EB (within a_tolerance) we consider it to be absorbed. In that case the output argument a_s determines the intersection point by x(s) = x0 + s*(x1-x0) where 
    x0 is the starting position (a_oldPos), x1 is the end position (a_newPos), and x(s) is the intersection point with the geometry. 
    @param[in]  a_impFunc    Implicit function. 
    @param[in]  a_oldPos     Particle starting position
    @param[in]  a_newPos     Particle end position
    @param[out] a_s          Relative length along the path 
    @return Returns true if the particle crossed into the EB.
    @param[in] a_tolerance Tolerance
  */
  static inline bool
  ebIntersectionRaycast(const RefCountedPtr<BaseIF>& a_impFunc,
                        const RealVect&              a_oldPos,
                        const RealVect&              a_newPos,
                        const Real&                  a_tolerance,
                        Real&                        a_s);

  /**
    @brief Move all particles from a_src into a_dst (per leaf), emptying a_src. SoA overload.
    @details Both containers must be defined over the same grids. Per patch, a_dst's leaf is replaced by
    a_src's (via clear + catenate -- an O(1) arena swap when a_dst's leaf is empty), and a_src is emptied.
    @param[in,out] a_dst Destination container (its leaves are overwritten).
    @param[in,out] a_src Source container (emptied).
    @tparam P      Payload type.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  static inline void
  copyDestructive(ParticleContainer<P, Traits>& a_dst, ParticleContainer<P, Traits>& a_src) noexcept;

  /**
    @brief Global sum of the container-owned weight column (SoA overload).
    @param[in] a_particles SoA particles.
    @return Global (MPI-reduced) sum of all valid particle weights.
    @tparam P      Payload type.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  static inline Real
  sum(const ParticleContainer<P, Traits>& a_particles) noexcept;

  /**
    @brief Set value function for SoA containers. Lets the user set particle parameters via a (leaf, index) callback.
    @param[in,out] a_particles Particles
    @param[in]    a_functor   Callback (leaf, index) for setting particle values via payload columns
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  static inline void
  setData(ParticleContainer<P, Traits>&                                    a_particles,
          const std::function<void(ParticleSoA<P, Traits>&, std::size_t)>& a_functor) noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_ParticleOpsImplem.H>

#endif
