/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBParticleMesh.H
  @brief  Declaration of a class for handling particle-mesh interpolation and deposition. 
  @author Robert Marskar
*/

#ifndef CD_EBParticleMesh_H
#define CD_EBParticleMesh_H

// Std includes
#include <type_traits>

// Chombo includes
#include <EBISBox.H>
#include <EBCellFAB.H>
#include <Vector.H>
#include <RealVect.H>

// Our includes
#include <CD_DepositionType.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief A class for depositing and interpolating particles on a single grid patch. 
  @details Currently, this class is cheap to construct and thus does not need to stored for each box during regrids.
  @note When depositing, this class will never divide by kappa. 
*/
class EBParticleMesh
{
public:
  /*!
    @brief Default constructor. Must subsequently call define. 
  */
  EBParticleMesh();

  /*!
    @brief Full constructor
    @param[in] a_domain        Computational domain
    @param[in] a_region        Cell-centered box
    @param[in] a_ebisbox       EBIS box 
    @param[in] a_dx            Resolution
    @param[in] a_probLo        Lower-left corner of computational domain
  */
  EBParticleMesh(const ProblemDomain& a_domain,
                 const Box&           a_region,
                 const EBISBox&       a_ebisbox,
                 const RealVect&      a_dx,
                 const RealVect&      a_probLo);

  /*!
    @brief Define function
    @param[in] a_domain        Computational domain
    @param[in] a_region        Cell-centered box
    @param[in] a_ebisbox       EBIS box 
    @param[in] a_dx            Resolution
    @param[in] a_probLo        Lower-left corner of computational domain
  */
  void
  define(const ProblemDomain& a_domain,
         const Box&           a_region,
         const EBISBox&       a_ebisbox,
         const RealVect&      a_dx,
         const RealVect&      a_probLo);

  /*!
    @brief Deposit a particle onto the mesh.
    @details The template parameters indicate the particle type and a pointer to a particle function. The user can scale
    the particle width by changing a_widthScale. If this is set to one, the standard width for CIC becomes the grid cell size,
    and the standard width for TSC becomes twice the grid size. 
    @param[inout] a_rho            Mesh data
    @param[in]    a_particles      Particles to be deposited    
    @param[in]    a_depositionType Deposition method
    @param[in]    a_widthScale     Particle scale multiplier.  
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, class Ret, Ret (P::*MemberFunc)() const>
  void
  deposit(EBCellFAB&           a_meshData,
          const List<P>&       a_particles,
          const DepositionType a_depositionType,
          const Real           a_widthScale,
          const bool           a_forceIrregNGP = false) const;

  /*!
    @brief Interpolate mesh data onto a particle. 
    @details The template parameters indicate the particle type and a pointer to a particle function. 
    @param[inout] a_particles      Particles to be deposited
    @param[in]    a_meshData       Mesh data    
    @param[in]    a_depositionType Deposition method
    @param[in]    a_widthScale     Particle scale multiplier.  
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, class Ret, Ret (P::*MemberFunc)()>
  void
  interpolate(List<P>&             a_particles,
              const EBCellFAB&     a_meshData,
              const DepositionType a_interpType,
              const bool           a_forceIrregNGP = false) const;

  /*!
    @brief Sanitize the return function type, ensuring it is either Real or RealVect.
    @return Returns the number of components corresponding to Ret
  */
  template <typename Ret>
  static inline constexpr int
  sanitize() noexcept
  {
    using NoRef = typename std::remove_reference<Ret>::type;
    using Base  = typename std::remove_cv<NoRef>::type;

    static_assert(std::is_same<Base, Real>::value || std::is_same<Base, RealVect>::value,
                  "Ret should be Real or RealVect");

    return std::is_same<Base, Real>::value ? 1 : SpaceDim;
  }

protected:
  /*!
    @brief Problem domain
  */
  ProblemDomain m_domain;

  /*!
    @brief Verbose or not
  */
  bool m_verbose;

  /*!
    @brief Cell-centered box, i.e. valid region. 
  */
  Box m_region;

  /*!
    @brief EBIS box
  */
  EBISBox m_ebisbox;

  /*!
    @brief Grid resolution
  */
  RealVect m_dx;

  /*!
    @brief Lower-left corner of computational domain.
  */
  RealVect m_probLo;

  /*!
    @brief Wrapper function for depositing a single particle with a standard NGP scheme. 
    @param[inout] a_rho            Mesh data
    @param[in]    a_position       Particle position
    @param[in]    a_volumeFactor   Volume factor to multiply by
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. 
    @param[in]    a_numComp        Number of components to deposit
  */
  inline void
  depositParticleNGP(EBCellFAB&      a_rho,
                     const RealVect& a_position,
                     const Real&     a_volumeFactor,
                     const Real*     a_strength,
                     const int&      a_numComp) const noexcept;

  /*!
    @brief Function for depositing a single particle with a CIC scheme.
    @details The user can input a differerent particle width than the standard width byu adjusting a_particleWidth. The input
    box a_particleBox must be large enough to contain the particle. E.g., if a_particleWidth = 1, then the box must be at least
    Box(-IntVect::Unit, IntVect::Unit).
    @param[inout] a_rho            Mesh data
    @param[in]    a_position       Particle position
    @param[in]    a_particleWidth  Particle width (relative to grid resolution)
    @param[in]    a_particleBox    Basic box that is guaranteed to contain the particle.
    @param[in]    a_volumeFactor   Volume factor to multiply by
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. 
    @param[in]    a_numCOmp        Number of components to deposit
    @param[in]    a_forceIrregNGP  Force NGP in cut-cells or not
  */
  inline void
  depositParticleCIC(EBCellFAB&      a_rho,
                     const RealVect& a_position,
                     const RealVect& a_particleWidth,
                     const Box&      a_particleBox,
                     const Real&     a_volumeFactor,
                     const Real*     a_strength,
                     const int&      a_numComp,
                     const bool      a_forceIrregNGP) const noexcept;

  /*!
    @brief Function for depositing a single particle with a TSC scheme.
    @details The user can input a differerent particle width than the standard width byu adjusting a_particleWidth. The input
    box (a_cicBox) must be large enough to contain a particle. E.g., if a_particleWidth = 2, then the box must be at least
    Box(-2*IntVect::Unit, 2*IntVect::Unit).
    @param[inout] a_rho            Mesh data
    @param[in]    a_position       Particle position
    @param[in]    a_particleWidth  Particle width (relative to grid resolution)
    @param[in]    a_particleBox    Basic box that is guaranteed to contain the particle.
    @param[in]    a_volumeFactor   Volume factor to multiply by
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. 
    @param[in]    a_numCOmp        Number of components to deposit
    @param[in]    a_forceIrregNGP  Force NGP in cut-cells or not
  */
  inline void
  depositParticleTSC(EBCellFAB&      a_rho,
                     const RealVect& a_position,
                     const RealVect& a_particleWidth,
                     const Box&      a_particleBox,
                     const Real&     a_volumeFactor,
                     const Real*     a_strength,
                     const int&      a_numComp,
                     const bool      a_forceIrregNGP) const noexcept;

  /*!
    @brief Particle interpolation function using an NGP scheme. 
    @param[out] a_particleField  Particle quantity to interpolate into.
    @param[in]  a_meshData       Mesh data
    @param[in]  a_position       Particle position
    @param[in]  a_numComp        Number of components to interpolate.
  */
  inline void
  interpolateParticleNGP(Real*            a_particleField,
                         const EBCellFAB& a_meshData,
                         const RealVect&  a_position,
                         const int&       a_numComp) const noexcept;

  /*!
    @brief Particle interpolation function using a CIC scheme. 
    @param[out] a_particleField  Particle quantity to interpolate into.
    @param[in]  a_meshData       Mesh data
    @param[in]  a_validRegion    Region where particle can live without kernels reaching outside the domain.
    @param[in]  a_particleBox    Box that is sufficiently large to contain the particle, computed from the clouds lower-left corner
    @param[in]  a_position       Particle position
    @param[in]  a_numComp        Number of components to interpolate.
    @param[in]  a_forceIrregNGP  Force NGP in cut-cells
  */
  inline void
  interpolateParticleCIC(Real*            a_particleField,
                         const EBCellFAB& a_meshData,
                         const Box&       a_validBox,
                         const Box&       a_particleBox,
                         const RealVect&  a_position,
                         const int&       a_numComp,
                         const bool       a_forceIrregNGP) const noexcept;

  /*!
    @brief Particle interpolation function using a TSC scheme. 
    @param[out] a_particleField  Particle quantity to interpolate into.
    @param[in]  a_meshData       Mesh data
    @param[in]  a_validRegion    Region where particle can live without kernels reaching outside the domain.
    @param[in]  a_particleBox    Box that is sufficiently large to contain the particle, computed from the clouds lower-left corner
    @param[in]  a_position       Particle position
    @param[in]  a_numComp        Number of components to interpolate.
    @param[in]  a_forceIrregNGP  Force NGP in cut-cells
  */
  inline void
  interpolateParticleTSC(Real*            a_particleField,
                         const EBCellFAB& a_meshData,
                         const Box&       a_validBox,
                         const Box&       a_particleBox,
                         const RealVect&  a_position,
                         const int&       a_numComp,
                         const bool       a_forceIrregNGP) const noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_EBParticleMeshImplem.H>

#endif
