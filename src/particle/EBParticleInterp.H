/*!
  @file   EBParticleInterp.H
  @brief  My clunky attempt at deposition particles in an EB context
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _EBMESHINTERP_
#define _EBMESHINTERP_

#include "Box.H"
#include "EBISBox.H"
#include "FArrayBox.H"
#include "REAL.H"
#include "Vector.H"
#include "RealVect.H"

#include "CD_NamespaceHeader.H"
  
/*!
  @brief Deposition type
*/
namespace DepositionType {
  enum Which{
    NGP,
    CIC,
    TSC,
    W4
  };
};

/*!
  @brief Currently a shameless rip-off of MeshInterp.H
*/
class EBParticleInterp{
public:

  /*!
    @brief ctor
  */
  EBParticleInterp();

  /*!
    @brief Full constructor
  */
  EBParticleInterp(const Box& a_domain, const EBISBox&  a_ebisbox, const RealVect& a_dx, const RealVect& a_prob_lo, const bool a_force_irreg_ngp);

  /*!
    @brief Define function
  */
  void define(const Box& a_domain, const EBISBox&  a_ebisbox, const RealVect& a_dx, const RealVect& a_prob_lo, const bool a_force_irreg_ngp);                     

  /*!
    @brief Deposit particle mass onto the mesh using a standard cloud width. 
    @details This reqires that P has a member function mass()
  */
  template <class P>
  void deposit(const List<P>& a_particleList, FArrayBox&  a_rho, DepositionType::Which a_deposition);

  /*!
    @brief Deposit particle mass onto the mesh using twice the cloud width. 
    @details This reqires that P has a member function mass()
  */
  template <class P>
  void deposit2(const List<P>& a_particleList, FArrayBox& a_rho, DepositionType::Which a_deposition);

  /*!
    @brief Deposit particle mass onto the mesh using 4x the cloud width. 
    @details This reqires that P has a member function mass()
  */
  template <class P>
  void deposit4(const List<P>& a_particleList, FArrayBox& a_rho, DepositionType::Which a_deposition);

  /*!
    @brief Interpolate the field 
    @details This requires that P has a member function setDiffusion() and that a_diffusionField has one component
  */
  template <class P>
  void interpolateDiffusion(List<P>&              a_particleList,
			    const FArrayBox&      a_diffusionField,
			    DepositionType::Which a_deposition);

  /*!
    @brief Interpolate the field 
    @details This requires that P has a member function setMobility() and that a_mobilityField has one component
  */
  template <class P>
  void interpolateMobility(List<P>&              a_particleList,
			   const FArrayBox&      a_mobilityField,
			   DepositionType::Which a_deposition);
  
  /*!
    @brief Interpolate a velocity field
    @details This requires that P has a member function setVelocity() and that a_velocityField has SpaceDim components
  */
  template <class P>
  void interpolateVelocity(List<P>&              a_particleList,
			   const FArrayBox&      a_velocityField,
			   DepositionType::Which a_deposition);

  /*!
    @brief Interpolate an acceleration field
    @details This requires that P has a member function setAcceleration() and that a_accelerationField has SpaceDim components
  */
  template <class P>
  void interpolateAcceleration(List<P>&              a_particleList,
			       const FArrayBox&      a_accelerationField,
			       DepositionType::Which a_deposition);
    
protected:

  /*!
    @brief Wrapper function for depositing a single particle
  */
  void depositParticle(FArrayBox&       a_rho,
		       const RealVect&  a_domainLeftEdge,
		       const RealVect&  a_dx,
		       const RealVect&  a_position,
		       const Real&      a_strength,
		       const DepositionType::Which a_interpType);

  /*!
    @brief Wrapper function for depositing a single particle which has twice the usual cloud width. 
  */
  void depositParticle2(FArrayBox&       a_rho,
			const RealVect&  a_domainLeftEdge,
			const RealVect&  a_dx,
			const RealVect&  a_position,
			const Real&      a_strength,
			const DepositionType::Which a_interpType);

  /*!
    @brief Wrapper function for depositing a single particle which has four times the usual cloud width.
  */
  void depositParticle4(FArrayBox&       a_rho,
			const RealVect&  a_domainLeftEdge,
			const RealVect&  a_dx,
			const RealVect&  a_position,
			const Real&      a_strength,
			const DepositionType::Which a_interpType);
  
  /*!
    @brief Wrapper function that interpolates a scalar field onto the particle position
  */
  void interpolateParticle(Real&             a_particleField,
			   const FArrayBox&  a_field,
			   const RealVect&   a_prob_lo,
			   const RealVect&   a_dx,
			   const RealVect&   a_position,
			   const DepositionType::Which a_deposition);

  /*!
    @brief Wrapper function that interpolates a vector field on the particle position
  */
  void interpolateParticle(RealVect&         a_particleField,
			   const FArrayBox&  a_field,
			   const RealVect&   a_prob_lo,
			   const RealVect&   a_dx,
			   const RealVect&   a_position,
			   const DepositionType::Which a_deposition);
  
  Box      m_domain;
  const EBISBox* m_ebisbox;
  RealVect m_dx;
  RealVect m_prob_lo;
  bool m_irr_ngp; // Force NGP deposition in cut-cells. 
};
#include "CD_NamespaceFooter.H"

#include "EBParticleInterpI.H"

#endif
