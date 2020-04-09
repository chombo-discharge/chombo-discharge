/*!
  @file   EBMeshInterp.H
  @brief  Currently a shameless rip-off of 
  @author Robert Marskar
  @date   April 2020
  @todo   Remove the dependence of MeshInterp and move the typedef into here when we're ready for that. 
*/

#ifndef _EBMESHINTERP_
#define _EBMESHINTERP_

#include "Box.H"
#include "FArrayBox.H"
#include "REAL.H"
#include "Vector.H"
#include "RealVect.H"

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
class EBMeshInterp{
public:

  /*!
    @brief ctor
  */
  EBMeshInterp();

  /*!
    @brief Full constructor
  */
  EBMeshInterp(const Box&      a_domain,
	       const RealVect& a_dx,
	       const RealVect& a_prob_lo);

  /*!
    @brief Define function
  */
  void define(const Box&      a_domain,
	      const RealVect& a_dx,
	      const RealVect& a_prob_lo);

  /*!
    @brief Deposit particle mass onto the mesh
    @details This reqires that P has a member function mass()
  */
  template <class P>
  void deposit(const List<P>&        a_particleList,
	       FArrayBox&            a_rho,
	       DepositionType::Which a_deposition);

  /*!
    @brief Interpolate the field 
    @details This requires that P has a member function setDiffusion() and that a_diffusionField has one component
  */
  template <class P>
  void interpolateDiffusion(List<P>&              a_particleList,
			    const FArrayBox&      a_diffusionField,
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
  RealVect m_dx;
  RealVect m_prob_lo;
};

#include "EBMeshInterpI.H"

#endif
