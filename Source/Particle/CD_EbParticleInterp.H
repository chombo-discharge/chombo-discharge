/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbParticleInterp.H
  @brief  Declaration of a class for handling particle-mesh interpolation and deposition. 
  @author Robert Marskar
*/

#ifndef CD_EbParticleInterp_H
#define CD_EbParticleInterp_H

// CHombo includes
#include <Box.H>
#include <EBISBox.H>
#include <FArrayBox.H>
#include <REAL.H>
#include <Vector.H>
#include <RealVect.H>

// Our includes
#include <CD_DepositionType.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief A class for depositing and interpolating particles. Contains various useful routines for interpolation of different fields in the
  base particle class (e.g., mass, "mobility"). This class also contains functions for depositing particles with 2/4 times the usual particle width,
  which is included because we want to support deposition of particles over refinement boundaries. 
*/
class EbParticleInterp{
public:

  /*!
    @brief Default constructor. Must subsequently call define. 
  */
  EbParticleInterp();

  /*!
    @brief Full constructor
    @param[in] a_region        Cell-centered box
    @param[in] a_ebisbox       EBIS box 
    @param[in] a_dx            Resolution
    @param[in] a_probLo        Lower-left corner of computational domain
    @param[in] a_forceIrregNGP If true, forces NGP deposition in cut-cells. 
  */
  EbParticleInterp(const Box& a_region, const EBISBox& a_ebisbox, const RealVect& a_dx, const RealVect& a_probLo, const bool a_forceIrregNGP);

  /*!
    @brief Define function
    @param[in] a_region        Cell-centered box
    @param[in] a_ebisbox       EBIS box 
    @param[in] a_dx            Resolution
    @param[in] a_probLo        Lower-left corner of computational domain
    @param[in] a_forceIrregNGP If true, forces NGP deposition in cut-cells. 
  */
  void define(const Box& a_region, const EBISBox&  a_ebisbox, const RealVect& a_dx, const RealVect& a_probLo, const bool a_forceIrregNGP);                     

  /*!
    @brief Deposit particle onto the mesh using a standard cloud width. 
    @details This is just like regular particle-mesh deposition. The mesh field should have exactly one component and the quantity to be deposited
    is fetched from the specified template parameter 'a_particleScalarField' This should be a pointer to a member function with signature 
    const Real& P::function() const. E.g. if we are depositing mass through a function const Real& P::mass() const we will specifiy
    
        interpolate<Particle, &Particle::mass>(a_particleList, a_rho, a_interpType)
	
    To deposit a different quantity, simply change the function pointer in the input parameter. 
    @param[in]    a_particleList   Particles to be deposited
    @param[inout] a_rho            Mesh data
    @param[in]    a_depositionType Deposition method
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, const Real&(P::*particleScalarField)() const>
  void deposit(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType) const;

  /*!
    @brief Deposit particle onto the mesh using twice the standard cloud width. 
    @details This is just like regular particle-mesh deposition. The mesh field should have exactly one component and the quantity to be deposited
    is fetched from the specified template parameter 'a_particleScalarField' This should be a pointer to a member function with signature 
    const Real& P::function() const. E.g. if we are depositing mass through a function const Real& P::mass() const we will specifiy
    
        interpolate<Particle, &Particle::mass>(a_particleList, a_rho, a_interpType)
	
    To deposit a different quantity, simply change the function pointer in the input parameter. 
    @param[in]    a_particleList   Particles to be deposited
    @param[inout] a_rho            Mesh data
    @param[in]    a_depositionType Deposition method
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, const Real&(P::*particleScalarField)() const>
  void deposit2(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType) const;

  /*!
    @brief Deposit particle onto the mesh using 4x the standard cloud width. 
    @details This is just like regular particle-mesh deposition. The mesh field should have exactly one component and the quantity to be deposited
    is fetched from the specified template parameter 'a_particleScalarField' This should be a pointer to a member function with signature 
    const Real& P::a_particleScalarField. E.g. if we are depositing mass through a function const Real& P::mass() we will specifiy
    
        interpolate<Particle, &Particle::mass>(a_particleList, a_rho, a_interpType)
	
    To deposit a different quantity, simply change the function pointer in the input parameter. 
    @param[in]    a_particleList   Particles to be deposited
    @param[inout] a_rho            Mesh data
    @param[in]    a_depositionType Deposition method
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, const Real&(P::*particleScalarField)() const>
  void deposit4(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType) const;

  /*!
    @brief Interpolate a vector field onto the particle position. 
    @details This is just like regular particle-mesh interpolation. The input field should have exactly one component and the
    the field will be interpolated onto the template parameter's input field. The template parameter 'particleScalar' should be a pointer
    to a member function which will set the particle field. E.g. the function must have a signature Real& P::particleScalarField(). A valid
    expression is e.g.
    
        interpolate<Particle, &Particle::mass>(a_particleList, a_meshScalarField, a_interpType)
	
    To interpolate onto a different field, e.g. the acceleration field, use interpolate<Particle, &Particle::acceleration>, etc.
    @param[inout] a_particleList    Particles to be interpolated. 
    @param[in]    a_meshScalarField Scalar field on the mesh 
    @param[in]    a_interpType      Interpolation type. 
  */
  template <class P, Real&(P::*particleScalarField)()>
  void interpolate(List<P>&             a_particleList,
		   const FArrayBox&     a_meshScalarField,
		   const DepositionType a_interpType);  

  /*!
    @brief Interpolate a vector field onto the particle position. 
    @details This is just like regular particle-mesh interpolation. The input field should have SpaceDim components and the
    the field will be interpolated onto the template parameter's input field. The template parameter 'particleVectorField' should be a pointer
    to a member function which will set the particle field. E.g. the function must have a signature RealVect& P::particleVectorField(). A valid 
    expression is e.g.
    
        interpolate<Particle, &Particle::velocity>(a_particleList, a_meshVectorField, a_interpType)
	
    To interpolate onto a different field, e.g. the acceleration field, use interpolate<Particle, &Particle::acceleration>, etc.
    @param[inout] a_particleList    Particles to be interpolated. 
    @param[in]    a_meshVectorField Vector field on the mesh
    @param[in]    a_interpType      Interpolation type. 
  */
  template <class P, RealVect&(P::*particleVectorField)()>
  void interpolate(List<P>&             a_particleList,
		   const FArrayBox&     a_meshVectorField,
		   const DepositionType a_interpType);  
  
  /*!
    @brief Deposit particle mass onto the mesh using a standard cloud width. 
    @details This reqires that P has a member function mass()
    @param[in]    a_particleList Particles to be deposited
    @param[inout] a_rho          Mesh data
    @param[in]    a_depositionType   Deposition method
    @note This routine will increment a_rho. 
  */
  template <class P>
  void deposit(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType) const;

  /*!
    @brief Deposit particle mass onto the mesh using twice the cloud width. 
    @details This reqires that P has a member function mass()
    @param[in]    a_particleList Particles to be deposited
    @param[inout] a_rho          Mesh data
    @param[in]    a_depositionType   Deposition method
    @note This routine will increment a_rho. 
  */
  template <class P>
  void deposit2(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType) const;

  /*!
    @brief Deposit particle mass onto the mesh using 4x the cloud width. 
    @details This reqires that P has a member function mass()
    @param[in]    a_particleList Particles to be deposited
    @param[inout] a_rho          Mesh data
    @param[in]    a_depositionType   Deposition method
    @note This routine will increment a_rho. 
  */
  template <class P>
  void deposit4(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType) const;

  /*!
    @brief Interpolate a field onto the particle position. 
    @details This is just like regular particle-mesh interpolation except that it requires that the particle type P 
    has a member function Real& P::diffusion(). The field a_diffusionField must have exactly one component. 
    @param[inout] a_particleList   Particles to be interpolated. 
    @param[in]    a_diffusionField Diffusion "field". I.e. mesh data. 
    @param[in]    a_interpType     Interpolation type. 
  */
  template <class P>
  void interpolateDiffusion(List<P>&             a_particleList,
			    const FArrayBox&     a_diffusionField,
			    const DepositionType a_interpType) const;

  /*!
    @brief Interpolate a field onto the particle position. 
    @details This is just like regular particle-mesh interpolation except that it requires that the particle type P 
    has a member function Real& P::mobility(). The field a_mobilityField must have exact one component. 
    @param[inout] a_particleList   Particles to be interpolated. 
    @param[in]    a_mobilityField  Mobility "field". I.e. mesh data. 
    @param[in]    a_interpType     Interpolation type. 
  */
  template <class P>
  void interpolateMobility(List<P>&             a_particleList,
			   const FArrayBox&     a_mobilityField,
			   const DepositionType a_interpType) const;
  
  /*!
    @brief Interpolate a field onto the particle position. 
    @details This is just like regular particle-mesh interpolation except that it requires that the particle type P 
    has a member function RealVect& P::velocity(). The field a_velocityField must have exactly SpaceDim components. 
    @param[inout] a_particleList   Particles to be interpolated. 
    @param[in]    a_velocityField  Velocity field
    @param[in]    a_interpType     Interpolation type. 
  */
  template <class P>
  void interpolateVelocity(List<P>&             a_particleList,
			   const FArrayBox&     a_velocityField,
			   const DepositionType a_interpType) const;

  /*!
    @brief Interpolate a field onto the particle position. 
    @details This is just like regular particle-mesh interpolation except that it requires that the particle type P 
    has a member function RealVect& P::acceleration(). The field a_accelerationField must have exact SpaceDim components. 
    @param[inout] a_particleList      Particles to be interpolated. 
    @param[in]    a_accelerationField Acceleration field
    @param[in]    a_interpType        Interpolation type. 
  */
  template <class P>
  void interpolateAcceleration(List<P>&             a_particleList,
			       const FArrayBox&     a_accelerationField,
			       const DepositionType a_interpType) const;


protected:

  /*!
    @brief Wrapper function for depositing a single particle.
    @param[inout] a_rho            Mesh data
    @param[in]    a_probLo         Lower-left corner of domain
    @param[in]    a_dx             Grid resolution
    @param[in]    a_position       Particle position
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. 
    @param[in]    a_depositionType Deposition type
  */
  inline
  void depositParticle(FArrayBox&           a_rho,
		       const RealVect&      a_probLo,
		       const RealVect&      a_dx,
		       const RealVect&      a_position,
		       const Real&          a_strength,
		       const DepositionType a_depositionType) const;

  /*!
    @brief Wrapper function for depositing a single particle which has twice the usual cloud width. 
    @param[inout] a_rho            Mesh data
    @param[in]    a_probLo         Lower-left corner of domain
    @param[in]    a_dx             Grid resolution
    @param[in]    a_position       Particle position
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. 
    @param[in]    a_depositionType Deposition type
  */
  inline
  void depositParticle2(FArrayBox&           a_rho,
			const RealVect&      a_probLo,
			const RealVect&      a_dx,
			const RealVect&      a_position,
			const Real&          a_strength,
			const DepositionType a_interpType) const;

  /*!
    @brief Wrapper function for depositing a single particle which has four times the usual cloud width.
    @param[inout] a_rho            Mesh data
    @param[in]    a_probLo         Lower-left corner of domain
    @param[in]    a_dx             Grid resolution
    @param[in]    a_position       Particle position
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. 
    @param[in]    a_depositionType Deposition type
  */
  inline
  void depositParticle4(FArrayBox&           a_rho,
			const RealVect&      a_probLo,
			const RealVect&      a_dx,
			const RealVect&      a_position,
			const Real&          a_strength,
			const DepositionType a_interpType) const;
  
  /*!
    @brief Wrapper function that interpolates a scalar field onto the particle position
    @param[out] a_particleField  Field evaluation at the particle position. 
    @param[in]  a_field          Field. Must have exactly one component. 
    @param[in]  a_probLo         Lower-left corner of domain
    @param[in]  a_dx             Grid resolution
    @param[in]  a_position       Particle position
    @param[in]  a_interpType     Interpolation type. 
  */
  inline
  void interpolateParticle(Real&                a_particleField,
			   const FArrayBox&     a_field,
			   const RealVect&      a_probLo,
			   const RealVect&      a_dx,
			   const RealVect&      a_position,
			   const DepositionType a_interpType) const;

  /*!
    @brief Wrapper function that interpolates a vector field onto the particle position
    @param[out] a_particleField  Field evaluation at the particle position. 
    @param[in]  a_field          Field. Must have exactly one component. 
    @param[in]  a_probLo         Lower-left corner of domain
    @param[in]  a_dx             Grid resolution
    @param[in]  a_position       Particle position
    @param[in]  a_interpType     Interpolation type. 
  */
  inline
  void interpolateParticle(RealVect&            a_particleField,
			   const FArrayBox&     a_field,
			   const RealVect&      a_probLo,
			   const RealVect&      a_dx,
			   const RealVect&      a_position,
			   const DepositionType a_interpType) const;

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
    @brief Special flag which (if true) forces NGP deposition in cut cells.
  */
  bool m_forceIrregNGP;
};

#include <CD_NamespaceFooter.H>

#include <CD_EbParticleInterpImplem.H>

#endif
