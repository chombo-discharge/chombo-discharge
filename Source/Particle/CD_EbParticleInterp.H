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
  @details Currently, this class is cheap to construct and thus does not need to stored for each box during regrids. 
  @note If the user sets the object to force ngp deposition in cut-cells, the result is not kappa-weighted in any form, i.e. the cut-cell will contain
  the full deposited mass (without multiplying or dividing by kappa). Likewise, it will force NGP interpolation in cut-cells. 
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
  */
  EbParticleInterp(const Box& a_region, const EBISBox& a_ebisbox, const RealVect& a_dx, const RealVect& a_probLo);

  /*!
    @brief Define function
    @param[in] a_region        Cell-centered box
    @param[in] a_ebisbox       EBIS box 
    @param[in] a_dx            Resolution
    @param[in] a_probLo        Lower-left corner of computational domain
    @param[in] a_forceIrregNGP If true, forces NGP deposition and interpolation in cut-cells. 
  */
  void define(const Box& a_region, const EBISBox&  a_ebisbox, const RealVect& a_dx, const RealVect& a_probLo);

  /*!
    @brief Deposit particle onto the mesh using a standard cloud width. 
    @details This is just like regular particle-mesh deposition. The mesh field should have exactly one component and the quantity to be deposited
    is fetched from the specified template parameter 'a_particleScalarField' This should be a pointer to a member function with signature 
    const Real& P::function() const. E.g. if we are depositing mass through a function const Real& P::mass() const we will specifiy
    
        deposit<Particle, &Particle::mass>(a_particleList, a_rho, a_depositionType, a_forceIrregNGP)
	
    To deposit a different quantity, simply change the function pointer in the input parameter. 
    @param[in]    a_particleList   Particles to be deposited
    @param[inout] a_rho            Mesh data
    @param[in]    a_depositionType Deposition method
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, const Real&(P::*particleScalarField)() const>
  void deposit(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType, const bool a_forceIrregNGP = false) const;

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
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, const Real&(P::*particleScalarField)() const>
  void deposit2(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType, const bool a_forceIrregNGP = false) const;

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
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, const Real&(P::*particleScalarField)() const>
  void deposit4(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType, const bool a_forceIrregNGP = false) const;

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
		   const DepositionType a_interpType,
		   const bool           a_forceIrregNGP = false) const;  

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
		   const DepositionType a_interpType,
		   const bool           a_forceIrregNGP = false) const;


protected:

  /*!
    @brief Wrapper function for depositing a single particle.
    @param[inout] a_rho            Mesh data
    @param[in]    a_probLo         Lower-left corner of domain
    @param[in]    a_dx             Grid resolution
    @param[in]    a_position       Particle position
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. 
    @param[in]    a_depositionType Deposition type
    @param[in]    a_forceIrregNGP  Force NGP in cut-cells
  */
  inline
  void depositParticle(FArrayBox&           a_rho,
		       const RealVect&      a_probLo,
		       const RealVect&      a_dx,
		       const RealVect&      a_position,
		       const Real&          a_strength,
		       const DepositionType a_depositionType,
		       const bool           a_forceIrregNGP) const;

  /*!
    @brief Wrapper function for depositing a single particle which has twice the usual cloud width. 
    @param[inout] a_rho            Mesh data
    @param[in]    a_probLo         Lower-left corner of domain
    @param[in]    a_dx             Grid resolution
    @param[in]    a_position       Particle position
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. 
    @param[in]    a_depositionType Deposition type
    @param[in]    a_forceIrregNGP  Force NGP in cut-cells
  */
  inline
  void depositParticle2(FArrayBox&           a_rho,
			const RealVect&      a_probLo,
			const RealVect&      a_dx,
			const RealVect&      a_position,
			const Real&          a_strength,
			const DepositionType a_depositionType,						
			const bool           a_forceIrregNGP) const;

  /*!
    @brief Wrapper function for depositing a single particle which has four times the usual cloud width.
    @param[inout] a_rho            Mesh data
    @param[in]    a_probLo         Lower-left corner of domain
    @param[in]    a_dx             Grid resolution
    @param[in]    a_position       Particle position
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. 
    @param[in]    a_depositionType Deposition type
    @param[in]    a_forceIrregNGP  Force NGP in cut-cells
  */
  inline
  void depositParticle4(FArrayBox&           a_rho,
			const RealVect&      a_probLo,
			const RealVect&      a_dx,
			const RealVect&      a_position,
			const Real&          a_strength,
			const DepositionType a_depositionType,			
			const bool           a_forceIrregNGP) const;			
  
  /*!
    @brief Wrapper function that interpolates a scalar field onto the particle position
    @param[out] a_particleField  Field evaluation at the particle position. 
    @param[in]  a_field          Field. Must have exactly one component. 
    @param[in]  a_probLo         Lower-left corner of domain
    @param[in]  a_dx             Grid resolution
    @param[in]  a_position       Particle position
    @param[in]  a_interpType     Interpolation type. 
    @param[in]  a_forceIrregNGP  Force NGP in cut-cells
  */
  inline
  void interpolateParticle(Real&                a_particleField,
			   const FArrayBox&     a_field,
			   const RealVect&      a_probLo,
			   const RealVect&      a_dx,
			   const RealVect&      a_position,
			   const DepositionType a_interpType,
			   const bool           a_forceIrregNGP) const;			  

  /*!
    @brief Wrapper function that interpolates a vector field onto the particle position
    @param[out] a_particleField  Field evaluation at the particle position. 
    @param[in]  a_field          Field. Must have exactly one component. 
    @param[in]  a_probLo         Lower-left corner of domain
    @param[in]  a_dx             Grid resolution
    @param[in]  a_position       Particle position
    @param[in]  a_interpType     Interpolation type. 
    @param[in]  a_forceIrregNGP  Force NGP in cut-cells
  */
  inline
  void interpolateParticle(RealVect&            a_particleField,
			   const FArrayBox&     a_field,
			   const RealVect&      a_probLo,
			   const RealVect&      a_dx,
			   const RealVect&      a_position,
			   const DepositionType a_interpType,
			   const bool           a_forceIrregNGP) const;
  
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
};

#include <CD_NamespaceFooter.H>

#include <CD_EbParticleInterpImplem.H>

#endif
