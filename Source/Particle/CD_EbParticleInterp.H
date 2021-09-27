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
    @brief Deposit particle mass onto the mesh using a standard cloud width. 
    @details This reqires that P has a member function mass()
    @param[in]    a_particleList Particles to be deposited
    @param[inout] a_rho          Mesh data
    @param[in]    a_depositionType   Deposition method
    @note This routine will increment a_rho. 
  */
  template <class P>
  void deposit(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType);

  /*!
    @brief Deposit particle mass onto the mesh using twice the cloud width. 
    @details This reqires that P has a member function mass()
    @param[in]    a_particleList Particles to be deposited
    @param[inout] a_rho          Mesh data
    @param[in]    a_depositionType   Deposition method
    @note This routine will increment a_rho. 
  */
  template <class P>
  void deposit2(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType);

  /*!
    @brief Deposit particle mass onto the mesh using 4x the cloud width. 
    @details This reqires that P has a member function mass()
    @param[in]    a_particleList Particles to be deposited
    @param[inout] a_rho          Mesh data
    @param[in]    a_depositionType   Deposition method
    @note This routine will increment a_rho. 
  */
  template <class P>
  void deposit4(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType);

  /*!
    @brief Interpolate a field onto the particle position. 
    @details This is just like regular particle-mesh interpolation except that it requires that the particle type P 
    has a member function setDiffusion() and that a_diffusionField has exactly one component. 
    @param[inout] a_particleList   Particles to be interpolated. 
    @param[in]    a_diffusionField Diffusion "field". I.e. mesh data. 
    @param[in]    a_interpType     Interpolation type. 
  */
  template <class P>
  void interpolateDiffusion(List<P>&             a_particleList,
			    const FArrayBox&     a_diffusionField,
			    const DepositionType a_interpType);

  /*!
    @brief Interpolate a field onto the particle position. 
    @details This is just like regular particle-mesh interpolation except that it requires that the particle type P 
    has a member function setMobility() and that a_mobilityField has exactly one component. 
    @param[inout] a_particleList   Particles to be interpolated. 
    @param[in]    a_mobilityField  Mobility "field". I.e. mesh data. 
    @param[in]    a_interpType     Interpolation type. 
  */
  template <class P>
  void interpolateMobility(List<P>&             a_particleList,
			   const FArrayBox&     a_mobilityField,
			   const DepositionType a_interpType);
  
  /*!
    @brief Interpolate a field onto the particle position. 
    @details This is just like regular particle-mesh interpolation except that it requires that the particle type P 
    has a member function setVelocity() and that a_velocityField has SpaceDim components. 
    @param[inout] a_particleList   Particles to be interpolated. 
    @param[in]    a_velocityField  Velocity field
    @param[in]    a_interpType     Interpolation type. 
  */
  template <class P>
  void interpolateVelocity(List<P>&             a_particleList,
			   const FArrayBox&     a_velocityField,
			   const DepositionType a_interpType);

  /*!
    @brief Interpolate a field onto the particle position. 
    @details This is just like regular particle-mesh interpolation except that it requires that the particle type P 
    has a member function setAcceleration() and that a_accelerationField has SpaceDim components. 
    @param[inout] a_particleList      Particles to be interpolated. 
    @param[in]    a_accelerationField Acceleration field
    @param[in]    a_interpType        Interpolation type. 
  */
  template <class P>
  void interpolateAcceleration(List<P>&             a_particleList,
			       const FArrayBox&     a_accelerationField,
			       const DepositionType a_interpType);
    
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
  void depositParticle(FArrayBox&           a_rho,
		       const RealVect&      a_probLo,
		       const RealVect&      a_dx,
		       const RealVect&      a_position,
		       const Real&          a_strength,
		       const DepositionType a_depositionType);

  /*!
    @brief Wrapper function for depositing a single particle which has twice the usual cloud width. 
    @param[inout] a_rho            Mesh data
    @param[in]    a_probLo         Lower-left corner of domain
    @param[in]    a_dx             Grid resolution
    @param[in]    a_position       Particle position
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. 
    @param[in]    a_depositionType Deposition type
  */  
  void depositParticle2(FArrayBox&           a_rho,
			const RealVect&      a_probLo,
			const RealVect&      a_dx,
			const RealVect&      a_position,
			const Real&          a_strength,
			const DepositionType a_interpType);

  /*!
    @brief Wrapper function for depositing a single particle which has four times the usual cloud width.
    @param[inout] a_rho            Mesh data
    @param[in]    a_probLo         Lower-left corner of domain
    @param[in]    a_dx             Grid resolution
    @param[in]    a_position       Particle position
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. 
    @param[in]    a_depositionType Deposition type
  */    
  void depositParticle4(FArrayBox&           a_rho,
			const RealVect&      a_probLo,
			const RealVect&      a_dx,
			const RealVect&      a_position,
			const Real&          a_strength,
			const DepositionType a_interpType);
  
  /*!
    @brief Wrapper function that interpolates a scalar field onto the particle position
    @param[out] a_particleField  Field evaluation at the particle position. 
    @param[in]  a_field          Field. Must have exactly one component. 
    @param[in]  a_probLo         Lower-left corner of domain
    @param[in]  a_dx             Grid resolution
    @param[in]  a_position       Particle position
    @param[in]  a_interpType     Interpolation type. 
  */
  void interpolateParticle(Real&                a_particleField,
			   const FArrayBox&     a_field,
			   const RealVect&      a_probLo,
			   const RealVect&      a_dx,
			   const RealVect&      a_position,
			   const DepositionType a_interpType);

  /*!
    @brief Wrapper function that interpolates a vector field onto the particle position
    @param[out] a_particleField  Field evaluation at the particle position. 
    @param[in]  a_field          Field. Must have exactly one component. 
    @param[in]  a_probLo         Lower-left corner of domain
    @param[in]  a_dx             Grid resolution
    @param[in]  a_position       Particle position
    @param[in]  a_interpType     Interpolation type. 
  */
  void interpolateParticle(RealVect&            a_particleField,
			   const FArrayBox&     a_field,
			   const RealVect&      a_probLo,
			   const RealVect&      a_dx,
			   const RealVect&      a_position,
			   const DepositionType a_interpType);

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
