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

// CHombo includes
#include <EBISBox.H>
#include <EBCellFAB.H>
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
  template <class P, const Real& (P::*particleScalarField)() const>
  void
  deposit(const List<P>&       a_particleList,
          EBCellFAB&           a_rho,
          const DepositionType a_depositionType,
          const bool           a_forceIrregNGP = false) const;

  /*!
    @brief Deposit particle onto the mesh using a standard cloud width. 
    @details Just like the version above, except that the function signature is Real P::particleScalarField() const.
    @param[in]    a_particleList   Particles to be deposited
    @param[inout] a_rho            Mesh data
    @param[in]    a_depositionType Deposition method
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, Real (P::*particleScalarField)() const>
  void
  deposit(const List<P>&       a_particleList,
          EBCellFAB&           a_rho,
          const DepositionType a_depositionType,
          const bool           a_forceIrregNGP = false) const;

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
  template <class P, const Real& (P::*particleScalarField)() const>
  void
  deposit2(const List<P>&       a_particleList,
           EBCellFAB&           a_rho,
           const DepositionType a_depositionType,
           const bool           a_forceIrregNGP = false) const;

  /*!
    @brief Deposit particle onto the mesh using twice the standard cloud width. 
    @details Justl ike the version above, except that the function signature is Real P::particleScalarField() const.
    @param[in]    a_particleList   Particles to be deposited
    @param[inout] a_rho            Mesh data
    @param[in]    a_depositionType Deposition method
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, Real (P::*particleScalarField)() const>
  void
  deposit2(const List<P>&       a_particleList,
           EBCellFAB&           a_rho,
           const DepositionType a_depositionType,
           const bool           a_forceIrregNGP = false) const;

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
  template <class P, const Real& (P::*particleScalarField)() const>
  void
  deposit4(const List<P>&       a_particleList,
           EBCellFAB&           a_rho,
           const DepositionType a_depositionType,
           const bool           a_forceIrregNGP = false) const;

  /*!
    @brief Deposit particle onto the mesh using 4x the standard cloud width. 
    @details Just like the version above except that the function signature is Real P::particleScalarField() const
    @param[in]    a_particleList   Particles to be deposited
    @param[inout] a_rho            Mesh data
    @param[in]    a_depositionType Deposition method
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, Real (P::*particleScalarField)() const>
  void
  deposit4(const List<P>&       a_particleList,
           EBCellFAB&           a_rho,
           const DepositionType a_depositionType,
           const bool           a_forceIrregNGP = false) const;

  /*!
    @brief Deposit particle onto the mesh using a standard cloud width. 
    @details This is just like regular particle-mesh deposition. The mesh field should have exactly one component and the quantity to be deposited
    is fetched from the specified template parameter 'a_particleVectorField' This should be a pointer to a member function with signature 
    const Real& P::function() const. E.g. if we are depositing mass through a function const Real& P::mass() const we will specifiy
    
    deposit<Particle, &Particle::mass>(a_particleList, a_rho, a_depositionType, a_forceIrregNGP)
	
    To deposit a different quantity, simply change the function pointer in the input parameter. 
    @param[in]    a_particleList   Particles to be deposited
    @param[inout] a_rho            Mesh data
    @param[in]    a_depositionType Deposition method
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, const RealVect& (P::*particleVectorField)() const>
  void
  deposit(const List<P>&       a_particleList,
          EBCellFAB&           a_rho,
          const DepositionType a_depositionType,
          const bool           a_forceIrregNGP = false) const;

  /*!
    @brief Deposit particle onto the mesh using a standard cloud width. 
    @details Just like the version above, except that the function signature is Real P::particleVectorField() const.
    @param[in]    a_particleList   Particles to be deposited
    @param[inout] a_rho            Mesh data
    @param[in]    a_depositionType Deposition method
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, RealVect (P::*particleVectorField)() const>
  void
  deposit(const List<P>&       a_particleList,
          EBCellFAB&           a_rho,
          const DepositionType a_depositionType,
          const bool           a_forceIrregNGP = false) const;

  /*!
    @brief Deposit particle onto the mesh using twice the standard cloud width. 
    @details This is just like regular particle-mesh deposition. The mesh field should have exactly one component and the quantity to be deposited
    is fetched from the specified template parameter 'a_particleVectorField' This should be a pointer to a member function with signature 
    const Real& P::function() const. E.g. if we are depositing mass through a function const Real& P::mass() const we will specifiy
    
    interpolate<Particle, &Particle::mass>(a_particleList, a_rho, a_interpType)
	
    To deposit a different quantity, simply change the function pointer in the input parameter. 
    @param[in]    a_particleList   Particles to be deposited
    @param[inout] a_rho            Mesh data
    @param[in]    a_depositionType Deposition method
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, const RealVect& (P::*particleVectorField)() const>
  void
  deposit2(const List<P>&       a_particleList,
           EBCellFAB&           a_rho,
           const DepositionType a_depositionType,
           const bool           a_forceIrregNGP = false) const;

  /*!
    @brief Deposit particle onto the mesh using twice the standard cloud width. 
    @details Justl ike the version above, except that the function signature is Real P::particleVectorField() const.
    @param[in]    a_particleList   Particles to be deposited
    @param[inout] a_rho            Mesh data
    @param[in]    a_depositionType Deposition method
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, RealVect (P::*particleVectorField)() const>
  void
  deposit2(const List<P>&       a_particleList,
           EBCellFAB&           a_rho,
           const DepositionType a_depositionType,
           const bool           a_forceIrregNGP = false) const;

  /*!
    @brief Deposit particle onto the mesh using 4x the standard cloud width. 
    @details This is just like regular particle-mesh deposition. The mesh field should have exactly one component and the quantity to be deposited
    is fetched from the specified template parameter 'a_particleVectorField' This should be a pointer to a member function with signature 
    const Real& P::a_particleVectorField. E.g. if we are depositing mass through a function const Real& P::mass() we will specifiy
    
    interpolate<Particle, &Particle::mass>(a_particleList, a_rho, a_interpType)
	
    To deposit a different quantity, simply change the function pointer in the input parameter. 
    @param[in]    a_particleList   Particles to be deposited
    @param[inout] a_rho            Mesh data
    @param[in]    a_depositionType Deposition method
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, const RealVect& (P::*particleVectorField)() const>
  void
  deposit4(const List<P>&       a_particleList,
           EBCellFAB&           a_rho,
           const DepositionType a_depositionType,
           const bool           a_forceIrregNGP = false) const;

  /*!
    @brief Deposit particle onto the mesh using 4x the standard cloud width. 
    @details Just like the version above except that the function signature is Real P::particleVectorField() const
    @param[in]    a_particleList   Particles to be deposited
    @param[inout] a_rho            Mesh data
    @param[in]    a_depositionType Deposition method
    @param[in]    a_forceIrregNGP  If true, force NGP in cut-cells
    @note This routine will INCREMENT a_rho. 
  */
  template <class P, RealVect (P::*particleVectorField)() const>
  void
  deposit4(const List<P>&       a_particleList,
           EBCellFAB&           a_rho,
           const DepositionType a_depositionType,
           const bool           a_forceIrregNGP = false) const;

  /*!
    @brief Interpolate a scalar field onto the particle position. 
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
  template <class P, Real& (P::*particleScalarField)()>
  void
  interpolate(List<P>&             a_particleList,
              const EBCellFAB&     a_meshScalarField,
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
  template <class P, RealVect& (P::*particleVectorField)()>
  void
  interpolate(List<P>&             a_particleList,
              const EBCellFAB&     a_meshVectorField,
              const DepositionType a_interpType,
              const bool           a_forceIrregNGP = false) const;

protected:
  /*!
    @brief Wrapper function for depositing a single particle.
    @param[inout] a_rho            Mesh data
    @param[in]    a_probLo         Lower-left corner of domain
    @param[in]    a_dx             Grid resolution
    @param[in]    a_position       Particle position
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. Must be of type Real[] (or &Real if single comopnent)
    @param[in]    a_components     Components to deposit
    @param[in]    a_depositionType Deposition type
    @param[in]    a_forceIrregNGP  Force NGP in cut-cells
  */
  inline void
  depositParticle(EBCellFAB&           a_rho,
                  const RealVect&      a_probLo,
                  const RealVect&      a_dx,
                  const RealVect&      a_position,
                  const Real*          a_strength,
                  const Interval       a_components,
                  const DepositionType a_depositionType,
                  const bool           a_forceIrregNGP) const;

  /*!
    @brief Wrapper function for depositing a single particle which has twice the usual cloud width. 
    @param[inout] a_rho            Mesh data
    @param[in]    a_probLo         Lower-left corner of domain
    @param[in]    a_dx             Grid resolution
    @param[in]    a_position       Particle position
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. Must be of type Real[] (or &Real if single comopnent)
    @param[in]    a_components     Components to deposit
    @param[in]    a_depositionType Deposition type
    @param[in]    a_forceIrregNGP  Force NGP in cut-cells
  */
  inline void
  depositParticle2(EBCellFAB&           a_rho,
                   const RealVect&      a_probLo,
                   const RealVect&      a_dx,
                   const RealVect&      a_position,
                   const Real*          a_strength,
                   const Interval       a_components,
                   const DepositionType a_depositionType,
                   const bool           a_forceIrregNGP) const;

  /*!
    @brief Wrapper function for depositing a single particle which has four times the usual cloud width.
    @param[inout] a_rho            Mesh data
    @param[in]    a_probLo         Lower-left corner of domain
    @param[in]    a_dx             Grid resolution
    @param[in]    a_position       Particle position
    @param[in]    a_strength       Mass (or similar quantity) to be deposited. Must be of type Real[] (or &Real if single comopnent)
    @param[in]    a_components     Components to deposit
    @param[in]    a_depositionType Deposition type
    @param[in]    a_forceIrregNGP  Force NGP in cut-cells
  */
  inline void
  depositParticle4(EBCellFAB&           a_rho,
                   const RealVect&      a_probLo,
                   const RealVect&      a_dx,
                   const RealVect&      a_position,
                   const Real*          a_strength,
                   const Interval       a_components,
                   const DepositionType a_depositionType,
                   const bool           a_forceIrregNGP) const;

  /*!
    @brief Wrapper function that interpolates a particle fields onto positions. 
    @details This switches between various possible kernels, NGP, CIC, TSC, W4. It can interpolate both vectors and scalars, the user inputs the data holder 
    as a Real[] and specifies the components to interpolate. 
    @param[out] a_particleField  Field evaluations at the particle position. Must be of type Real[]. 
    @param[in]  a_field          Field. 
    @param[in]  a_validRegion    Region where particle can live without kernels reaching outside the domain.
    @param[in]  a_probLo         Lower-left corner of domain
    @param[in]  a_dx             Grid resolution
    @param[in]  a_position       Particle position
    @param[in]  a_interval       Components to interpolate
    @param[in]  a_interpType     Interpolation type. 
    @param[in]  a_forceIrregNGP  Force NGP in cut-cells
  */
  inline void
  interpolateParticle(Real*                a_particleField,
                      const EBCellFAB&     a_field,
                      const Box&           a_validBox,
                      const RealVect&      a_probLo,
                      const RealVect&      a_dx,
                      const RealVect&      a_position,
                      const Interval&      a_interval,
                      const DepositionType a_interpType,
                      const bool           a_forceIrregNGP) const;

  /*!
    @brief Problem domain
  */
  ProblemDomain m_domain;

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

#include <CD_EBParticleMeshImplem.H>

#endif
