/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbParticleInterpImplem.H
  @brief  Implementatoin of CD_EbParticleInterp.H
  @author Robert Marskar
*/

#ifndef CD_EbParticleInterpImplem_H
#define CD_EbParticleInterpImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EbParticleInterp.H>
#include <CD_EbParticleInterpF_F.H>
#include <CD_NamespaceHeader.H>

inline
void EbParticleInterp::depositParticle(FArrayBox&           a_rho,
				       const RealVect&      a_probLo,
				       const RealVect&      a_dx,
				       const RealVect&      a_position,
				       const Real&          a_mass,
				       const DepositionType a_depositionType) const {

  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition. 

  // Check if the cell is irregular. 
  const RealVect rv = (a_position - a_probLo)/a_dx;
  const IntVect  iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));
  const bool isIrregular = m_ebisbox.isIrregular(iv);

  // We can force NGP deposition in cut-cells if we want. 
  if(m_ebisbox.isIrregular(iv) && m_forceIrregNGP){
    FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
			    CHF_CONST_REALVECT(a_probLo),
			    CHF_CONST_REALVECT(a_dx),
			    CHF_CONST_REALVECT(a_position),
			    CHF_CONST_REAL(a_mass));
  }
  else{
    switch (a_depositionType){
    case DepositionType::NGP:
      {
	FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position),
				CHF_CONST_REAL(a_mass));
	break;
      }
    case DepositionType::CIC:
      {
	FORT_CIC_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position),
				CHF_CONST_REAL(a_mass));
	break;
      }
    case DepositionType::TSC:
      {
	FORT_TSC_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position),
				CHF_CONST_REAL(a_mass));
	break;
      }
    case DepositionType::W4:
      {
	FORT_W4_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
			       CHF_CONST_REALVECT(a_probLo),
			       CHF_CONST_REALVECT(a_dx),
			       CHF_CONST_REALVECT(a_position),
			       CHF_CONST_REAL(a_mass));
	break;
      }
    default:
      {
	MayDay::Error("EbParticleInterp::depositParticle - logic bust, unknown particle deposition.");
	break;
      }
    }
  }
}

inline
void EbParticleInterp::depositParticle2(FArrayBox&           a_rho,
					const RealVect&      a_probLo,
					const RealVect&      a_dx,
					const RealVect&      a_position,
					const Real&          a_mass,
					const DepositionType a_depositionType) const {

  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition.
  //
  //       Note that this is the version which deposits with particle widths that are 4 times the "usual" width. Currently, only NGP
  //       and CIC is supported. 
  
  const RealVect rv = (a_position - a_probLo)/a_dx;
  const IntVect iv  = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

  // Force NGP in cut-cells if we want. 
  if(m_ebisbox.isIrregular(iv) && m_forceIrregNGP){
    FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
			    CHF_CONST_REALVECT(a_probLo),
			    CHF_CONST_REALVECT(a_dx),
			    CHF_CONST_REALVECT(a_position),
			    CHF_CONST_REAL(a_mass));
  }
  else{
    switch (a_depositionType){
    case DepositionType::NGP:
      {
	FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position),
				CHF_CONST_REAL(a_mass));
	break;
      }
    case DepositionType::CIC:
      {
	FORT_CIC_DEPOSIT_SCALAR2(CHF_FRA1(a_rho, 0),
				 CHF_CONST_REALVECT(a_probLo),
				 CHF_CONST_REALVECT(a_dx),
				 CHF_CONST_REALVECT(a_position),
				 CHF_CONST_REAL(a_mass));
	break;
      }
    default:
      {
	MayDay::Error("EbParticleInterp::depositParticle2 - Invalid deposition type - only NGP and CIC supported for this deposition method.");
	break;
      }
    }
  }
}

inline
void EbParticleInterp::depositParticle4(FArrayBox&           a_rho,
					const RealVect&      a_probLo,
					const RealVect&      a_dx,
					const RealVect&      a_position,
					const Real&          a_mass,
					const DepositionType a_depositionType) const {

  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition.
  //
  //       Note that this is the version which deposits with particle widths that are 4 times the "usual" width. Currently, only NGP
  //       and CIC is supported. 

  const RealVect rv = (a_position - a_probLo)/a_dx;
  const IntVect iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

  // We can force NGP in cut-cells if we want. 
  if(m_ebisbox.isIrregular(iv) && m_forceIrregNGP){
    FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
			    CHF_CONST_REALVECT(a_probLo),
			    CHF_CONST_REALVECT(a_dx),
			    CHF_CONST_REALVECT(a_position),
			    CHF_CONST_REAL(a_mass));
  }
  else{
    switch (a_depositionType){
    case DepositionType::NGP:
      {
	FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position),
				CHF_CONST_REAL(a_mass));
	break;
      }
    case DepositionType::CIC:
      {
	FORT_CIC_DEPOSIT_SCALAR4(CHF_FRA1(a_rho, 0),
				 CHF_CONST_REALVECT(a_probLo),
				 CHF_CONST_REALVECT(a_dx),
				 CHF_CONST_REALVECT(a_position),
				 CHF_CONST_REAL(a_mass));
	break;
      }
    default:
      {
	MayDay::Error("EbParticleInterp::depositParticle4 - Invalid deposition type - only NGP and CIC supported for this deposition method.");
	break;
      }
    }
  }
}

inline
void EbParticleInterp::interpolateParticle(Real&                a_particleField,
					   const FArrayBox&     a_field,
					   const RealVect&      a_probLo,
					   const RealVect&      a_dx,
					   const RealVect&      a_position,
					   const DepositionType a_interpType) const {
  const RealVect rv = (a_position - a_probLo)/a_dx;
  const IntVect  iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

  // Irregular cells always do an NGP deposit to prevent clouds leaking into the other side. 
  if(m_ebisbox.isIrregular(iv) && m_forceIrregNGP){
    FORT_NGP_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				CHF_CONST_FRA1(a_field,0),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position));
  }
  else{
    switch (a_interpType) {
    case DepositionType::NGP:
      FORT_NGP_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				  CHF_CONST_FRA1(a_field,0),
				  CHF_CONST_REALVECT(a_probLo),
				  CHF_CONST_REALVECT(a_dx),
				  CHF_CONST_REALVECT(a_position));

      break;
    case DepositionType::CIC:
      FORT_CIC_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				  CHF_CONST_FRA1(a_field,0),
				  CHF_CONST_REALVECT(a_probLo),
				  CHF_CONST_REALVECT(a_dx),
				  CHF_CONST_REALVECT(a_position));
      break;
    case DepositionType::TSC:
      FORT_TSC_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				  CHF_CONST_FRA1(a_field,0),
				  CHF_CONST_REALVECT(a_probLo),
				  CHF_CONST_REALVECT(a_dx),
				  CHF_CONST_REALVECT(a_position));
      break;
    case DepositionType::W4:
      FORT_W4_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				 CHF_CONST_FRA1(a_field,0),
				 CHF_CONST_REALVECT(a_probLo),
				 CHF_CONST_REALVECT(a_dx),
				 CHF_CONST_REALVECT(a_position));
      break;
    default:
      MayDay::Error("EbParticleInterp::interpolateParticle(Real) - Invalid interpolation type requested.");
    }
  }
}

inline
void EbParticleInterp::interpolateParticle(RealVect&            a_particleField,
					   const FArrayBox&     a_field,
					   const RealVect&      a_probLo,
					   const RealVect&      a_dx,
					   const RealVect&      a_position,
					   const DepositionType a_interpType) const {
  const RealVect rv = (a_position - a_probLo)/a_dx;
  const IntVect  iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

  // Irregular cells always do an NGP deposit to prevent clouds leaking into the other side.
  if(m_ebisbox.isIrregular(iv) && m_forceIrregNGP){    
    FORT_NGP_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				CHF_CONST_FRA(a_field),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position));
  }
  else{
    switch (a_interpType) {
    case DepositionType::NGP:
      FORT_NGP_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				  CHF_CONST_FRA(a_field),
				  CHF_CONST_REALVECT(a_probLo),
				  CHF_CONST_REALVECT(a_dx),
				  CHF_CONST_REALVECT(a_position));

      break;
    case DepositionType::CIC:
      FORT_CIC_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				  CHF_CONST_FRA(a_field),
				  CHF_CONST_REALVECT(a_probLo),
				  CHF_CONST_REALVECT(a_dx),
				  CHF_CONST_REALVECT(a_position));
      break;
    case DepositionType::TSC:
      FORT_TSC_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				  CHF_CONST_FRA(a_field),
				  CHF_CONST_REALVECT(a_probLo),
				  CHF_CONST_REALVECT(a_dx),
				  CHF_CONST_REALVECT(a_position));
      break;
    case DepositionType::W4:
      FORT_W4_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				 CHF_CONST_FRA(a_field),
				 CHF_CONST_REALVECT(a_probLo),
				 CHF_CONST_REALVECT(a_dx),
				 CHF_CONST_REALVECT(a_position));
      break;
    default:
      MayDay::Error("EbParticleInterp::interpolateParticle(RealVect) - invalid interpolation type requested.");
    }
  }
}

template <class P, const Real&(P::*particleScalarField)() const>
void EbParticleInterp::deposit(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType) const {
  CH_TIME("EbParticleInterp::deposit");

  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();    
    
    this->depositParticle(a_rho, 
			  m_probLo, 
			  m_dx, 
			  curPosition, 
			  curStrength, 
			  a_depositionType);
  }  
}

template <class P, const Real&(P::*particleScalarField)() const>
void EbParticleInterp::deposit2(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType) const {
  CH_TIME("EbParticleInterp::deposit");

  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();    

    this->depositParticle2(a_rho, 
			   m_probLo, 
			   m_dx, 
			   curPosition, 
			   curStrength, 
			   a_depositionType);
  }  
}

template <class P, const Real&(P::*particleScalarField)() const>
void EbParticleInterp::deposit4(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType) const {
  CH_TIME("EbParticleInterp::deposit");

  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();    
    
    this->depositParticle4(a_rho, 
			   m_probLo, 
			   m_dx, 
			   curPosition, 
			   curStrength, 
			   a_depositionType);
  }  
}

template <class P, Real&(P::*particleScalarField)()>
void EbParticleInterp::interpolate(List<P>&             a_particleList,
				   const FArrayBox&     a_meshScalarField,
				   const DepositionType a_interpType){

  CH_assert(a_meshScalarField.nComp() == 1);

  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    P&              curParticle  = lit();
    const RealVect& curPosition  = curParticle.position();

    // Note: We MUST initialize = 0 because the higher-order kernels do not reset the field. It doesn't matter if we do it here
    //       or in the Fortran kernel.     
    Real& curParticleField = (curParticle.*particleScalarField)();
    curParticleField = 0.0;
    this->interpolateParticle(curParticleField,
			      a_meshScalarField,
			      m_probLo,
			      m_dx,
			      curPosition,
			      a_interpType);
  }
}

template <class P, RealVect&(P::*particleVectorField)()>
void EbParticleInterp::interpolate(List<P>&             a_particleList,
				   const FArrayBox&     a_meshVectorField,
				   const DepositionType a_interpType){

  CH_assert(a_meshVectorField.nComp() == SpaceDim);

  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    P&              curParticle  = lit();
    const RealVect& curPosition  = curParticle.position();

    // Note: We MUST initialize = 0 because the higher-order kernels do not reset the field. It doesn't matter if we do it here
    //       or in the Fortran kernel.     
    RealVect& curParticleField = (curParticle.*particleVectorField)();
    curParticleField = RealVect::Zero;
    this->interpolateParticle(curParticleField,
			      a_meshVectorField,
			      m_probLo,
			      m_dx,
			      curPosition,
			      a_interpType);
  }
}
  
template <class P>
void EbParticleInterp::deposit(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType) const {
  CH_TIME("EbParticleInterp::deposit");
  
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curMass     = curParticle.mass();
    
    this->depositParticle(a_rho, 
			  m_probLo, 
			  m_dx, 
			  curPosition, 
			  curMass, 
			  a_depositionType);
  }
}

template <class P>
void EbParticleInterp::deposit2(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType) const {
  CH_TIME("EbParticleInterp::deposit2");
  
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curMass     = curParticle.mass();
    
    this->depositParticle2(a_rho, 
			   m_probLo, 
			   m_dx, 
			   curPosition, 
			   curMass, 
			   a_depositionType);
  }
}

template <class P>
void EbParticleInterp::deposit4(const List<P>& a_particleList, FArrayBox& a_rho, const DepositionType a_depositionType) const {
  CH_TIME("EbParticleInterp::deposit4");
  
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curMass     = curParticle.mass();
    
    this->depositParticle4(a_rho, 
			   m_probLo, 
			   m_dx, 
			   curPosition, 
			   curMass, 
			   a_depositionType);
  }
}



template <class P>
void EbParticleInterp::interpolateDiffusion(List<P>&             a_particleList,
					    const FArrayBox&     a_diffusionField,
					    const DepositionType a_interpType) const {
  CH_TIME("EbParticleInterp::interpolationDiffusion");

  CH_assert(a_diffusionField.nComp() == 1);
  
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    P&              curParticle  = lit();
    Real&           curDiffusion = curParticle.diffusion();
    const RealVect& curPosition  = curParticle.position();

    // Note: We MUST initialize = 0 because the higher-order kernels do not reset the field. It doesn't matter if we do it here
    //       or in the Fortran kernel.     
    curDiffusion = 0.0;
    this->interpolateParticle(curDiffusion,
			      a_diffusionField,
			      m_probLo,
			      m_dx,
			      curPosition,
			      a_interpType);
  }
}

template <class P>
void EbParticleInterp::interpolateMobility(List<P>&             a_particleList,
					   const FArrayBox&     a_mobilityField,
					   const DepositionType a_interpType) const {
  CH_TIME("EbParticleInterp::interpolateMobility");
  
  CH_assert(a_mobilityField.nComp() == 1);
  
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    P&              curParticle = lit();
    Real&           curMobility = curParticle.mobility();
    const RealVect& curPosition = curParticle.position();

    // Note: We MUST initialize = 0 because the higher-order kernels do not reset the field. It doesn't matter if we do it here
    //       or in the Fortran kernel. 
    curMobility = 0.0;    
    this->interpolateParticle(curMobility,
			      a_mobilityField,
			      m_probLo,
			      m_dx,
			      curPosition,
			      a_interpType);
  }
}
  
template <class P>
void EbParticleInterp::interpolateVelocity(List<P>&             a_particleList,
					   const FArrayBox&     a_velocityField,
					   const DepositionType a_interpType) const {
  CH_TIME("EbParticleInterp::interpolateVelocity");

  CH_assert(a_velocityField.nComp() == SpaceDim);
  
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    P&              curParticle = lit();
    RealVect&       curVelocity = curParticle.velocity();
    const RealVect& curPosition = curParticle.position();

    // Note: We MUST initialize = 0 because the higher-order kernels do not reset the field. It doesn't matter if we do it here
    //       or in the Fortran kernel.         
    curVelocity = RealVect::Zero;
    this->interpolateParticle(curVelocity,
			      a_velocityField,
			      m_probLo,
			      m_dx,
			      curPosition,
			      a_interpType);
  }
}

template <class P>
void EbParticleInterp::interpolateAcceleration(List<P>&             a_particleList,
					       const FArrayBox&     a_accelerationField,
					       const DepositionType a_interpType) const {
  CH_TIME("EbParticleInterp::interpolateAcceleration");

  CH_assert(a_accelerationField.nComp() == SpaceDim);
  
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    P&              curParticle     = lit();
    RealVect&       curAcceleration = curParticle.acceleration();
    const RealVect& curPosition     = curParticle.position();

    // Note: We MUST initialize = 0 because the higher-order kernels do not reset the field. It doesn't matter if we do it here
    //       or in the Fortran kernel.             
    curAcceleration = RealVect::Zero;
    this->interpolateParticle(curAcceleration,
			      a_accelerationField,
			      m_probLo,
			      m_dx,
			      curPosition,
			      a_interpType);
  }
}

#include <CD_NamespaceFooter.H>

#endif
