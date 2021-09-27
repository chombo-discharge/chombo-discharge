/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbParticleInterp.cpp
  @brief  Implementationof CD_EbParticleInterp.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EbParticleInterp.H>
#include <CD_EbParticleInterpF_F.H>
#include <CD_NamespaceHeader.H>

EbParticleInterp::EbParticleInterp(){
  CH_TIME("EbParticleInterp::EbParticleInterp");
}

EbParticleInterp::EbParticleInterp(const Box&      a_region,
				   const EBISBox&  a_ebisbox,
				   const RealVect& a_dx,
				   const RealVect& a_probLo,
				   const bool      a_forceIrregNGP){
  CH_TIME("EbParticleInterp::EbParticleInterp");
  
  this->define(a_region, a_ebisbox, a_dx, a_probLo, a_forceIrregNGP);
}

void EbParticleInterp::define(const Box&      a_region,
			      const EBISBox&  a_ebisbox,
			      const RealVect& a_dx,
			      const RealVect& a_probLo,
			      const bool      a_forceIrregNGP){
  CH_TIME("EbParticleInterp::define");
  
  m_region        = a_region;
  m_ebisbox       = a_ebisbox;
  m_probLo        = a_probLo;
  m_dx            = a_dx;
  m_forceIrregNGP = a_forceIrregNGP;
}

void EbParticleInterp::depositParticle(FArrayBox&           a_rho,
				       const RealVect&      a_probLo,
				       const RealVect&      a_dx,
				       const RealVect&      a_position,
				       const Real&          a_mass,
				       const DepositionType a_depositionType){

  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition. 

  // Check if the cell is irregular. 
  const RealVect rv = (a_position - a_probLo)/a_dx;
  const IntVect  iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));
  const bool isIrregular = m_ebisbox.isIrregular(iv);

  // We can force irregular cells to always do an NGP deposit to prevent clouds leaking into the other side. 
  if(m_forceIrregNGP && isIrregular){
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

void EbParticleInterp::depositParticle2(FArrayBox&           a_rho,
					const RealVect&      a_probLo,
					const RealVect&      a_dx,
					const RealVect&      a_position,
					const Real&          a_mass,
					const DepositionType a_depositionType){

  const RealVect rv = (a_position - a_probLo)/a_dx;
  const IntVect iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

  if(m_ebisbox.isMultiValued(iv)){
    MayDay::Abort("EbParticleInterp::depositParticle2 - multivalued cells not supported (yet)");
  }

  // Irregular cells always do an NGP deposit to prevent clouds leaking into the other side. 
  if(m_ebisbox.isIrregular(iv)){
    FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
			    CHF_CONST_REALVECT(a_probLo),
			    CHF_CONST_REALVECT(a_dx),
			    CHF_CONST_REALVECT(a_position),
			    CHF_CONST_REAL(a_mass));
  }
  else{
    switch (a_depositionType){
    case DepositionType::NGP:
      FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
			      CHF_CONST_REALVECT(a_probLo),
			      CHF_CONST_REALVECT(a_dx),
			      CHF_CONST_REALVECT(a_position),
			      CHF_CONST_REAL(a_mass));
      break;
    case DepositionType::CIC:
      FORT_CIC_DEPOSIT_SCALAR2(CHF_FRA1(a_rho, 0),
			       CHF_CONST_REALVECT(a_probLo),
			       CHF_CONST_REALVECT(a_dx),
			       CHF_CONST_REALVECT(a_position),
			       CHF_CONST_REAL(a_mass));
      break;
    default:
      MayDay::Error("Invalid interpolation type in MeshInterp::depositParticle2 - you may have to use CIC or NGP");
    }
  }
}

void EbParticleInterp::depositParticle4(FArrayBox&           a_rho,
					const RealVect&      a_probLo,
					const RealVect&      a_dx,
					const RealVect&      a_position,
					const Real&          a_mass,
					const DepositionType a_depositionType){

  const RealVect rv = (a_position - a_probLo)/a_dx;
  const IntVect iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

  if(m_ebisbox.isMultiValued(iv)){
    MayDay::Abort("EbParticleInterp::depositParticle2 - multivalued cells not supported (yet)");
  }

  // Irregular cells always do an NGP deposit to prevent clouds leaking into the other side. 
  if(m_ebisbox.isIrregular(iv)){
    FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
			    CHF_CONST_REALVECT(a_probLo),
			    CHF_CONST_REALVECT(a_dx),
			    CHF_CONST_REALVECT(a_position),
			    CHF_CONST_REAL(a_mass));
  }
  else{
    switch (a_depositionType){
    case DepositionType::NGP:
      FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
			      CHF_CONST_REALVECT(a_probLo),
			      CHF_CONST_REALVECT(a_dx),
			      CHF_CONST_REALVECT(a_position),
			      CHF_CONST_REAL(a_mass));
      break;
    case DepositionType::CIC:
      FORT_CIC_DEPOSIT_SCALAR4(CHF_FRA1(a_rho, 0),
			       CHF_CONST_REALVECT(a_probLo),
			       CHF_CONST_REALVECT(a_dx),
			       CHF_CONST_REALVECT(a_position),
			       CHF_CONST_REAL(a_mass));
      break;
    default:
      MayDay::Error("Invalid interpolation type in MeshInterp::depositParticle2 - you may have to use CIC or NGP");
    }
  }
}

void EbParticleInterp::interpolateParticle(Real&                a_particleField,
					   const FArrayBox&     a_field,
					   const RealVect&      a_probLo,
					   const RealVect&      a_dx,
					   const RealVect&      a_position,
					   const DepositionType a_interpType){
  const RealVect rv = (a_position - a_probLo)/a_dx;
  const IntVect iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

  if(m_ebisbox.isMultiValued(iv)){
    MayDay::Abort("EbParticleInterp::depositParticle - multivalued cells not supported (yet)");
  }

  // Irregular cells always do an NGP deposit to prevent clouds leaking into the other side. 
  if(m_ebisbox.isIrregular(iv)){
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
      MayDay::Error("EbParticleInterp::interpolateParticle(RealVect) - Invalid interpolation type.");
    }
  }
}

void EbParticleInterp::interpolateParticle(RealVect&            a_particleField,
					   const FArrayBox&     a_field,
					   const RealVect&      a_probLo,
					   const RealVect&      a_dx,
					   const RealVect&      a_position,
					   const DepositionType a_interpType){
  const RealVect rv = (a_position - a_probLo)/a_dx;
  const IntVect iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

  if(m_ebisbox.isMultiValued(iv)){
    MayDay::Abort("EbParticleInterp::depositParticle - multivalued cells not supported (yet)");
  }

  // Irregular cells always do an NGP deposit to prevent clouds leaking into the other side.
  if(m_ebisbox.isIrregular(iv)){
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
      MayDay::Error("EbParticleInterp::interpolateParticle(RealVect) - Invalid interpolation type.");
    }
  }
}

#include <CD_NamespaceFooter.H>
