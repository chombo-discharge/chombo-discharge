/*!
  @file   EBMeshInterp.cpp
  @brief  Implementationof EBMeshInterp.H
  @author Robert Marskar
  @date   April 2020
*/

#include "EBMeshInterp.H"
#include "EBMeshInterpF_F.H"


EBMeshInterp::EBMeshInterp(){
}

EBMeshInterp::EBMeshInterp(const Box&      a_domain,
			   const RealVect& a_dx,
			   const RealVect& a_prob_lo){
  this->define(a_domain, a_dx, a_prob_lo);
}

void EBMeshInterp::define(const Box&      a_domain,
			  const RealVect& a_dx,
			  const RealVect& a_prob_lo){
  m_domain  = a_domain;
  m_prob_lo = a_prob_lo;
  m_dx      = a_dx;
}

void EBMeshInterp::depositParticle(FArrayBox&       a_rho,
				   const RealVect&  a_prob_lo,
				   const RealVect&  a_dx,
				   const RealVect&  a_position,
				   const Real&      a_mass,
				   const DepositionType::Which a_interpType){
  switch (a_interpType){
  case DepositionType::NGP:
    FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
			    CHF_CONST_REALVECT(a_prob_lo),
			    CHF_CONST_REALVECT(a_dx),
			    CHF_CONST_REALVECT(a_position),
			    CHF_CONST_REAL(a_mass));
    break;
  case DepositionType::CIC:
    FORT_CIC_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
			    CHF_CONST_REALVECT(a_prob_lo),
			    CHF_CONST_REALVECT(a_dx),
			    CHF_CONST_REALVECT(a_position),
			    CHF_CONST_REAL(a_mass));
    break;
  case DepositionType::TSC:
    FORT_TSC_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
			    CHF_CONST_REALVECT(a_prob_lo),
			    CHF_CONST_REALVECT(a_dx),
			    CHF_CONST_REALVECT(a_position),
			    CHF_CONST_REAL(a_mass));
    break;
  case DepositionType::W4:
    FORT_W4_DEPOSIT_SCALAR(CHF_FRA1(a_rho, 0),
			   CHF_CONST_REALVECT(a_prob_lo),
			   CHF_CONST_REALVECT(a_dx),
			   CHF_CONST_REALVECT(a_position),
			   CHF_CONST_REAL(a_mass));
    break;
  default:
    MayDay::Error("Invalid interpolation type in MeshInterp::depositParticle");
  }
}

void EBMeshInterp::interpolateParticle(Real&             a_particleField,
				       const FArrayBox&  a_field,
				       const RealVect&   a_prob_lo,
				       const RealVect&   a_dx,
				       const RealVect&   a_position,
				       const DepositionType::Which a_interpType){
  switch (a_interpType) {
  case DepositionType::NGP:
    FORT_NGP_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				CHF_CONST_FRA1(a_field,0),
				CHF_CONST_REALVECT(a_prob_lo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position));

    break;
  case DepositionType::CIC:
    FORT_CIC_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				CHF_CONST_FRA1(a_field,0),
				CHF_CONST_REALVECT(a_prob_lo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position));
    break;
  case DepositionType::TSC:
    FORT_TSC_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				CHF_CONST_FRA1(a_field,0),
				CHF_CONST_REALVECT(a_prob_lo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position));
    break;
  case DepositionType::W4:
    FORT_W4_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
			       CHF_CONST_FRA1(a_field,0),
			       CHF_CONST_REALVECT(a_prob_lo),
			       CHF_CONST_REALVECT(a_dx),
			       CHF_CONST_REALVECT(a_position));
    break;
  default:
    MayDay::Error("EBMeshInterp::interpolateParticle(RealVect) - Invalid interpolation type.");
  }
}

void EBMeshInterp::interpolateParticle(RealVect&         a_particleField,
				       const FArrayBox&  a_field,
				       const RealVect&   a_prob_lo,
				       const RealVect&   a_dx,
				       const RealVect&   a_position,
				       const DepositionType::Which a_interpType){
  switch (a_interpType) {
  case DepositionType::NGP:
    FORT_NGP_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				CHF_CONST_FRA(a_field),
				CHF_CONST_REALVECT(a_prob_lo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position));

    break;
  case DepositionType::CIC:
    FORT_CIC_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				CHF_CONST_FRA(a_field),
				CHF_CONST_REALVECT(a_prob_lo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position));
    break;
  case DepositionType::TSC:
    FORT_TSC_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				CHF_CONST_FRA(a_field),
				CHF_CONST_REALVECT(a_prob_lo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position));
    break;
  case DepositionType::W4:
    FORT_W4_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
			       CHF_CONST_FRA(a_field),
			       CHF_CONST_REALVECT(a_prob_lo),
			       CHF_CONST_REALVECT(a_dx),
			       CHF_CONST_REALVECT(a_position));
    break;
  default:
    MayDay::Error("EBMeshInterp::interpolateParticle(RealVect) - Invalid interpolation type.");
  }
}
  
