/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ConductivityEddingtonSP1DomainBc.cpp
  @brief  Implementation of CD_ConductivityEddingtonSP1DomainBc.H
  @author Robert Marskar
*/

// Chombo includes
#include <DirichletConductivityDomainBC.H>
#include <NeumannConductivityDomainBC.H>

// Our includes
#include <CD_ConductivityEddingtonSP1DomainBc.H>
#include <CD_RobinConductivityDomainBc.H>
#include <CD_NamespaceHeader.H>

ConductivityEddingtonSP1DomainBc::ConductivityEddingtonSP1DomainBc(const EddingtonSP1DomainBc& a_domainBc,
								   const LarsenMap&            a_coefficients,
								   const RealVect              a_probLo){
  m_hasCoeff = false;

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const EddingtonSP1DomainBc::Wall curWall       = std::make_pair(dir, sit());
      const EddingtonSP1DomainBc::BcType&     bcType = a_domainBc.getBc(curWall).first;
      const EddingtonSP1DomainBc::BcFunction& bcFunc = a_domainBc.getBc(curWall).second;

      // Generate a BC object using Chombos BC data holders. 
      switch (bcType) {
      case EddingtonSP1DomainBc::BcType::Dirichlet:
	m_conductivityBaseDomainBcObjects.emplace(curWall, std::make_shared<DirichletConductivityDomainBC>());
	break;
      case EddingtonSP1DomainBc::BcType::Neumann:
	m_conductivityBaseDomainBcObjects.emplace(curWall, std::make_shared<NeumannConductivityDomainBC>());
	break;
      case EddingtonSP1DomainBc::BcType::Robin:
	m_conductivityBaseDomainBcObjects.emplace(curWall, std::make_shared<RobinConductivityDomainBc>());
	break;
      default:
	MayDay::Abort("ConductivityEddingtonSP1DomainBc -- unsupported boundary condition passed into ConductivityEddingtonSP1DomainBc");
      }

      // Strange but true thing -- RefCountedPtr is an abomination from before the days of C++11. In this case we pass in
      // std::shared_ptr through EddingtonSP1DomainBc (like God intended), but the Chombo interface uses RefCountedPtr. Both
      // classes will try to delete the bare pointer, so we call neverDelete() to make sure that only std::shared_ptr gets this responsibility. 
      RefCountedPtr<BaseBCFuncEval> func = RefCountedPtr<BaseBCFuncEval> (new EddingtonSP1DomainBcFuncEval(bcFunc, a_probLo));
      func.neverDelete(); 

      std::shared_ptr<ConductivityBaseDomainBC>& curObject = m_conductivityBaseDomainBcObjects.at(curWall);
      curObject->setFunction(func);
    }
  }
}

ConductivityEddingtonSP1DomainBc::~ConductivityEddingtonSP1DomainBc(){

}

void ConductivityEddingtonSP1DomainBc::getFaceFlux(BaseFab<Real>&        a_faceFlux,
						   const BaseFab<Real>&  a_phi,
						   const RealVect&       a_probLo,
						   const RealVect&       a_dx,
						   const int&            a_idir,
						   const Side::LoHiSide& a_side,
						   const DataIndex&      a_dit,
						   const Real&           a_time,
						   const bool&           a_useHomogeneous) {

  if(!m_hasCoeff) this->setCoefficients();
  
  auto& bcPtr = m_conductivityBaseDomainBcObjects.at(std::make_pair(a_idir, a_side));

  bcPtr->getFaceFlux(a_faceFlux,
		     a_phi,
		     a_probLo,
		     a_dx,
		     a_idir,
		     a_side,
		     a_dit,
		     a_time,
		     a_useHomogeneous);
}

void ConductivityEddingtonSP1DomainBc::getFaceFlux(Real&                 a_faceFlux,
						   const VolIndex&       a_vof, 
						   const int&            a_comp, 
						   const EBCellFAB&      a_phi, 
						   const RealVect&       a_probLo, 
						   const RealVect&       a_dx, 
						   const int&            a_idir, 
						   const Side::LoHiSide& a_side, 
						   const DataIndex&      a_dit, 
						   const Real&           a_time, 
						   const bool&           a_useHomogeneous) {

  if(!m_hasCoeff) this->setCoefficients();

  auto& bcPtr = m_conductivityBaseDomainBcObjects.at(std::make_pair(a_idir, a_side));

  bcPtr->getFaceFlux(a_faceFlux,
		     a_vof,
		     a_comp,
		     a_phi,
		     a_probLo,
		     a_dx,
		     a_idir,
		     a_side,
		     a_dit,
		     a_time,
		     a_useHomogeneous);
}
  
void ConductivityEddingtonSP1DomainBc::getFaceGradPhi(Real&                 a_faceFlux,
						      const FaceIndex&      a_face,
						      const int&            a_comp,
						      const EBCellFAB&      a_phi,
						      const RealVect&       a_probLo,
						      const RealVect&       a_dx,
						      const int&            a_idir,
						      const Side::LoHiSide& a_side,
						      const DataIndex&      a_dit,
						      const Real&           a_time,
						      const bool&           a_useAreaFrac,
						      const RealVect&       a_centroid,
						      const bool&           a_useHomogeneous) {

  MayDay::Abort("ConductivityEddingtonSP1DomainBc::getFaceGradPhi -- calling this is an error. How did it get called?");
}


void ConductivityEddingtonSP1DomainBc::fillPhiGhost(FArrayBox&     a_phi,
						    const Box&     a_valid,
						    const Box&     a_domain,
						    Real           a_dx,
						    bool           a_homogeneous) {
  MayDay::Abort("ConductivityEddingtonSP1DomainBc::fillPhiGhost -- calling this is an error. How did it get called?");
}

void ConductivityEddingtonSP1DomainBc::setCoefficients(){
  for (auto& cond : m_conductivityBaseDomainBcObjects){
    cond.second->setCoef(this->m_eblg, this->m_beta, this->m_bcoef);
  }

  m_hasCoeff = true;
}

ConductivityEddingtonSP1DomainBc::EddingtonSP1DomainBcFuncEval::EddingtonSP1DomainBcFuncEval(const EddingtonSP1DomainBc::BcFunction a_bcFunc, const RealVect a_probLo){
  m_bcFunc = a_bcFunc;
  m_probLo = a_probLo;
}

ConductivityEddingtonSP1DomainBc::EddingtonSP1DomainBcFuncEval::~EddingtonSP1DomainBcFuncEval(){
}

Real ConductivityEddingtonSP1DomainBc::EddingtonSP1DomainBcFuncEval::value(const RealVect& a_point, const int& a_comp) const {
  constexpr Real dummyDt = 0.0;
  return m_bcFunc(m_probLo + a_point, dummyDt); // Might seem weird but EbHelmholtzOp does not have access to lower-left corner so we pass that in here. 
}

Real ConductivityEddingtonSP1DomainBc::EddingtonSP1DomainBcFuncEval::derivative(const RealVect& a_point, const int& a_comp, const int& a_derivDir) const {
  MayDay::Abort("EddingtonSP1DomainBcFuncEval::derivative -- calling this is an error. How did you get here?");
}

#include <CD_NamespaceFooter.H>
