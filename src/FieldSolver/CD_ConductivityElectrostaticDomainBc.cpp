/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ConductivityElectrostaticDomainBc.cpp
  @brief  Implementation of CD_ConductivityElectrostaticDomainBc.H
  @author Robert Marskar
*/

// Chombo includes
#include <DirichletConductivityDomainBC.H>
#include <NeumannConductivityDomainBC.H>

// Our includes
#include <CD_ConductivityElectrostaticDomainBc.H>
#include <CD_NamespaceHeader.H>

ConductivityElectrostaticDomainBc::ConductivityElectrostaticDomainBc(const ElectrostaticDomainBc& a_domainBc, const RealVect a_probLo){
  m_hasCoeff = false;

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const ElectrostaticDomainBc::Wall curWall       = std::make_pair(dir, sit());
      const ElectrostaticDomainBc::BcType&     bcType = a_domainBc.getBc(curWall).first;
      const ElectrostaticDomainBc::BcFunction& bcFunc = a_domainBc.getBc(curWall).second;

      // Generate a BC object using Chombos BC data holders. 
      switch (bcType) {
      case ElectrostaticDomainBc::BcType::Dirichlet:
	m_conductivityBaseDomainBcObjects.emplace(curWall, std::make_shared<DirichletConductivityDomainBC>());
	break;
      case ElectrostaticDomainBc::BcType::Neumann:
	m_conductivityBaseDomainBcObjects.emplace(curWall, std::make_shared<NeumannConductivityDomainBC>());
	break;
      default:
	MayDay::Abort("ConductivityElectrostaticDomainBc -- unsupported boundary condition passed into ConductivityElectrostaticDomainBc");
      }

      // Strange but true thing -- RefCountedPtr is an abomination from before the days of C++11. In this case we pass in
      // std::shared_ptr through ElectrostaticDomainBc (like God intended), but the Chombo interface uses RefCountedPtr. Both
      // classes will try to delete the bare pointer, so we call neverDelete() to make sure that only std::shared_ptr gets this responsibility. 
      RefCountedPtr<BaseBCFuncEval> func = RefCountedPtr<BaseBCFuncEval> (new ElectrostaticDomainBcFuncEval(bcFunc, a_probLo));
      func.neverDelete(); 

      std::shared_ptr<ConductivityBaseDomainBC>& curObject = m_conductivityBaseDomainBcObjects.at(curWall);
      curObject->setFunction(func);
    }
  }
}

ConductivityElectrostaticDomainBc::~ConductivityElectrostaticDomainBc(){

}

void ConductivityElectrostaticDomainBc::getFaceFlux(BaseFab<Real>&        a_faceFlux,
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

void ConductivityElectrostaticDomainBc::getFaceFlux(Real&                 a_faceFlux,
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
  
void ConductivityElectrostaticDomainBc::getFaceGradPhi(Real&                 a_faceFlux,
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

  MayDay::Abort("ConductivityElectrostaticDomainBc::getFaceGradPhi -- calling this is an error. How did it get called?");
}


void ConductivityElectrostaticDomainBc::fillPhiGhost(FArrayBox&     a_phi,
						     const Box&     a_valid,
						     const Box&     a_domain,
						     Real           a_dx,
						     bool           a_homogeneous) {
  MayDay::Abort("ConductivityElectrostaticDomainBc::fillPhiGhost -- calling this is an error. How did it get called?");
}

void ConductivityElectrostaticDomainBc::setCoefficients(){
  for (auto& cond : m_conductivityBaseDomainBcObjects){
    cond.second->setCoef(this->m_eblg, this->m_beta, this->m_bcoef);
  }

  m_hasCoeff = true;
}

ConductivityElectrostaticDomainBc::ElectrostaticDomainBcFuncEval::ElectrostaticDomainBcFuncEval(const ElectrostaticDomainBc::BcFunction a_bcFunc, const RealVect a_probLo){
  m_bcFunc = a_bcFunc;
  m_probLo = a_probLo;
}

ConductivityElectrostaticDomainBc::ElectrostaticDomainBcFuncEval::~ElectrostaticDomainBcFuncEval(){
}

Real ConductivityElectrostaticDomainBc::ElectrostaticDomainBcFuncEval::value(const RealVect& a_point, const int& a_comp) const {
  constexpr Real dummyDt = 0.0;
  return m_bcFunc(a_point, dummyDt); // If everything has been done correctly, m_bcFunc will have captured m_time in FieldSolver by reference and not require anything else. 
}

Real ConductivityElectrostaticDomainBc::ElectrostaticDomainBcFuncEval::derivative(const RealVect& a_point, const int& a_comp, const int& a_derivDir) const {
  MayDay::Abort("ElectrostaticDomainBcFuncEval::derivative -- calling this is an error. How did you get here?");
}

#include <CD_NamespaceFooter.H>
