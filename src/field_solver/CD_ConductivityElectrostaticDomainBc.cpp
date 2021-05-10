/*!
  @file   CD_ConductivityElectrostaticDomainBc.cpp
  @brief  Implementation of CD_ConductivityElectrostaticDomainBc.H
  @author Robert Marskar
  @date   May 2021
*/

#include <DirichletConductivityDomainBC.H>
#include <NeumannConductivityDomainBC.H>

#include "CD_ConductivityElectrostaticDomainBc.H"
#include "CD_NamespaceHeader.H"

ConductivityElectrostaticDomainBc::ConductivityElectrostaticDomainBc(const WallBcFuncs& a_bcFunctions, const WallBcTypes& a_bcTypes){
  m_hasCoeff = false;

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const ElectrostaticDomainBc::Wall curWall = std::make_pair(dir, sit());

      // Make the bounadry condition type. 
      const auto& bcType = a_bcTypes.at(curWall);
      switch (bcType) {
      case ElectrostaticDomainBc::BcType::Dirichlet:
	m_conductivityBaseDomainBcObjects.emplace(curWall, std::shared_ptr<ConductivityBaseDomainBC>(new DirichletConductivityDomainBC()));
	break;
      case ElectrostaticDomainBc::BcType::Neumann:
	m_conductivityBaseDomainBcObjects.emplace(curWall, std::shared_ptr<ConductivityBaseDomainBC>(new NeumannConductivityDomainBC()));
	break;
      default:
	MayDay::Abort("ConductivityElectrostaticDomainBc -- unsupported boundary condition passed into ConductivityElectrostaticDomainBc");
      }

      // Set the associated pointer.
      const std::shared_ptr<ElectrostaticDomainBcFuncEval>& curFunc = a_bcFunctions.at(curWall);
      std::shared_ptr<ConductivityBaseDomainBC>& curObject    = m_conductivityBaseDomainBcObjects.at(curWall);

      curObject->setFunction(RefCountedPtr<BaseBCFuncEval> (&(*curFunc)));
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

#include "CD_NamespaceFooter.H"
