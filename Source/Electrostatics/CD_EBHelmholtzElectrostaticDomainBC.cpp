/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzElectrostaticDomainBC.cpp
  @brief  Implementation of CD_EBHelmholtzElectrostaticDomainBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzElectrostaticDomainBC.H>
#include <CD_EBHelmholtzDirichletDomainBC.H>
#include <CD_EBHelmholtzNeumannDomainBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzElectrostaticDomainBC::EBHelmholtzElectrostaticDomainBC(const ElectrostaticDomainBc& a_electrostaticBCs){
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const ElectrostaticDomainBc::Wall        curWall = std::make_pair(dir, sit());
      const ElectrostaticDomainBc::BcType&     bcType  = a_electrostaticBCs.getBc(curWall).first;
      const ElectrostaticDomainBc::BcFunction& bcFunc  = a_electrostaticBCs.getBc(curWall).second;

      // Generate the required BC type
      switch(bcType) {
      case ElectrostaticDomainBc::BcType::Dirichlet:
	m_bcObjects.emplace(curWall, std::make_shared<EBHelmholtzDirichletDomainBC>(1.0));
	break;
      case ElectrostaticDomainBc::BcType::Neumann:
	m_bcObjects.emplace(curWall, std::make_shared<EBHelmholtzNeumannDomainBC>(0.0));
	break;
      default:
	MayDay::Error("EBHelmholtzElectrostaticDomainBC::EBHelmholtzElectrostaticDomainBC - unsupported boundary condition passed into constructor!");
      }
    }
  }
}

EBHelmholtzElectrostaticDomainBC::~EBHelmholtzElectrostaticDomainBC(){

}

void EBHelmholtzElectrostaticDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
						   const BaseFab<Real>&  a_phi,
						   const int&            a_dir,
						   const Side::LoHiSide& a_side,
						   const DataIndex&      a_dit,
						   const bool            a_useHomogeneous) const {

  const auto& bcPtr = m_bcObjects.at(std::make_pair(a_dir, a_side));

  bcPtr->getFaceFlux(a_faceFlux,
		     a_phi,
		     a_dir,
		     a_side,
		     a_dit,
		     a_useHomogeneous);
}

Real EBHelmholtzElectrostaticDomainBC::getFaceFlux(const VolIndex&       a_vof,
						   const EBCellFAB&      a_phi,
						   const int&            a_dir,
						   const Side::LoHiSide& a_side,
						   const DataIndex&      a_dit,
						   const bool            a_useHomogeneous) const {

  const auto& bcPtr = m_bcObjects.at(std::make_pair(a_dir, a_side));

  return bcPtr->getFaceFlux(a_vof, a_phi, a_dir, a_side, a_dit, a_useHomogeneous);
}

#include <CD_NamespaceFooter.H>
