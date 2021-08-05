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

      // Make a lambda which allows us to pass in the function using EBHelmholtzDomainBC API, which takes a std::function<Real(const RealVect a_position)>
      // type of function.
      //
      // This is the function type that the EBHelmholtzOp API requires, and it is a design choice mandated by our choice to make that operator
      // time-independent. Although this might seem weird, the time dependence is nonetheless passed in because a_electrostaticBCs are passed in
      // from FieldSolver, and in that solver we capture FieldSolver::m_time by reference.
      // 
      // This might seem clunky, but I can't see any other way of doing it properly without changing EBHelmholtzOp to a time-dependent operator (which
      // I really don't want to do).
      //
      // Here is the function
      auto func = [&bcFunc](const RealVect& a_position) -> Real {
	constexpr Real dummyDt = 0.0;
	return bcFunc(a_position, dummyDt);
      };

      // Create either a Dirichlet or Neumann boundary condition object for the current domain edge (face in 3D).
      switch(bcType) {
      case ElectrostaticDomainBc::BcType::Dirichlet:
	m_bcObjects.emplace(curWall, std::make_shared<EBHelmholtzDirichletDomainBC>(func));
	break;
      case ElectrostaticDomainBc::BcType::Neumann:
	m_bcObjects.emplace(curWall, std::make_shared<EBHelmholtzNeumannDomainBC>(func));
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
