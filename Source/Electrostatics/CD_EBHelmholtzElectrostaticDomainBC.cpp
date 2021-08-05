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
#include <CD_NamespaceHeader.H>

EBHelmholtzElectrostaticDomainBC::EBHelmholtzElectrostaticDomainBC(const ElectrostaticDomainBc& a_electrostaticBCs){
  m_electrostaticBCs = a_electrostaticBCs;
}

EBHelmholtzElectrostaticDomainBC::~EBHelmholtzElectrostaticDomainBC(){

}

void EBHelmholtzElectrostaticDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
						   const BaseFab<Real>&  a_phi,
						   const int&            a_dir,
						   const Side::LoHiSide& a_side,
						   const DataIndex&      a_dit,
						   const bool            a_useHomogeneous) const {

}

Real EBHelmholtzElectrostaticDomainBC::getFaceFlux(const VolIndex&       a_vof,
						   const EBCellFAB&      a_phi,
						   const int&            a_dir,
						   const Side::LoHiSide& a_side,
						   const DataIndex&      a_dit,
						   const bool            a_useHomogeneous) const {

}

#include <CD_NamespaceFooter.H>
