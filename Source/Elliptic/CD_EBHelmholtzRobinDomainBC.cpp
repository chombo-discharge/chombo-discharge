/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzRobinDomainBC.cpp
  @brief  Implementation of CD_EBHelmholtzRobinDomainBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzRobinDomainBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzRobinDomainBC::EBHelmholtzRobinDomainBC(){
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzRobinDomainBC::~EBHelmholtzRobinDomainBC(){

}


void EBHelmholtzRobinDomainBC::setCoefficients(const Real a_A, const Real a_B, const Real a_C){
  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
}


void EBHelmholtzRobinDomainBC::setCoefficients(const std::function<Real(const RealVect& a_pos) >& a_A,
					       const std::function<Real(const RealVect& a_pos) >& a_B,
					       const std::function<Real(const RealVect& a_pos) >& a_C){
  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;

}
  

void EBHelmholtzRobinDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
					   const BaseFab<Real>&  a_phi,
					   const int&            a_dir,
					   const Side::LoHiSide& a_side,
					   const DataIndex&      a_dit,
					   const bool            a_useHomogeneous) const {

  
}
  


Real EBHelmholtzRobinDomainBC::getFaceFlux(const VolIndex&       a_vof,
					   const EBCellFAB&      a_phi,
					   const int&            a_dir,
					   const Side::LoHiSide& a_side,
					   const DataIndex&      a_dit,
					   const bool            a_useHomogeneous) const {

}

#include <CD_NamespaceFooter.H>
