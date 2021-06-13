/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletDomainBC.cpp
  @brief  Implementation of CD_EBHelmholtzDirichletDomainBC.cpp
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzDirichletDomainBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzDirichletDomainBC::EBHelmholtzDirichletDomainBC(){
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzDirichletDomainBC::EBHelmholtzDirichletDomainBC(const Real a_value){
  this->setValue(a_value);
}

EBHelmholtzDirichletDomainBC::EBHelmholtzDirichletDomainBC(const std::function<Real(const RealVect& a_pos)>& a_value){
  this->setValue(a_value);
}

EBHelmholtzDirichletDomainBC::~EBHelmholtzDirichletDomainBC(){

}

void EBHelmholtzDirichletDomainBC::setValue(const Real a_value){
  m_useConstant = true;
  m_useFunction = false;
  
  m_constantValue = a_value;
}

void EBHelmholtzDirichletDomainBC::setValue(const std::function<Real(const RealVect& a_pos)>& a_value){
  m_useConstant = false;
  m_useFunction = true;
  
  m_functionValue = a_value;
}

void EBHelmholtzDirichletDomainBC::getFaceFlux(FArrayBox& a_faceFlux, const FArrayBox& a_phi, const int& dir, const Side::LoHiSide& a_side) const {
  
}

void EBHelmholtzDirichletDomainBC::getFaceFlux(Real&                 a_faceFlux,
					       const VolIndex&       a_vof,
					       const EBCellFAB&      a_phi,
					       const int&            a_dir,
					       const Side::LoHiSide& a_side,
					       const DataIndex&      a_dit) const {

}

#include <CD_NamespaceFooter.H>
