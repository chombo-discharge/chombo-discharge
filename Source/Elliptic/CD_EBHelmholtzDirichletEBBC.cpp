/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletEBBC.cpp
  @brief  Implementation of CD_EBHelmholtzDiricheltEBBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_LeastSquares.H>
#include <CD_EBHelmholtzDirichletEBBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzDirichletEBBC::EBHelmholtzDirichletEBBC(){
  m_order       = -1;
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzDirichletEBBC::~EBHelmholtzDirichletEBBC(){

}

void EBHelmholtzDirichletEBBC::setOrder(const int a_order){
  m_order = a_order;
}

void EBHelmholtzDirichletEBBC::setValue(const int a_value){
  m_useConstant = true;
  m_useFunction = false;
  
  m_constantValue = a_value;
}

void EBHelmholtzDirichletEBBC::setValue(const std::function<Real(const RealVect& a_pos)>& a_value){
  m_useConstant = false;
  m_useFunction = true;
  
  m_functionValue = a_value;
}
  
void EBHelmholtzDirichletEBBC::define() {
  if(!(m_order == 1  || m_order == 2 )) MayDay::Error("EBHelmholtzDirichletEBBC - order is not 1 or 2!");
  if(!(m_useConstant || m_useFunction)) MayDay::Error("EBHelmholtzDirichletEBBC - not using constant or function!");
}

void EBHelmholtzDirichletEBBC::applyEBFlux(EBCellFAB&         a_Lphi,
					   const EBCellFAB&   a_phi,
					   const VolIndex&    a_vof,
					   const Real&        a_factor,
					   const bool&        a_useHomogeneous) const {
  return;
}

#include <CD_NamespaceFooter.H>
