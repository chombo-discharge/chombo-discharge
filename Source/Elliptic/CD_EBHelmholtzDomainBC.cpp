/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDomainBC.cpp
  @brief  Implementation of CD_EBHelmholtzDomainBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

constexpr int EBHelmholtzDomainBC::m_comp;
constexpr int EBHelmholtzDomainBC::m_nComp;

EBHelmholtzDomainBC::EBHelmholtzDomainBC(){

}

EBHelmholtzDomainBC::~EBHelmholtzDomainBC(){

}

void EBHelmholtzDomainBC::define(const EBLevelGrid& a_eblg, const RefCountedPtr<LevelData<EBFluxFAB> >& a_Bcoef, const RealVect& a_probLo, const Real a_dx){
  m_eblg    = a_eblg;
  m_Bcoef   = a_Bcoef;
  m_probLo  = a_probLo;
  m_dx      = a_dx;
}

RealVect EBHelmholtzDomainBC::getBoundaryPosition(const VolIndex& a_vof, const int& a_dir, const Side::LoHiSide& a_side) const {

  RealVect pos = m_probLo + (0.5*RealVect::Unit + RealVect(a_vof.gridIndex()))*m_dx + 0.5*m_dx*RealVect(BASISV(a_dir))*sign(a_side);

  return pos;
}

#include <CD_NamespaceFooter.H>
