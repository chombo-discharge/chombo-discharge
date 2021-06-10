/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzEBBC.cpp
  @brief  Implementation of CD_EBHelmholtzEBBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

constexpr int EBHelmholtzEBBC::m_comp;
constexpr int EBHelmholtzEBBC::m_nComp;

EBHelmholtzEBBC::EBHelmholtzEBBC(){
}

EBHelmholtzEBBC::~EBHelmholtzEBBC(){

}

void EBHelmholtzEBBC::applyEBFlux(VoFIterator&       a_vofit,
				  EBCellFAB&         a_Lphi,
				  const EBCellFAB&   a_phi,
				  const DataIndex&   a_dit,
				  const Real&        a_beta) const {
  for (a_vofit.reset(); a_vofit.ok(); ++a_vofit){
    this->applyEBFlux(a_Lphi, a_phi, a_vofit(), a_dit, a_beta);
  }
}

void EBHelmholtzEBBC::define(const EBLevelGrid& a_eblg, const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_Bcoef, const RealVect& a_probLo, const Real& a_dx){
  m_Bcoef  = a_Bcoef;
  m_eblg   = a_eblg;
  m_probLo = a_probLo;
  m_dx     = a_dx;
  
  this->define();
}

const LayoutData<BaseIVFAB<VoFStencil> >& EBHelmholtzEBBC::getKappaDivFStencils() const{
  return m_kappaDivFStencils;
}

RealVect EBHelmholtzEBBC::getBoundaryPosition(const VolIndex& a_vof, const DataIndex& a_dit) const {
  const EBISBox&  ebisbox    = m_eblg.getEBISL()[a_dit];
  const RealVect& ebCentroid = ebisbox.bndryCentroid(a_vof);

  RealVect position = m_probLo + (0.5*RealVect::Unit + RealVect(a_vof.gridIndex()) + ebCentroid)*m_dx;

  return position;
}

#include <CD_NamespaceFooter.H>
