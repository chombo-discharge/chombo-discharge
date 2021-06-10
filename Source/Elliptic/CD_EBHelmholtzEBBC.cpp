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

EBHelmholtzEBBC::EBHelmholtzEBBC(){
  m_comp  = 0;
  m_nComp = 1;
}

EBHelmholtzEBBC::~EBHelmholtzEBBC(){

}

void EBHelmholtzEBBC::applyEBFlux(VoFIterator&       a_vofit,
				  EBCellFAB&         a_Lphi,
				  const EBCellFAB&   a_phi,
				  const Real&        a_beta) const {
  for (a_vofit.reset(); a_vofit.ok(); ++a_vofit){
    this->applyEBFlux(a_Lphi, a_phi, a_vofit(), a_beta);
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

RealVect EBHelmholtzEBBC::position(VolIndex& a_vof, const EBCellFAB& a_phi) const {
  const EBISBox&  ebisbox    = a_phi.getEBISBox();
  const RealVect& ebCentroid = ebisbox.bndryCentroid(a_vof);

  RealVect position = m_probLo + (0.5*RealVect::Unit + a_vof.gridIndex())*m_dx + ebCentroid*m_dx;
}

#include <CD_NamespaceFooter.H>
