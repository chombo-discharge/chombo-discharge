/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzEBBC.cpp
  @brief  Implementation of CD_EBHelmholtzEBBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBCellFactory.H>

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

constexpr int EBHelmholtzEBBC::m_comp;
constexpr int EBHelmholtzEBBC::m_nComp;

EBHelmholtzEBBC::EBHelmholtzEBBC(){
}

EBHelmholtzEBBC::~EBHelmholtzEBBC(){

}

void EBHelmholtzEBBC::define(const EBLevelGrid&                                 a_eblg,
			     const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_Bcoef,
			     const RealVect&                                    a_probLo,
			     const Real&                                        a_dx,
			     const int                                          a_ghostCF,
			     const int                                          a_refRat){
  m_Bcoef   = a_Bcoef;
  m_eblg    = a_eblg;
  m_probLo  = a_probLo;
  m_dx      = a_dx;
  m_ghostCF = a_ghostCF;

  if(a_refRat > 1){
    m_hasFine = true;
    m_refRat  = a_refRat;

    refine(m_eblgFine, m_eblg, a_refRat);

    m_fineData.define(m_eblgFine.getDBL(), m_nComp, IntVect::Zero, EBCellFactory(m_eblgFine.getEBISL()));
  }
  
  this->define();
}

void EBHelmholtzEBBC::copyFineDataToBuffer(const LevelData<EBCellFAB>& a_fineData) const {
  a_fineData.copyTo(m_fineData);
}

const LayoutData<BaseIVFAB<VoFStencil> >& EBHelmholtzEBBC::getKappaDivFStencils() const {
  return m_kappaDivFStencils;
}

RealVect EBHelmholtzEBBC::getBoundaryPosition(const VolIndex& a_vof, const DataIndex& a_dit) const {
  const EBISBox&  ebisbox    = m_eblg.getEBISL()[a_dit];
  const RealVect& ebCentroid = ebisbox.bndryCentroid(a_vof);

  RealVect position = m_probLo + (0.5*RealVect::Unit + RealVect(a_vof.gridIndex()) + ebCentroid)*m_dx;

  return position;
}

#include <CD_NamespaceFooter.H>
