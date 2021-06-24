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
#include <NeighborIterator.H>

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

constexpr int EBHelmholtzEBBC::m_comp;
constexpr int EBHelmholtzEBBC::m_nComp;

EBHelmholtzEBBC::EBHelmholtzEBBC(){
  m_isMGLevel = false;
}

EBHelmholtzEBBC::~EBHelmholtzEBBC(){

}

void EBHelmholtzEBBC::setMG(const bool a_isMGLevel){
  m_isMGLevel = a_isMGLevel;
}

void EBHelmholtzEBBC::define(const EBLevelGrid&                                 a_eblg,
			     const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_Bcoef,
			     const RealVect&                                    a_probLo,
			     const Real&                                        a_dx,
			     const int                                          a_ghostCF){
  m_Bcoef   = a_Bcoef;
  m_eblg    = a_eblg;
  m_probLo  = a_probLo;
  m_dx      = a_dx;
  m_ghostCF = a_ghostCF;
  
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

bool EBHelmholtzEBBC::isStencilValidCF(const VoFStencil& a_stencil, const DataIndex& a_dit) const {
  const DisjointBoxLayout& dbl = m_eblg.getDBL();
  const ProblemDomain& domain  = m_eblg.getDomain();

  // Construct boxes that contain all the valid cells for this stencil. 
  std::vector<Box> validBoxes;
  
  Box curBox = dbl[a_dit];
  curBox.grow(m_ghostCF);
  curBox &= domain;
  validBoxes.emplace_back(curBox);

  NeighborIterator nit(dbl);
  for (nit.begin(a_dit); nit.ok(); ++nit){
    Box neighBox = dbl[nit()];
    neighBox.grow(m_ghostCF);
    neighBox &= domain;

    validBoxes.emplace_back(neighBox);
  }

  // Now check that the stencil if the stencil. We set valid = false
  // if any of the stencil points reaches out of the ghosted boxes. 
  bool valid = true;
  
  for (int i = 0; i < a_stencil.size(); i++){
    const VolIndex& vof = a_stencil.vof(i);
    const IntVect& iv   = vof.gridIndex();
    
    bool insideOneBox = false;
    for (const auto& b : validBoxes){
      if(b.contains(iv)) insideOneBox = true;
    }

    if(!insideOneBox) {
      valid = false;
      break;
    }
  }

  return valid;
}

#include <CD_NamespaceFooter.H>
