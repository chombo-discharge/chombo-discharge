/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzEBBCImplem.H
  @brief  Implementation of CD_EBHelmholtzEBBC.H
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzEBBCImplem_H
#define CD_EBHelmholtzEBBCImplem_H

// Chombo includes
#include <NeighborIterator.H>

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

inline
RealVect EBHelmholtzEBBC::getBoundaryPosition(const VolIndex& a_vof, const DataIndex& a_dit) const {

  // TLDR: Return the physical coordinates on the EB centroid.
  
  const EBISBox&  ebisbox    = m_eblg.getEBISL()[a_dit];
  const RealVect& ebCentroid = ebisbox.bndryCentroid(a_vof);

  RealVect position = m_probLo + (0.5*RealVect::Unit + RealVect(a_vof.gridIndex()) + ebCentroid)*m_dx;

  return position;
}

inline
Real EBHelmholtzEBBC::applyStencil(const VoFStencil& a_stencil, const EBCellFAB& a_phi) const {

  // TLDR: We just apply a stencil here. 
  Real ret = 0.0;
  for (int i = 0; i < a_stencil.size(); i++){
    ret += a_stencil.weight(i) * a_phi(a_stencil.vof(i), m_comp);
  }

  return ret;
}

inline
bool EBHelmholtzEBBC::isStencilValidCF(const VoFStencil& a_stencil, const DataIndex& a_dit) const {

  // TLDR: In this routine we check if the input stencil has a range that reaches into ghost cells
  //       that will not be filled by interpolator. We just get the neighboring boxes around the
  //       input grid patch with index a_dit. We then check that each vof in the input stencil reaches
  //       either into valid cells, or ghost cells up to a specified range.

  CH_assert(m_ghostCF >= 0);
  
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

#endif
