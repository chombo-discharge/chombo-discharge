/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_JumpBCImplem.H
  @brief  Implementation of CD_JumpBC.H
  @author Robert Marskar
*/

#ifndef CD_JumpBCImplem_H
#define CD_JumpBCImplem_H

// Chombo includes
#include <NeighborIterator.H>

// Our includes
#include <CD_JumpBC.H>
#include <CD_NamespaceHeader.H>


inline
bool JumpBC::isStencilValidCF(const VoFStencil& a_stencil, const DataIndex& a_dit) const {
  const DisjointBoxLayout& dbl = m_mflg.getGrids();
  const ProblemDomain& domain  = m_mflg.getDomain();

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
