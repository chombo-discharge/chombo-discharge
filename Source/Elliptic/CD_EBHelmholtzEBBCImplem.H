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

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

inline
Real EBHelmholtzEBBC::applyStencil(const VoFStencil& a_stencil, const EBCellFAB& a_phi) const {
  Real ret = 0.0;
  for (int i = 0; i < a_stencil.size(); i++){
    ret += a_stencil.weight(i)*a_phi(a_stencil.vof(i), m_comp);
  }

  return ret;
}

#include <CD_NamespaceFooter.H>

#endif
