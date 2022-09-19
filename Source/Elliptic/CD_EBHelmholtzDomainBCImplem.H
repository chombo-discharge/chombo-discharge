/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDomainBCImplem.H
  @brief  Implementation of CD_EBHelmholtzDomainBC.H
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzDomainBCImplem_H
#define CD_EBHelmholtzDomainBCImplem_H

// Our includes
#include <CD_EBHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

inline RealVect
EBHelmholtzDomainBC::getBoundaryPosition(const IntVect& a_iv, const int& a_dir, const Side::LoHiSide& a_side) const
{

  const RealVect pos = m_probLo + (0.5 * RealVect::Unit + RealVect(a_iv)) * m_dx +
                       0.5 * m_dx * RealVect(BASISV(a_dir)) * sign(a_side);

  return pos;
}

#include <CD_NamespaceFooter.H>

#endif
