/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  @file   CD_EBHelmholtzDomainBCImplem.H
  @brief  Implementation of CD_EBHelmholtzDomainBC.H
  @author Robert Marskar
*/

#ifndef CD_EBHELMHOLTZDOMAINBCIMPLEM_H
#define CD_EBHELMHOLTZDOMAINBCIMPLEM_H

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
