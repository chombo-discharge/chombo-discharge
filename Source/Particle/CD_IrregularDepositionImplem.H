/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_IrregularDepositionImplem.H
 * @brief  Implementation of CD_IrregularDeposition.H
 * @author Robert Marskar
 */

#ifndef CD_IRREGULARDEPOSITIONIMPLEM_H
#define CD_IRREGULARDEPOSITIONIMPLEM_H

// Our includes
#include <CD_IrregularDeposition.H>
#include <CD_NamespaceHeader.H>

inline IrregularDeposition
irregularDepositionFromString(const std::string& a_str) noexcept
{
  if (a_str == "native") {
    return IrregularDeposition::Native;
  }
  else if (a_str == "ngp") {
    return IrregularDeposition::NGP;
  }
  else if (a_str == "mirror") {
    return IrregularDeposition::Mirror;
  }
  else if (a_str == "redistribute") {
    return IrregularDeposition::Redistribute;
  }
  else if (a_str == "redistribute_blended") {
    return IrregularDeposition::RedistributeBlended;
  }
  else {
    MayDay::Abort(("irregularDepositionFromString - unknown cut-cell deposition '" + a_str + "'").c_str());

    return IrregularDeposition::Native;
  }
}

#include <CD_NamespaceFooter.H>

#endif
