/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ScanShopImplem.H
 * @brief  Implementation of CD_ScanShop.H
 * @author Robert Marskar
 */

#ifndef CD_SCANSHOPIMPLEM_H
#define CD_SCANSHOPIMPLEM_H

// Std includes
#include <algorithm>

// Our includes
#include <CD_ScanShop.H>
#include <CD_NamespaceHeader.H>

inline int
ScanShop::getScanSkip(Real a_slack, Real a_dx) const
{
  // A signed distance function changes by at most the distance travelled, so no cell centre closer than
  // a_slack to the evaluated point can cross the threshold. Those cells are passed over. The guard keeps
  // the skipped distance strictly below a_slack, so a cell exactly a_slack away is still evaluated.
  if (!m_distanceSkip || a_slack <= 0.0) {
    return 0;
  }

  int skip = static_cast<int>(a_slack / a_dx);

  if (skip * a_dx >= a_slack) {
    skip--;
  }

  return std::max(skip, 0);
}

inline bool
ScanShop::isRegular(Box a_box, const RealVect& a_probLo, Real a_dx) const
{
  CH_TIME("ScanShop::isRegular(Box, RealVect, Real)");

  // Short-circuit search over the box: a direct index loop with early return rather than a BoxLoops::loop,
  // since BoxLoops::loop cannot break early and m_baseIF->value(...) is expensive.
  const Real halfDiagonal = 0.5 * a_dx * sqrt(SpaceDim);

#if CH_SPACEDIM == 3
  for (int k = a_box.smallEnd(2); k <= a_box.bigEnd(2); k++) {
#endif
    for (int j = a_box.smallEnd(1); j <= a_box.bigEnd(1); j++) {
      for (int i = a_box.smallEnd(0); i <= a_box.bigEnd(0); i++) {
        const RealVect point = a_probLo + a_dx * (0.5 * RealVect::Unit + RealVect(IntVect(D_DECL(i, j, k))));
        const Real     value = m_baseIF->value(point);

        if (value >= -halfDiagonal) {
          return false;
        }

        i += this->getScanSkip(-value - halfDiagonal, a_dx);
      }
    }
#if CH_SPACEDIM == 3
  }
#endif

  return true;
}

inline bool
ScanShop::isCovered(Box a_box, const RealVect& a_probLo, Real a_dx) const
{
  CH_TIME("ScanShop::isCovered(Box, RealVect, Real)");

  // Short-circuit search over the box: a direct index loop with early return rather than a BoxLoops::loop,
  // since BoxLoops::loop cannot break early and m_baseIF->value(...) is expensive.
  const Real halfDiagonal = 0.5 * a_dx * sqrt(SpaceDim);

#if CH_SPACEDIM == 3
  for (int k = a_box.smallEnd(2); k <= a_box.bigEnd(2); k++) {
#endif
    for (int j = a_box.smallEnd(1); j <= a_box.bigEnd(1); j++) {
      for (int i = a_box.smallEnd(0); i <= a_box.bigEnd(0); i++) {
        const RealVect point = a_probLo + a_dx * (0.5 * RealVect::Unit + RealVect(IntVect(D_DECL(i, j, k))));
        const Real     value = m_baseIF->value(point);

        if (value <= halfDiagonal) {
          return false;
        }

        i += this->getScanSkip(value - halfDiagonal, a_dx);
      }
    }
#if CH_SPACEDIM == 3
  }
#endif

  return true;
}

inline std::vector<std::pair<Box, int>>
ScanShop::getSortedBoxesAndTypes(const Vector<Box>& a_boxes, const Vector<int>& a_types) const
{

  std::vector<std::pair<Box, int>> sortedBoxesAndTypes;

  for (int i = 0; i < a_boxes.size(); i++) {
    sortedBoxesAndTypes.emplace_back(std::make_pair(a_boxes[i], a_types[i]));
  }

  auto comparator = [](const std::pair<Box, int>& a, const std::pair<Box, int>& b) -> bool {
    return a.first < b.first;
  };

  std::sort(sortedBoxesAndTypes.begin(), sortedBoxesAndTypes.end(), comparator);

  return sortedBoxesAndTypes;
}

#endif

#include <CD_NamespaceFooter.H>
