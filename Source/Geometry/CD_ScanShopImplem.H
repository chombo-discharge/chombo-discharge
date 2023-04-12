/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @brief  CD_ScanShopImplem.H
  @brief  Implementation of CD_ScanShop.H
  @author Robert Marskar
*/

#ifndef CD_ScanShopImplem_H
#define CD_ScanShopImplem_H

// Our includes
#include <CD_ScanShop.H>
#include <CD_NamespaceHeader.H>

inline bool
ScanShop::isRegular(const Box a_box, const RealVect a_probLo, const Real a_dx) const
{
  CH_TIME("ScanShop::isRegular(Box, RealVect, Real)");

  bool ret = true;

  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const RealVect a_point = a_probLo + a_dx * (0.5 * RealVect::Unit + RealVect(bit()));
    if (m_baseIF->value(a_point) >= -0.5 * a_dx * sqrt(SpaceDim)) {
      ret = false;

      break;
    }
  }

  return ret;
}

inline bool
ScanShop::isCovered(const Box a_box, const RealVect a_probLo, const Real a_dx) const
{
  CH_TIME("ScanShop::isCovered(Box, RealVect, Real)");

  bool ret = true;

  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const RealVect a_point = a_probLo + a_dx * (0.5 * RealVect::Unit + RealVect(bit()));
    if (m_baseIF->value(a_point) <= 0.5 * a_dx * sqrt(SpaceDim)) {
      ret = false;

      break;
    }
  }

  return ret;
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
