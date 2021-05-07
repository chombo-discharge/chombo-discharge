/*!
  @file   dcel_algorithms.H
  @brief  Summary of some enum classes for navigating dcel functionality
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_ALGORITHMS_H_
#define _DCEL_ALGORITHMS_H_

#include "CD_NamespaceHeader.H"

namespace dcel {

  enum class InsideOutsideAlgorithm {
    SubtendedAngle,
    CrossingNumber,
    WindingNumber
  };

  enum class VertexNormalWeight {
    None,
    Angle,
  };

  enum class SearchAlgorithm{
    Direct,
    Direct2,
    KdTree
  };
}
#include "CD_NamespaceFooter.H"

#endif
