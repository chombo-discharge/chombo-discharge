/*!
  @file   dcel_algorithms.H
  @brief  Summary of some enum classes for navigating dcel functionality
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_ALGORITHMS_H_
#define _DCEL_ALGORITHMS_H_

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

#endif
