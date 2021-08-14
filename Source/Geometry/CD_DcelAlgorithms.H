/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DcelAlgorithms.H
  @brief  Summary of some enum classes for navigating dcel functionality
  @author Robert Marskar
*/

#ifndef CD_DcelAlgorithms_H
#define CD_DcelAlgorithms_H

// Our includes
#include <CD_NamespaceHeader.H>

namespace Dcel {

  /*!
    @brief Supported algorithms for performing inside/outside tests when checking if a point projects to the inside or outside of a polygon face. 
  */
  enum class InsideOutsideAlgorithm {
    SubtendedAngle,
    CrossingNumber,
    WindingNumber
  };

  /*!
    @brief How to weight vertex normal
  */
  enum class VertexNormalWeight {
    None,
    Angle,
  };

  /*!
    @brief Possible search algorithms for BVH. 
    @details Direct means compute the signed distance for all primitives, Direct2 means compute the squared signed distance for all primitives,
    and KdTree means 'use a bounding volume hiearchy'.
  */
  enum class SearchAlgorithm{
    Direct,
    Direct2,
    KdTree
  };
}

#include <CD_NamespaceFooter.H>

#endif
