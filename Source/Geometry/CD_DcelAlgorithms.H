/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
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

#include <CD_NamespaceFooter.H>

#endif
