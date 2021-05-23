/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PolyUtils.H
  @brief  Agglomeration of some useful algebraic/polynomial routines
  @author Robert Marskar
*/

#ifndef CD_PolyUtils_H
#define CD_PolyUtils_H

// Chombo includes
#include <RefCountedPtr.H>
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Agglomeration class for various useful polynomial routines
*/
namespace PolyUtils {
  RealVect brentRootFinder(const RefCountedPtr<BaseIF>& a_impFunc, const RealVect& a_point1, const RealVect& a_point2);
};

#include <CD_NamespaceFooter.H>

#endif
