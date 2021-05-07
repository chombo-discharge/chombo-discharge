/*!
  @file   poly.H
  @brief  Agglomeration of some useful algebraic/polynomial routines
  @author Robert Marskar
*/

#ifndef _POLY_
#define _POLY_

#include <RefCountedPtr.H>
#include <BaseIF.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Agglomeration class for various useful polynomial routines
*/
namespace poly {
  RealVect brent_root_finder(const RefCountedPtr<BaseIF>& a_impFunc, const RealVect& a_point1, const RealVect& a_point2);
};

#include "CD_NamespaceFooter.H"

#endif
