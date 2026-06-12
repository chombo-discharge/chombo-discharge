/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_PolyUtils.H
  @brief  Agglomeration of some useful algebraic/polynomial routines
  @author Robert Marskar
*/

#ifndef CD_POLYUTILS_H
#define CD_POLYUTILS_H

// Chombo includes
#include <RefCountedPtr.H>
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Agglomeration class for various useful polynomial routines
*/
namespace PolyUtils {

  /**
    @brief Compute the root of a function between two points. This is a 1D problem along the line
    @param[in] a_impFunc Implicit function
    @param[in] a_point1  First end-point
    @param[in] a_point2  Other end-point
    @note The root MUST be bracketed, i.e. f(x1)*f(x2) < 0.0. Anything else is an error.
    @return Return value
  */
  RealVect
  brentRootFinder(const RefCountedPtr<BaseIF>& a_impFunc, const RealVect& a_point1, const RealVect& a_point2);

  /**
    @brief Compute the root of a function between two points. This is a 1D problem.
    @return Returns the root (if it could be found). Returns std::numeric_limits<Real>::max() otherwise.
    @param[in] a_point1 Point1
    @param[in] a_point2 Point2
    @param[in] a_func Func
  */
  Real
  brentSolve(const Real a_point1, const Real a_point2, const std::function<Real(const Real x)>& a_func) noexcept;
} // namespace PolyUtils

#include <CD_NamespaceFooter.H>

#endif
