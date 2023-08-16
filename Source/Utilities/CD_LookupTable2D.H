/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LookupTable1D.H
  @brief  Declaration of a table for looking up coefficients etc. 
  @author Robert Marskar
*/

#ifndef CD_LookupTable2D_H
#define CD_LookupTable2D_H

// Std includes
#include <iostream>
#include <vector>
#include <array>
#include <tuple>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_NamespaceHeader.H>
#include <CD_LookupTable.H>
#include <CD_LookupTable1D.H>

/*!
  @brief Class for interpolation of data f = (x,y) data.
*/
template <int N>
class LookupTable2D
{
public:
  static_assert(N >= 1, "LookupTable2D<N> must have N >= 1");

  /*!
    @brief Default constructor. Creates a table without any entires.
  */
  LookupTable2D() noexcept;

  /*!
    @brief Copy constructor
  */
  LookupTable2D(const LookupTable2D& a_otherTable) noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~LookupTable2D() noexcept;

  /*!
    @brief Interpolation function. 
    @details Returns all dependent variables using bilinear interpolation. Out-of-range strategies are respected.
    @param[in] a_x x-value
    @param[in] a_y y-value
  */
  inline std::array<Real, N>
  interpolate(const Real& a_x, const Real& a_y) const noexcept;

  /*!
    @brief Interpolation function for a specific variable. 
    @details Returns specified dependent variables using bilinear interpolation. Out-of-range strategies are respected.
    @param[in] a_x x-value
    @param[in] a_y y-value
  */
  template <size_t K>
  inline Real
  interpolate(const Real& a_x, const Real& a_y) const noexcept;

protected:
  /*!
    @brief Vector of one-dimensional tables. 
  */
  std::vector < LookupTable1D<N + 1>
                  /*!
    @brief Table spacing strategy along x
  */
                  TableSpacing m_tableSpacingX;

  /*!
    @brief Table spacing strategy along y
  */
  TableSpacing m_tableSpacingY;

  /*!
    @brief Out-of-range strategy along x
  */
  OutOfRangeStrategy m_strategyLoX;

  /*!
    @brief Out-of-range strategy along x
  */
  OutOfRangeStrategy m_strategyHiX;

  /*!
    @brief Out-of-range strategy along x
  */
  OutOfRangeStrategy m_strategyLoY;

  /*!
    @brief Out-of-range strategy along x
  */
  OutOfRangeStrategy m_strategyHiY;

  /*!
    @brief Spacing along first variable x
  */
  Real m_dx;

  /*!
    @brief Spacing along second variable y
  */
  Real m_dy;
};

#include <CD_NamespaceFooter.H>

//#include <CD_LookupTable1DImplem.H>

#endif
