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

/*!
  @brief Class for interpolation of data f = (x,y) data
*/
template <int Nx, int Ny>
class LookupTable2D
{
public:
  static_assert(Nx > 1, "LookupTable2D<Nx,Ny> must have Nx > 1");
  static_assert(Ny > 1, "LookupTable2D<Nx,Ny> must have Ny > 1");

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
  ~LookupTable2D() noexcept;

protected:
};

#include <CD_NamespaceFooter.H>

//#include <CD_LookupTable1DImplem.H>

#endif
