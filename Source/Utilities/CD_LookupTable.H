/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LookupTable.H
  @brief  Declaration of a table for looking up coefficients etc. 
  @author Robert Marskar
*/

#ifndef CD_LookupTable_H
#define CD_LookupTable_H

// Std includes
#include <memory>
#include <vector>
#include <type_traits>
#include <iostream>

namespace LookupTable {
  /*!
  @brief Strategy for obtaining data when requesting data that is out of range. 
  We can use either constant (i.e. virtually padding the data with the first/last entry) or 
  with linear interpolation (i.e. virtually padding the data with a slope connecting the first
  two or last two entries in the table). 
*/
  enum class OutOfRangeStrategy
  {
    Constant,
    Interpolate
  };

  /*!
  @brief Enum for classifying the coordinate system.
*/
  enum class Spacing
  {
    Uniform,
    Exponential
  };

} // namespace LookupTable

#include <CD_LookupTable1D.H>

#endif
