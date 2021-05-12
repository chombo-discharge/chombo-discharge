/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_BoxSorting.H
  @brief  Declaration of various load-balancing routines
  @author Robert Marskar
*/

#ifndef CD_BoxSorting_H
#define CD_BoxSorting_H

#include <CD_NamespaceHeader.H>

/*!
  @brief Enum for sorting boxes
*/
enum class BoxSorting {
  None,
  Std,
  Shuffle,
  Morton
};

#include <CD_NamespaceFooter.H>

#endif
