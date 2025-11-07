/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BoxSorting.H
  @brief  Declaration of an enum for identifying various ways of sorting grid patches. 
  @author Robert Marskar
*/

#ifndef CD_BoxSorting_H
#define CD_BoxSorting_H

#include <CD_NamespaceHeader.H>

/*!
  @brief Enum for sorting boxes
*/
enum class BoxSorting
{
  None,
  Std,
  Shuffle,
  Morton,
  Hilbert
};

#include <CD_NamespaceFooter.H>

#endif
