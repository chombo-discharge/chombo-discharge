/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DepositionType.H
  @brief  Declaration of deposition types. 
  @author Robert Marskar
*/

#ifndef CD_DepositionType_H
#define CD_DepositionType_H

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Deposition types. 
  @details NGP = Nearest grid point, CIC = Cloud-in-cell, TSC = Triangle-shaped cloud
*/
enum class DepositionType
{
  NGP,
  CIC,
  TSC
};

#include <CD_NamespaceFooter.H>

#endif
