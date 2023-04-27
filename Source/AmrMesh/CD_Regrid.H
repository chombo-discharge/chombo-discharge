/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Regrid.H
  @brief  Declaration of averaging methods.
  @author Robert Marskar
*/

#ifndef CD_Regrid_H
#define CD_Regrid_H

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Various regridding methods.
*/
enum class Regrid
{
  Arithmetic,
  Conservative,
  MinMod
};

#include <CD_NamespaceFooter.H>

#endif
