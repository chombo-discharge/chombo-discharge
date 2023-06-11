/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CopyStrategy.H
  @brief  Simple enum for distinguishing copying strategies
  @author Robert Marskar
*/

#ifndef CD_CopyStrategy_H
#define CD_CopyStrategy_H

#include <CD_NamespaceHeader.H>

/*!
  @brief Enum for distinguishing how we copy data
  Valid => valid region
  ValidGhost => valid+ghost region
*/
enum class CopyStrategy
{
  Valid,
  ValidGhost
};
#include <CD_NamespaceFooter.H>

#endif
