/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Phases.H
  @brief  Declaration of phases
  @author Robert Marskar
*/

#ifndef CD_Phases_H
#define CD_Phases_H

#include <CD_NamespaceHeader.H>

/*!
  @brief Enum for holding phase identifiers. 
*/
enum class Phase {
  Gas,
  Solid
};

#include <CD_NamespaceFooter.H>

#endif
