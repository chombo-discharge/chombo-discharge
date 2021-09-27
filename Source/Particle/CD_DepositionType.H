/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DepositionType.H
  @brief  Declaration of deposition types. 
  @author Robert Marskar
*/

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Deposition types. 
*/
enum class DepositionType {
  NGP,
  CIC,
  TSC,
  W4
};

#include <CD_NamespaceFooter.H>
