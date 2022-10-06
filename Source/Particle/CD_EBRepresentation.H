/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBRepresentation.H
  @brief  Enum for putting some logic into EB intersection tests. 
  @author Robert Marskar
*/

#ifndef CD_EBRepresentation_H
#define CD_EBRepresentation_H

#include <CD_NamespaceHeader.H>

/*!
  @brief Enum for putting some logic into how we think about EBs. This is just a simply supporting class for various solvers/algorithms that can distinguish
  between ways of think about the EB: 
*/
enum class EBRepresentation
{
  ImplicitFunction,
  Discrete,
  Voxel,
  Levelset
};

#include <CD_NamespaceFooter.H>

#endif
