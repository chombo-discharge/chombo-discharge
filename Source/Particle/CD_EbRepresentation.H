/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbRepresentation.H
  @brief  Enum for putting some logic into EB intersection tests. 
  @author Robert Marskar
*/

#ifndef CD_EbRepresentation_H
#define CD_EbRepresentation_H

#include <CD_NamespaceHeader.H>

/*!
  @brief Enum for putting some logic into how we think about EBs. This is just a simply supporting class for various solvers/algorithms that can distinguish
  between ways of think about the EB: 
*/
enum class EbRepresentation
{
  ImplicitFunction,
  Discrete,
  Voxel,
  Levelset
};

#include <CD_NamespaceFooter.H>

#endif
