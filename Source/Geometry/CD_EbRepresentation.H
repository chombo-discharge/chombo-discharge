/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_EbRepresentation.H
  @brief  Enum for putting some logic into EB
  @author Robert Marskar
*/

#ifndef CD_EbRepresentation_H
#define CD_EbRepresentation_H

#include <CD_NamespaceHeader.H>

/*!
  @brief Enum for putting some logic into how we thought about EBs
*/
enum class EbRepresentation {
  ImplicitFunction,
  Discrete,
  Voxel,
  Levelset
};

#include <CD_NamespaceFooter.H>

#endif
