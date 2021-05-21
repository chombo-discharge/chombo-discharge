/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_RoughRod.H
  @brief  Rod geometry with surface roughness at tip
  @author Robert Marskar
*/

#ifndef CD_RoughRod_H
#define CD_RoughRod_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that describes a rod electrode with surface roughness
*/
class RoughRod : public ComputationalGeometry {
public:

  RoughRod();
  ~RoughRod();
};

#include <CD_NamespaceFooter.H>

#endif
