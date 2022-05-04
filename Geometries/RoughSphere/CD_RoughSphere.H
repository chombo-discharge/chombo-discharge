/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RoughSphere.H
  @brief  Rough sphere geometry
  @author Robert Marskar
*/

#ifndef CD_RoughSphere_H
#define CD_RoughSphere_H

#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief This is a class that a rough sphere (either electrode or dielectric)
*/
class RoughSphere : public ComputationalGeometry
{
public:
  RoughSphere();
  ~RoughSphere();
};

#include <CD_NamespaceFooter.H>

#endif
