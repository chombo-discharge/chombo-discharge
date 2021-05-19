/*!
  @file   rough_sphere.H
  @brief  Sphere geometry
  @author Robert Marskar
  @date   Oct. 2018
*/

#ifndef _ROUGH_SPHERE_GEOM_
#define _ROUGH_SPHERE_GEOM_

#include <CD_ComputationalGeometry.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief This is a class that a rough sphere (either electrode or dielectric)
*/
class rough_sphere : public ComputationalGeometry {
public:

  rough_sphere();
  ~rough_sphere();
};
#include "CD_NamespaceFooter.H"

#endif
