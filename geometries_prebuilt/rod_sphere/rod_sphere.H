/*!
  @file   rod_sphere.H
  @brief  Rod-sphere geometry (one dielectric sphere and one rod electrode)
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _ROD_SPHERE_GEOM_
#define _ROD_SPHERE_GEOM_

#include "computational_geometry.H"

#include "real_box.H"

/*!
  @brief This is a class that defines one electrode sphere and a dielectric sphere. Both can be turned off through class options.
*/
class rod_sphere : public computational_geometry {
public:

  rod_sphere();
  ~rod_sphere();
};

#endif
