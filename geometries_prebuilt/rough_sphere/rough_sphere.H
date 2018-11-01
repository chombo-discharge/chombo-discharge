/*!
  @file   rough_sphere.H
  @brief  Sphere geometry
  @author Robert Marskar
  @date   Oct. 2018
*/

#ifndef _ROUGH_SPHERE_GEOM_
#define _ROUGH_SPHERE_GEOM_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines a wire-wire geometry
*/
class rough_sphere : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  rough_sphere();

  /*!
    @brief Destructor
  */
  ~rough_sphere();
};

#endif
