/*!
  @file sphere_sphere_geometry.H
  @brief Sphere-sphere geometry (one dielectric sphere and one electrode)
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _SPHERE_SPHERE_GEOM_
#define _SPHERE_SPHERE_GEOM_

#include "ComputationalGeometry.H"

/*!
  @brief This is a class that defines one electrode sphere and one dielectric sphere. 
*/
class sphere_sphere_geometry : public ComputationalGeometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  sphere_sphere_geometry();

  /*!
    @brief Destructor
  */
  ~sphere_sphere_geometry();

  /*!
    @brief Dump script for this class
  */
  static void dumpScript(std::string a_filename = "sphere_sphere_geometry.input");
};

#endif
