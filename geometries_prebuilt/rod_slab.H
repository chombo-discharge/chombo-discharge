/*!
  @file rod_slab_geometry.H
  @brief Sphere-sphere geometry (one dielectric sphere and one electrode)
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _ROD_SLAB_GEOM_
#define _ROD_SLAB_GEOM_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines one electrode sphere and one dielectric sphere. 
*/
class rod_slab : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  rod_slab();

  /*!
    @brief Destructor
  */
  ~rod_slab();
};

#endif
