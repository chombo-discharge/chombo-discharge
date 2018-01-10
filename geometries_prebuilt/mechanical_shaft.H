/*!
  @file mechanical_shaft_geometry.H
  @brief Sphere-sphere geometry (one dielectric sphere and one electrode)
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _MECHANICAL_SHAFT_
#define _MECHANICAL_SHAFT_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines one electrode sphere and one dielectric sphere. 
*/
class mechanical_shaft : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  mechanical_shaft();

  /*!
    @brief Destructor
  */
  ~mechanical_shaft();
};

#endif
