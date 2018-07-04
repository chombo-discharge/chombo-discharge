/*!
  @file   coaxial_cable.H
  @brief  Coaxial cable geometry
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _COAXIAL_CABLE_
#define _COAXIAL_CABLE_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines a coaxial cable geometry 
*/
class coaxial_cable : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  coaxial_cable();

  /*!
    @brief Destructor
  */
  ~coaxial_cable();
};

#endif
