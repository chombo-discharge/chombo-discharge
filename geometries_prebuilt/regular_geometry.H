/*!
  @file   regular_geometry.H
  @brief  Regular geometry (no EBs)
  @author Robert Marskar
  @date   July 2018
*/

#ifndef _REGULAR_GEOMETRY_GEOM_
#define _REGULAR_GEOMETRY_GEOM_

#include "computational_geometry.H"

/*!
  @brief This class implements an empty geometry
*/
class regular_geometry : public computational_geometry {
public:

  /*!
    @brief Constructor. First sphere is the electrode. 
  */
  regular_geometry();

  /*!
    @brief Destructor
  */
  ~regular_geometry();
};

#endif
