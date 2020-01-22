/*!
  @file   fast_rod_plane.H
  @brief  Rod-plane geometry. The dielectric plane can have surface roughness
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _FAST_ROD_PLANE_GEOM_
#define _FAST_ROD_PLANE_GEOM_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines one electrode rod and a dielectric plane along one of the cardinal directions
*/
class fast_rod_plane : public computational_geometry {
public:

  fast_rod_plane();
  ~fast_rod_plane();

protected:


};

#endif
