/*!
  @file   rod_plane_profile.H
  @brief  Rod-plane geometry. The dielectric plane can have surface roughness
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _ROD_PLANE_PROFILE_GEOM_
#define _ROD_PLANE_PROFILE_GEOM_

#include "computational_geometry.H"

#include "CD_NamespaceHeader.H"

/*!
  @brief Geometry class which defines an electrode blade and a dielectric slab with various surface profiles. 
*/
class rod_plane_profile : public computational_geometry {
public:

  enum class profile {
    circle,
    square,
  };

  rod_plane_profile();
  ~rod_plane_profile();

protected:

  profile m_profile;

  BaseIF* getBaseIF();
  BaseIF* getBaseIF_circle();
  BaseIF* getBaseIF_square();
};
#include "CD_NamespaceFooter.H"

#endif
