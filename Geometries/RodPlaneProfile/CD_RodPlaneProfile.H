/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RodPlaneProfile.H
  @brief  Rod-plane geometry. The dielectric plane can have surface roughness
  @author Robert Marskar
*/

#ifndef CD_RodPlaneProfile_H
#define CD_RodPlaneProfile_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Geometry class which defines an electrode blade and a dielectric slab with various surface profiles. 
*/
class RodPlaneProfile : public ComputationalGeometry {
public:

  enum class profile {
    circle,
    square,
  };

  RodPlaneProfile();
  ~RodPlaneProfile();

protected:

  profile m_profile;

  BaseIF* getBaseIF();
  BaseIF* getBaseIFCircle();
  BaseIF* getBaseIFSquare();
};

#include <CD_NamespaceFooter.H>

#endif
