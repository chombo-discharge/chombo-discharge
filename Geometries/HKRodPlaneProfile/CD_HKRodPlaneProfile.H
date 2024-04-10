/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_HKRodPlaneProfile.H
  @brief  Rod-plane geometry. The dielectric plane can have surface roughness
  @author Robert Marskar
*/

#ifndef CD_HKRodPlaneProfile_H
#define CD_HKRodPlaneProfile_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Geometry class which defines an electrode blade and a dielectric slab with various surface profiles. 
*/
class HKRodPlaneProfile : public ComputationalGeometry
{
public:
  enum class profile
  {
    circle,
    square,
    shroom,
  };

  HKRodPlaneProfile();
  ~HKRodPlaneProfile();

protected:
  profile m_profile;

  BaseIF*
  getBaseIF();
  BaseIF*
  getBaseIFCircle();
  BaseIF*
  getBaseIFSquare();
};

#include <CD_NamespaceFooter.H>

#endif
