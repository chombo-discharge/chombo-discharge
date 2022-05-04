/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RodDielectric.H
  @brief  Definition of a electrode rod and dielectric dielectric geometry
  @author Robert Marskar
*/

#ifndef CD_RodDielectric_H
#define CD_RodDielectric_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that defines an electrode rod and a dielectric (some shapes supported)
*/
class RodDielectric : public ComputationalGeometry
{
public:
  RodDielectric();
  ~RodDielectric();

protected:
  void
  defineElectrode();
  void
  defineInsulator();

  RefCountedPtr<BaseIF>
  getBox();
  RefCountedPtr<BaseIF>
  getPlane();
  RefCountedPtr<BaseIF>
  getSphere();
  RefCountedPtr<BaseIF>
  getPerlinBox();
};

#include <CD_NamespaceFooter.H>

#endif
