/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_MechanicalShaft.H
  @brief  Declaration of the mechanical shaft geometry
  @author Robert Marskar
*/

#ifndef CD_MechanicalShaft_H
#define CD_MechanicalShaft_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Mechanical shaft class
*/
class MechanicalShaft : public ComputationalGeometry {
public:

  MechanicalShaft();
  ~MechanicalShaft();

protected:

  void defineElectrode();
  void defineDielectric();
  
  RefCountedPtr<BaseIF> getPolygon();
  RefCountedPtr<BaseIF> getCylinder();
  RefCountedPtr<BaseIF> getCylinderProfile();
};

#include <CD_NamespaceFooter.H>

#endif
