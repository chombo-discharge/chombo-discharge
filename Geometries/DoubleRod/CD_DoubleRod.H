/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_DoubleRod.H
  @brief  Double rod geometry.
  @author Robert Marskar
*/

#ifndef CD_DoubleRod_H
#define CD_DoubleRod_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief This is a class that defines two dielectric rods.
*/
class DoubleRod : public ComputationalGeometry {
public:

  DoubleRod();
  ~DoubleRod();
};

#include <CD_NamespaceFooter.H>

#endif
