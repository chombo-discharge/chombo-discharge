/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_CoaxialCable.H
  @brief  Coaxial cable geometry
  @author Robert Marskar
*/

#ifndef CD_CoaxialCable_H
#define CD_CoaxialCable_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief This is a class that defines a coaxial cable geometry 
*/
class CoaxialCable : public ComputationalGeometry {
public:
  CoaxialCable();
  ~CoaxialCable();
};

#include <CD_NamespaceFooter.H>

#endif
