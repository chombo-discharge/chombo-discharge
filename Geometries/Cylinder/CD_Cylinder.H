/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Cylinder.H
  @brief  Cylinder geometry
  @author Robert Marskar
*/

#ifndef CD_Cylinder_H
#define CD_Cylinder_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief This is a class that defines a cylinder geometry. 
*/
class Cylinder : public ComputationalGeometry
{
public:
  Cylinder();
  ~Cylinder();
};

#include <CD_NamespaceFooter.H>

#endif
