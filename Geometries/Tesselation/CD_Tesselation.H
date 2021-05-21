/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_Tesselation.H
  @brief  Geometry for surface tesselation
  @author Robert Marskar
*/

#ifndef CD_Tesselation_H
#define CD_Tesselation_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

class Tesselation : public ComputationalGeometry {
public:

  Tesselation();
  ~Tesselation();
};

#include <CD_NamespaceFooter.H>

#endif
