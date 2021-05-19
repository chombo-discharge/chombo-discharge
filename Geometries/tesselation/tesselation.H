/*!
  @file   tesselation.H
  @brief  Geometry for single TESSELATION-based mesh
  @author Robert Marskar
  @date   Jan .2019
*/

#ifndef _TESSELATION__
#define _TESSELATION_

#include <CD_ComputationalGeometry.H>

#include "CD_NamespaceHeader.H"

class tesselation : public ComputationalGeometry {
public:

  tesselation();
  ~tesselation();
};
#include "CD_NamespaceFooter.H"

#endif
