/*!
  @file   tesselation.H
  @brief  Geometry for single TESSELATION-based mesh
  @author Robert Marskar
  @date   Jan .2019
*/

#ifndef _TESSELATION__
#define _TESSELATION_

#include "computational_geometry.H"

#include "CD_NamespaceHeader.H"

class tesselation : public computational_geometry {
public:

  tesselation();
  ~tesselation();
};
#include "CD_NamespaceFooter.H"

#endif
