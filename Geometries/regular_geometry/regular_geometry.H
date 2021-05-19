/*!
  @file   regular_geometry.H
  @brief  Regular geometry (no EBs)
  @author Robert Marskar
  @date   July 2018
*/

#ifndef _REGULAR_GEOMETRY_GEOM_
#define _REGULAR_GEOMETRY_GEOM_

#include <CD_ComputationalGeometry.H>

#include "CD_NamespaceHeader.H"
  
/*!
  @brief This class implements an empty geometry
*/
class regular_geometry : public ComputationalGeometry {
public:

  regular_geometry();
  ~regular_geometry();
};
#include "CD_NamespaceFooter.H"

#endif
