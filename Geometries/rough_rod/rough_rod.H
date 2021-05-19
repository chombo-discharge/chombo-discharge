/*!
  @file   rough_rod.H
  @brief  Rod geometry
  @author Robert Marskar
  @date   Oct. 2018
*/

#ifndef _ROUGH_ROD_GEOM_
#define _ROUGH_ROD_GEOM_

#include <CD_ComputationalGeometry.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Class that describes a rod electrode with surface roughness
*/
class rough_rod : public ComputationalGeometry {
public:

  rough_rod();
  ~rough_rod();
};
#include "CD_NamespaceFooter.H"

#endif
