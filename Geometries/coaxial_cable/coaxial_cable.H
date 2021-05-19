/*!
  @file   coaxial_cable.H
  @brief  Coaxial cable geometry
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _COAXIAL_CABLE_
#define _COAXIAL_CABLE_

#include <CD_ComputationalGeometry.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief This is a class that defines a coaxial cable geometry 
*/
class coaxial_cable : public ComputationalGeometry {
public:
  coaxial_cable();
  ~coaxial_cable();
};
#include "CD_NamespaceFooter.H"

#endif
