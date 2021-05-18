/*!
  @file   coaxial_cable.H
  @brief  Coaxial cable geometry
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _COAXIAL_CABLE_
#define _COAXIAL_CABLE_

#include "computational_geometry.H"

#include "CD_NamespaceHeader.H"

/*!
  @brief This is a class that defines a coaxial cable geometry 
*/
class coaxial_cable : public computational_geometry {
public:
  coaxial_cable();
  ~coaxial_cable();
};
#include "CD_NamespaceFooter.H"

#endif
