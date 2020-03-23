/*!
  @file   rod_multisphere.H
  @brief  Declaration of a rod electrode and multiple spheres
  @author Robert Marskar
  @date   Aug. 2018
*/

#ifndef _ROD_MULTISPHERE_
#define _ROD_MULTISPHERE_

#include "computational_geometry.H"

/*!
  @brief Computational geometry class that describes a rod electrode and multiple dielectric spheres
*/
class rod_multisphere : public computational_geometry {
public:
  rod_multisphere();

  ~rod_multisphere();
};

#endif
