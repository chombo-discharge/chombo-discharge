/*!
  @file   rod_slab.H
  @brief  Definition of a electrode rod and dielectric slab geometry
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _ROD_SLAB_
#define _ROD_SLAB_

#include "computational_geometry.H"

/*!
  @brief Class that defines an electrode rod and a dielectric slab
*/
class rod_slab : public computational_geometry {
public:

  rod_slab();
  ~rod_slab();
};

#endif
