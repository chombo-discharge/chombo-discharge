/*!
  @file   aerosol.H
  @brief  Sphere geometry
  @author Robert Marskar
  @date   Oct. 2018
*/

#ifndef _AEROSOL_GEOM_
#define _AEROSOL_GEOM_

#include "computational_geometry.H"

namespace ChomboDischarge {

  class aerosol : public computational_geometry {
  public:

    aerosol();
    ~aerosol();
  };
}

#endif
