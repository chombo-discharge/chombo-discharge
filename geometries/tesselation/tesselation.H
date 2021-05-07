/*!
  @file   tesselation.H
  @brief  Geometry for single TESSELATION-based mesh
  @author Robert Marskar
  @date   Jan .2019
*/

#ifndef _TESSELATION__
#define _TESSELATION_

#include "computational_geometry.H"

namespace ChomboDischarge {

  class tesselation : public computational_geometry {
  public:

    tesselation();
    ~tesselation();
  };
}

#endif
