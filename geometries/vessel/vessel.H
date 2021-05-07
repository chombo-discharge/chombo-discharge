/*!
  @file vessel.H
  @brief High-voltage vessel geometry
  @author Robert Marskar
  @date   2019
*/

#ifndef _VESSEL_H_
#define _VESSEL_H_

#include "computational_geometry.H"

namespace ChomboDischarge {

  class vessel : public computational_geometry {
  public:
    vessel();
    ~vessel();
  protected:

    RealVect m_rod_center;
    RealVect m_shroom_center;
    Real m_rod_R;
    Real m_shroom_R;
    Real m_shroom_r;
    Real m_shroom_d;
    Real m_shroom_c;

    bool m_rod_live;
    bool m_shroom_live;

    bool m_shroom;
    bool m_rod;
  };
}

#endif
