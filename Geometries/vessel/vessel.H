/*!
  @file vessel.H
  @brief High-voltage vessel geometry
  @author Robert Marskar
  @date   2019
*/

#ifndef _VESSEL_H_
#define _VESSEL_H_

#include <CD_ComputationalGeometry.H>

#include "CD_NamespaceHeader.H"

class vessel : public ComputationalGeometry {
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
#include "CD_NamespaceFooter.H"

#endif
