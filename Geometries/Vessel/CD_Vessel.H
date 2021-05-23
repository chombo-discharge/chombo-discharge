/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Vessel.H
  @brief  High-voltage vessel geometry
  @author Robert Marskar
*/

#ifndef CD_Vessel_H
#define CD_Vessel_H

// our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

class Vessel : public ComputationalGeometry {
public:
  Vessel();
  ~Vessel();
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

#include <CD_NamespaceFooter.H>

#endif
