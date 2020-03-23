/*!
  @file   vessel.cpp
  @brief  High-voltage vessel geometry
  @author Robert Marskar
  @date   2019
*/

#include "vessel.H"
#include "mushroom_if.H"
#include "rod_if.H"
#include "cylinder_if.H"

#include <ParmParse.H>

vessel::vessel(){
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  ParmParse pp("vessel");

  std::string str;
  Vector<Real> vec(SpaceDim);

  pp.getarr("rod_point",    vec, 0, SpaceDim); m_rod_center    = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("shroom_point", vec, 0, SpaceDim); m_shroom_center = RealVect(D_DECL(vec[0], vec[1], vec[2]));

  pp.get("rod_radius",  m_rod_R);
  pp.get("shroom_R",    m_shroom_R);
  pp.get("shroom_r",    m_shroom_r);
  pp.get("shroom_d",    m_shroom_d);
  pp.get("shroom_curv", m_shroom_c);

  pp.get("turn_off_rod",    str); m_rod    = (str == "true") ? false : true;
  pp.get("turn_off_shroom", str); m_shroom = (str == "true") ? false : true;

  pp.get("live_rod",    str); m_rod_live    = (str == "true") ? true : false;
  pp.get("live_shroom", str); m_shroom_live = (str == "true") ? true : false;

  const RealVect v  = RealVect(BASISV(SpaceDim-1));
  if(m_rod){
    auto rod = RefCountedPtr<BaseIF> (new rod_if(m_rod_center + 100*v, m_rod_center, m_rod_R, false));
    m_electrodes.push_back(electrode(rod, m_rod_live));
  }
  if(m_shroom){
    auto shroom = RefCountedPtr<BaseIF> (new mushroom_if(m_shroom_center,
							 m_shroom_R,
							 m_shroom_r,
							 1.E4,
							 m_shroom_d,
							 m_shroom_c,
							 true));
    m_electrodes.push_back(electrode(shroom, m_shroom_live));
  }
}

vessel::~vessel(){

}

