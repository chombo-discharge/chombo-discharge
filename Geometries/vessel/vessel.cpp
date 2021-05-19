/*!
  @file   vessel.cpp
  @brief  High-voltage vessel geometry
  @author Robert Marskar
  @date   2019
*/

#include "vessel.H"
#include <CD_MushroomIF.H>
#include "rod_if.H"
#include <CD_CylinderSdf.H>

#include <ParmParse.H>

#include "CD_NamespaceHeader.H"

vessel::vessel(){
  ParmParse pp("vessel");

  std::string str;
  Vector<Real> vec(SpaceDim);

  pp.get("rod_radius",  m_rod_R);
  pp.get("shroom_R",    m_shroom_R);
  pp.get("shroom_r",    m_shroom_r);
  pp.get("shroom_d",    m_shroom_d);
  pp.get("shroom_curv", m_shroom_c);
  pp.get("use_rod",     m_rod);
  pp.get("use_shroom",  m_shroom);
  pp.get("live_rod",    m_rod_live);
  pp.get("live_shroom", m_shroom_live);
  pp.getarr("rod_point",    vec, 0, SpaceDim); m_rod_center    = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("shroom_point", vec, 0, SpaceDim); m_shroom_center = RealVect(D_DECL(vec[0], vec[1], vec[2]));

  const RealVect v  = RealVect(BASISV(SpaceDim-1));
  
  if(m_rod){
    auto rod = RefCountedPtr<BaseIF> (new rod_if(m_rod_center + 100*v, m_rod_center, m_rod_R, false));
    m_electrodes.push_back(Electrode(rod, m_rod_live));
  }
  if(m_shroom){
    auto shroom = RefCountedPtr<BaseIF> (new MushroomIF(m_shroom_center,
							 m_shroom_R,
							 m_shroom_r,
							 1.E4,
							 m_shroom_d,
							 m_shroom_c,
							 false));
    m_electrodes.push_back(Electrode(shroom, m_shroom_live));
  }
}

vessel::~vessel(){

}
#include "CD_NamespaceFooter.H"
