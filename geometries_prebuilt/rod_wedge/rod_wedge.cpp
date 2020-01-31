/*!
  @file   rod_wedge.cpp
  @brief  Implementation of rod_wedge.H
  @author Robert Marskar
  @date   March 2019
*/

#include "rod_wedge.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

#include "rod_if.H"
#include "wedge_if.H"


rod_wedge::rod_wedge(){
  std::string str;
  Vector<Real> vec(SpaceDim);

  bool rod_live, has_rod, has_wedge;
  int wedge_dir;
  Real eps0, radius, wedge_eps, wedge_angle, wedge_curv;
  RealVect center1, center2, wedge_point;

  ParmParse pp("rod_wedge");
  pp.get   ("eps0",           eps0);
  pp.get   ("rod_radius",     radius);
  pp.get   ("wedge_eps",      wedge_eps);
  pp.get   ("wedge_dir",      wedge_dir);
  pp.get   ("wedge_angle",    wedge_angle);
  pp.get   ("wedge_curv",     wedge_curv);
  pp.get   ("turn_off_rod",   str);              has_rod     = (str == "true") ? false : true;
  pp.get   ("turn_off_wedge", str);              has_wedge   = (str == "true") ? false : true;
  pp.getarr("rod_center1",    vec, 0, SpaceDim); center1     = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("rod_center2",    vec, 0, SpaceDim); center2     = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("wedge_point",    vec, 0, SpaceDim); wedge_point = RealVect(D_DECL(vec[0], vec[1], vec[2]));

  // build geometries
  if(has_rod){
    m_electrodes.resize(1);
    RefCountedPtr<BaseIF> rod  = RefCountedPtr<BaseIF> (new rod_if(center1, center2, radius, false));
    m_electrodes[0].define(rod, rod_live);
  }
  if(has_wedge){
    m_dielectrics.resize(1);
    RefCountedPtr<BaseIF> wedge = RefCountedPtr<BaseIF> (new wedge_if(wedge_dir, wedge_angle, wedge_curv, wedge_point, true));
    m_dielectrics[0].define(wedge, wedge_eps);
  }
  
  this->set_eps0(eps0);
}

rod_wedge::~rod_wedge(){
  
}
