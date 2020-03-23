/*!
  @file   rod_plane_profile.cpp
  @brief  Implementation of rod_plane_profile.H
  @author Robert Marskar
  @date   Nov. 2017
*/

#include "rod_plane_profile.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

#include "profile_plane_if.H"
#include "rod_if.H"

rod_plane_profile::rod_plane_profile(){
  if(SpaceDim == 3) MayDay::Abort("rod_plane_profile::rod_plane_profile - this is currently for 2D only");
  std::string str;
  Vector<Real> vec(SpaceDim);

  bool rod_live, has_rod, has_plane;
  Real eps, eps0, rod_rad, rad, offset, curv, dist, shift, width;
  int numl, numr;
  RealVect center1, center2, point, normal;

  ParmParse pp("rod_plane_profile");
  
  pp.get   ("rod_radius",        rod_rad);
  pp.get   ("profile_num_left",  numl);
  pp.get   ("profile_num_right", numr);
  pp.get   ("profile_radius",    rad);
  pp.get   ("profile_offset",    offset);
  pp.get   ("profile_curv",      curv);
  pp.get   ("profile_dist",      dist);
  pp.get   ("profile_shift",     shift);
  pp.get   ("plane_width",       width);
  pp.get   ("plane_eps",         eps);
  pp.get   ("rod_live",          str);              rod_live  = (str == "true") ? true  : false;
  pp.get   ("turn_off_rod",      str);              has_rod   = (str == "true") ? false : true;
  pp.get   ("turn_off_plane",    str);              has_plane = (str == "true") ? false : true;
  pp.getarr("rod_center1",       vec, 0, SpaceDim); center1   = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("rod_center2",       vec, 0, SpaceDim); center2   = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("plane_point",       vec, 0, SpaceDim); point     = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("plane_normal",      vec, 0, SpaceDim); normal    = RealVect(D_DECL(vec[0], vec[1], vec[2]));

  // Build geometries
  m_dielectrics.resize(0);
  m_electrodes.resize(0);
  if(has_rod){
    m_electrodes.resize(1);
    RefCountedPtr<BaseIF> rod  = RefCountedPtr<BaseIF> (new rod_if(center1, center2, rod_rad, false));
    m_electrodes[0].define(rod, rod_live);
  }
  if(has_plane){
    m_dielectrics.resize(1);
    RefCountedPtr<BaseIF> plane = RefCountedPtr<BaseIF> (new profile_plane_if(point,width,numl,numr,rad,offset,shift,dist,curv, true));
    m_dielectrics[0].define(plane, eps);
  }

  this->set_eps0(1.0);
}

rod_plane_profile::~rod_plane_profile(){
  
}
