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
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  bool rod_live  = true;
  bool has_rod   = true;
  bool has_plane = true;

  Real eps     = 4.0;
  Real eps0    = 1.0;
  Real rod_rad = 1.E-2;
  Real rad     = 1.E-2;
  Real offset  = 0.0;
  Real curv    = 1.E-4;
  Real dist    = 1.E-2;
  Real shift   = 0.0;
  Real width   = 1.0;
  
  int numl = 1;
  int numr = 1;

  RealVect center1 = -RealVect::Unit;
  RealVect center2 = RealVect::Unit;
  RealVect point   = RealVect::Zero;
  RealVect normal  = RealVect::Unit;

  ParmParse pp("rod_plane_profile");
  
  pp.get("rod_radius",        rod_rad);
  pp.get("profile_num_left",  numl);
  pp.get("profile_num_right", numr);
  pp.get("profile_radius",    rad);
  pp.get("profile_offset",    offset);
  pp.get("profile_curv",      curv);
  pp.get("profile_dist",      dist);
  pp.get("profile_shift",     shift);
  pp.get("plane_width",       width);
  pp.get("plane_eps",         eps);

  // Rod center, plane point, and plane normal
  pp.getarr("rod_center1",   vec, 0, SpaceDim); center1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("rod_center2",   vec, 0, SpaceDim); center2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("plane_point",   vec, 0, SpaceDim); point   = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("plane_normal",  vec, 0, SpaceDim); normal  = RealVect(D_DECL(vec[0], vec[1], vec[2]));

  // Rod on/off
  pp.get("rod_live",       str); rod_live  = (str == "true") ? true  : false;
  pp.get("turn_off_rod",   str); has_rod   = (str == "true") ? false : true;
  pp.get("turn_off_plane", str); has_plane = (str == "true") ? false : true;

  // Build geometries
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
