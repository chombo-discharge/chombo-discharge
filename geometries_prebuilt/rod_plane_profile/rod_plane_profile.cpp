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

  // Check if we should use the new gshop
  bool new_gshop = false;
  pp.get("use_new_gshop", str);
  const bool use_new_gshop = (str == "true") ? true : false;
  if(use_new_gshop){
    computational_geometry::s_use_new_gshop = true;
    

    
    RealVect lo, hi;

    // Make lo, hi cover the needle
    for (int dir = 0; dir < SpaceDim; dir++){
      lo[dir] = Min(center1[dir], center2[dir]);
      hi[dir] = Max(center1[dir], center2[dir]);
    }

    const real_box reg_gas(-100*RealVect::Unit, 100*RealVect::Unit);
    const real_box rod_gas(lo-rod_rad*RealVect::Unit, hi+rod_rad*RealVect::Unit);

    // Make lo, hi cover the plane
    lo = point - 0.5*width*RealVect(BASISV(0));
    hi = point + 0.5*width*RealVect(BASISV(0));
    hi += rad*RealVect::Unit;
    lo -= rad*RealVect::Unit;

    lo -= 100*RealVect(BASISV(1));

    m_regular_voxels_gas.push_back(reg_gas);
    if(has_rod)  m_bounded_voxels_gas.push_back(rod_gas);
    if(has_plane) m_bounded_voxels_gas.push_back(real_box(lo,hi));

    if(has_plane){
      m_covered_voxels_sol.push_back(reg_gas);

      //      std::cout << lo << "\t" << hi << std::endl;
      m_bounded_voxels_sol.push_back(real_box(lo,hi));
    }
  }
}

rod_plane_profile::~rod_plane_profile(){
  
}
