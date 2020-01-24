/*!
  @file   fast_rod_plane.cpp
  @brief  Implementation of fast_rod_plane.H
  @author Robert Marskar
  @date   Nov. 2017
*/

#include "fast_rod_plane.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

#include "perlin_plane_if.H"
#include "rod_if.H"

fast_rod_plane::fast_rod_plane(){
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  std::string str;
  std::string normal;
  Vector<Real> vec(SpaceDim);
  Real eps0, eps, radius, noise_amplitude, noise_persistence;
  RealVect center1, center2, plane_normal, ortho_normal, plane_point, noise_frequency;
  bool live, has_electrode, has_dielectric, noise_reseed, use_new_gshop;
  int noise_octaves;
  
  ParmParse pp("fast_rod_plane");
  pp.get   ("eps0",                eps0);
  pp.get   ("permittivity",        eps);
  pp.get   ("electrode_radius",    radius);
  pp.get   ("plane_normal",        normal);
  pp.get   ("noise_amplitude",     noise_amplitude);
  pp.get   ("noise_persistence",   noise_persistence);
  pp.get   ("noise_octaves",       noise_octaves);
  pp.get   ("noise_reseed",        str);              noise_reseed    = (str == "true") ? true : false;
  pp.get   ("use_new_gshop",       str);              use_new_gshop   = (str == "true") ? true : false;
  pp.get   ("turn_off_electrode",  str);              has_electrode   = (str == "true") ? true : false;
  pp.get   ("turn_off_dielectric", str);              has_dielectric  = (str == "true") ? true : false;
  pp.get   ("electrode_live",      str);              live            = (str == "true") ? true : false;
  pp.getarr("noise_frequency",     vec, 0, SpaceDim); noise_frequency = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("plane_point",         vec, 0, SpaceDim); plane_point     = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("electrode_center1",   vec, 0, SpaceDim); center1         = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("electrode_center2",   vec, 0, SpaceDim); center2         = RealVect(D_DECL(vec[0], vec[1], vec[2]));  

  set_eps0(eps0);

  m_electrodes.resize(0);
  m_dielectrics.resize(0);

  // Set the plane normal
  if(normal == "+x"){
    plane_normal = RealVect(D_DECL(1,0,0));
    ortho_normal = RealVect(D_DECL(0,1,1));
  }
  else if(normal == "-x"){
    plane_normal = -RealVect(D_DECL(1,0,0));
    ortho_normal =  RealVect(D_DECL(0,1,1));
  }
  else if(normal == "+y"){
    plane_normal = RealVect(D_DECL(0,1,0));
    ortho_normal = RealVect(D_DECL(1,0,1));
  }
  else if(normal == "-y"){
    plane_normal = -RealVect(D_DECL(0,1,0));
    ortho_normal =  RealVect(D_DECL(1,0,1));
  }
  else if(normal == "+z"){
    plane_normal = RealVect(D_DECL(0,0,1));
    ortho_normal = RealVect(D_DECL(1,1,0));
  }
  else if(normal == "-z"){
    plane_normal = -RealVect(D_DECL(0,0,1));
    ortho_normal =  RealVect(D_DECL(1,1,0));
  }

  
  RefCountedPtr<BaseIF> rod   = RefCountedPtr<BaseIF> (new rod_if(center1, center2, radius, false));
  RefCountedPtr<BaseIF> plane = RefCountedPtr<BaseIF> (new perlin_plane_if(plane_normal,
									   plane_point,
									   true,
									   noise_amplitude,
									   noise_frequency,
									   noise_persistence,
									   noise_octaves,
									   noise_reseed));

  if(!has_electrode){
    m_electrodes.push_back(electrode(rod,   live));
  }
  if(!has_dielectric){
    m_dielectrics.push_back(dielectric(plane, eps));
  }


  if(use_new_gshop){
    computational_geometry::s_use_new_gshop = true;

    RealVect lo, hi;
    for (int dir = 0; dir < SpaceDim; dir++){
      lo[dir] = Min(center1[dir], center2[dir]);
      hi[dir] = Max(center1[dir], center2[dir]);
    }

    const Real     thresh      = 1.E-6;
    const RealVect point_above = plane_point + noise_octaves*(noise_amplitude + thresh)*plane_normal;
    const RealVect point_below = plane_point - noise_octaves*(noise_amplitude + thresh)*plane_normal;

    const real_box gas_box(point_above - 1.E10*ortho_normal,                 point_above + 1.E10*RealVect::Unit);
    const real_box sol_box(point_below - 1.E10*(ortho_normal+plane_normal),  point_below + 1.E10*ortho_normal);
    const real_box rod_gas(lo-radius*RealVect::Unit,                         hi+radius*RealVect::Unit);
    const real_box diel_surf(point_below - 1.10*ortho_normal, point_above + 1.E10*ortho_normal);
    
    m_regular_voxels_gas.push_back(gas_box);
    m_covered_voxels_gas.push_back(sol_box);
    m_bounded_voxels_gas.push_back(rod_gas);
    m_bounded_voxels_gas.push_back(diel_surf);

    m_regular_voxels_sol.push_back(sol_box);
    m_covered_voxels_sol.push_back(gas_box);
    m_bounded_voxels_sol.push_back(diel_surf);
  }
} 

fast_rod_plane::~fast_rod_plane(){
  
}
