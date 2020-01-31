/*!
  @file   rod_plane.cpp
  @brief  Implementation of rod_plane.H
  @author Robert Marskar
  @date   Nov. 2017
*/

#include "rod_plane.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

#include "perlin_plane_if.H"
#include "rod_if.H"

rod_plane::rod_plane(){
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  Vector<Real> vec(SpaceDim);
  std::string str;
  bool live, has_electrode, has_dielectric, noise_reseed;
  Real eps0, radius, plane_permittivity, noise_amplitude, noise_persistence;
  RealVect center1, center2, plane_normal, plane_point, noise_frequency;
  int noise_octaves;

  ParmParse pp("rod_plane");
  pp.get   ("eps0",                    eps0);
  pp.get   ("electrode_radius",        radius);
  pp.get   ("noise_octaves",           noise_octaves);
  pp.get   ("noise_amplitude",         noise_amplitude);
  pp.get   ("noise_persistence",       noise_persistence);
  pp.get   ("plane_permittivity",      plane_permittivity);
  pp.get   ("electrode_live",          str);                live            = (str == "true") ? true  : false;
  pp.get   ("noise_reseed",            str);                noise_reseed    = (str == "true") ? true  : false;
  pp.get   ("turn_off_electrode",      str);                has_electrode   = (str == "true") ? false : true;
  pp.get   ("turn_off_dielectric",     str);                has_dielectric  = (str == "true") ? false : true;
  pp.getarr("electrode_center1",       vec, 0, SpaceDim);   center1         = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("electrode_center2",       vec, 0, SpaceDim);   center2         = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("plane_point",             vec, 0, SpaceDim);   plane_point     = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("plane_normal",            vec, 0, SpaceDim);   plane_normal    = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("noise_frequency",         vec, 0, SpaceDim);   noise_frequency = RealVect(D_DECL(vec[0], vec[1], vec[2]));


  // Make geometry
  this->set_eps0(eps0);
  if(has_electrode){
    m_electrodes.resize(1);
    RefCountedPtr<BaseIF> rod  = RefCountedPtr<BaseIF> (new rod_if(center1, center2, radius, false));
    m_electrodes[0].define(rod,   live);
  }
  if(has_dielectric){
    m_dielectrics.resize(1);
    RefCountedPtr<BaseIF> plane = RefCountedPtr<BaseIF> (new perlin_plane_if(plane_normal,
									     plane_point,
									     true,
									     noise_amplitude,
									     noise_frequency,
									     noise_persistence,
									     noise_octaves,
									     noise_reseed));
    m_dielectrics[0].define(plane, plane_permittivity);
  }
}

rod_plane::~rod_plane(){}
