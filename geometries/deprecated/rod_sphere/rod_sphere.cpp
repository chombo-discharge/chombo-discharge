/*!
  @file   rod_sphere.cpp
  @brief  Implementation of rod_sphere.H
  @author Robert Marskar
  @date   Nov. 2017
*/

#include "rod_sphere.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

#include "rounded_box_if.H"
#include "rod_if.H"
#include "new_sphere_if.H"

rod_sphere::rod_sphere(){
  m_dielectrics.resize(0);
  m_electrodes.resize(0);
  

  std::string str;
  Vector<Real> vec(SpaceDim);
  bool live, has_electrode, has_dielectric;
  Real eps0, elec_radius, radius, eps;
  RealVect center1, center2, centerD;

  ParmParse pp("rod_sphere");
  pp.get("eps0",eps0);
  pp.get("dielectric_permittivity", eps);
  pp.get("dielectric_radius", radius);
  pp.get("electrode_radius", elec_radius);
  pp.getarr("electrode_center1", vec, 0, SpaceDim); center1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("electrode_center2", vec, 0, SpaceDim); center2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("dielectric_center", vec, 0, SpaceDim); centerD = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.get("turn_off_electrode", str);                has_electrode  = (str == "true") ? false : true;
  pp.get("turn_off_dielectric", str);               has_dielectric = (str == "true") ? false : true;
  pp.get("electrode_live", str);                    live = (str == "true") ? true : false;

  set_eps0(eps0);
  if(has_electrode){
    m_electrodes.resize(1);
    RefCountedPtr<BaseIF> rod  = RefCountedPtr<BaseIF> (new rod_if(center1, center2, elec_radius, false));
    m_electrodes[0].define(rod,   live);
  }
  if(has_dielectric){
    m_dielectrics.resize(1);
    RefCountedPtr<BaseIF> slab = RefCountedPtr<BaseIF> (new new_sphere_if(centerD, radius, false));
    m_dielectrics[0].define(slab, eps);
  }
}

rod_sphere::~rod_sphere(){
}

