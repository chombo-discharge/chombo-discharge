/*!
  @file coaxial_cable.cpp
  @brief Implementation of coaxial_cable.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "coaxial_cable.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

#include "cylinder_if.H"
#include "perlin_sphere_if.H"
#include "box_if.H"
#include "new_sphere_if.H"
#include "cylinder_if.H"

coaxial_cable::coaxial_cable(){
  if(SpaceDim==3){
  MayDay::Abort("coaxial_cable::coaxial_cable - Only 2D currently supported. If you want 3D, you need to modify this class");
  }

  Real eps0 = 1.0;
  Real eps1 = 4.0;
  Real R0   = 0.9;
  Real R1   = 0.6;
  Real R2   = 0.3;

  RealVect c0 = RealVect(D_DECL(0.0, 0.0, -1.0));
  RealVect c1 = RealVect(D_DECL(0.0, 0.0, -1.0));
  RealVect c2 = RealVect(D_DECL(0.0, 0.0, -1.0));

  bool turn_off_outer      = false;
  bool turn_off_dielectric = false;
  bool turn_off_inner      = false;
  bool inner_live          = true;
  bool outer_live          = false;
  
  
  { // Get parameters from input script
    ParmParse pp("coaxial_cable");

    std::string str;
    Vector<Real> v0(SpaceDim);
    Vector<Real> v1(SpaceDim);
    Vector<Real> v2(SpaceDim);
    
    pp.query("eps0",                 eps0);
    pp.query("eps1",                 eps1);
    pp.query("outer_radius",         R0);
    pp.query("dielectric_radius",    R1);
    pp.query("inner_radius",         R2);

    pp.queryarr("outer_center",      v0, 0, SpaceDim);
    pp.queryarr("dielectric_center", v1, 0, SpaceDim);
    pp.queryarr("inner_center",      v2, 0, SpaceDim);
    c0 = RealVect(D_DECL(v0[0], v0[1], v0[2]));
    c1 = RealVect(D_DECL(v1[0], v1[1], v1[2]));
    c2 = RealVect(D_DECL(v2[0], v2[1], v2[2]));

    if(pp.contains("turn_off_outer")){
      pp.get("turn_off_outer", str);
      turn_off_outer = (str == "true") ? true : false;
    }
    if(pp.contains("turn_off_dielectric")){
      pp.get("turn_off_dielectric", str);
      turn_off_dielectric = (str == "true") ? true : false;
    }
    if(pp.contains("turn_off_inner")){
      pp.get("turn_off_inner", str);
      turn_off_inner = (str == "true") ? true : false;
    }
    if(pp.contains("outer_live")){
      pp.get("outer_live", str);
      outer_live = (str == "true") ? true : false;
    }
    if(pp.contains("inner_live")){
      pp.get("inner_live", str);
      inner_live = (str == "true") ? true : false;
    }
  }

  // Create geometry
  this->set_eps0(eps0);
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  if(!turn_off_outer){
    RefCountedPtr<BaseIF> outer = RefCountedPtr<BaseIF> (new new_sphere_if(c0, R0, true));
    m_electrodes.push_back(electrode(outer, outer_live, 1.0));
  }
  if(!turn_off_dielectric){
    RefCountedPtr<BaseIF> middle = RefCountedPtr<BaseIF> (new new_sphere_if(c1, R1, false));
    m_dielectrics.push_back(dielectric(middle, eps1));
  }
  if(!turn_off_inner){
    RefCountedPtr<BaseIF> inner = RefCountedPtr<BaseIF> (new new_sphere_if(c2, R2, false));
    m_electrodes.push_back(electrode(inner, inner_live, 1.0));
  }
  
}

coaxial_cable::~coaxial_cable(){
  
}
