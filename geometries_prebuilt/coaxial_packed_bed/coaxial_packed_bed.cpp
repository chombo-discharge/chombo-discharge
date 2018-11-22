/*!
  @file coaxial_packed_bed.cpp
  @brief Implementation of coaxial_packed_bed.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "coaxial_packed_bed.H"

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

coaxial_packed_bed::coaxial_packed_bed(){

  Real eps0 = 1.0;
  Real eps1 = 4.0;
  Real R0   = 0.9;
  Real R2   = 0.3;

  RealVect c0  = RealVect(D_DECL(0.0, 0.0, -10.0));
  RealVect c00 = RealVect(D_DECL(0.0, 0.0,  10.0));
  RealVect c1  = RealVect(D_DECL(0.0, 0.0, -10.0));
  RealVect c11 = RealVect(D_DECL(0.0, 0.0,  10.0));

  bool turn_off_outer      = false;
  bool turn_off_inner      = false;
  bool inner_live          = true;
  bool outer_live          = false;
  
  
  { // Get parameters from input script
    ParmParse pp("coaxial_packed_bed");

    std::string str;
    Vector<Real> v0(SpaceDim), v00(SpaceDim);
    Vector<Real> v1(SpaceDim), v11(SpaceDim);
    Vector<Real> v2(SpaceDim), v22(SpaceDim);
    
    pp.query("eps0",                 eps0);
    pp.query("eps1",                 eps1);
    pp.query("outer_radius",         R0);
    pp.query("inner_radius",         R2);

    pp.queryarr("outer_center1",      v0,  0, SpaceDim);
    pp.queryarr("outer_center2",      v00, 0, SpaceDim);
    pp.queryarr("inner_center1",      v2, 0, SpaceDim);
    pp.queryarr("inner_center2",      v22, 0, SpaceDim);

    c0  = RealVect(D_DECL(v0[0],  v0[1],  v0[2]));
    c00 = RealVect(D_DECL(v00[0], v00[1], v00[2]));
    c1  = RealVect(D_DECL(v2[0],  v2[1],  v2[2]));
    c11 = RealVect(D_DECL(v22[0], v22[1], v22[2]));

    if(pp.contains("turn_off_outer")){
      pp.get("turn_off_outer", str);
      turn_off_outer = (str == "true") ? true : false;
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
    RefCountedPtr<BaseIF> outer = RefCountedPtr<BaseIF> (new cylinder_if(c0, c00, R0, true));
    m_electrodes.push_back(electrode(outer, outer_live, 1.0));
  }
  if(!turn_off_inner){
    RefCountedPtr<BaseIF> inner = RefCountedPtr<BaseIF> (new cylinder_if(c1, c11, R2, false));
    m_electrodes.push_back(electrode(inner, inner_live, 1.0));
  }

#if 1 // Debug
  {
  const RealVect center = RealVect(0.5, 0.0, 0.0);
  RefCountedPtr<BaseIF> middle = RefCountedPtr<BaseIF> (new new_sphere_if(center, 0.2, false));
  m_dielectrics.push_back(dielectric(middle, eps1));
  }

  {
  const RealVect center = RealVect(0.0, 0.5, 0.0);
  RefCountedPtr<BaseIF> middle = RefCountedPtr<BaseIF> (new new_sphere_if(center, 0.2, false));
  m_dielectrics.push_back(dielectric(middle, eps1));
  }
  {
  const RealVect center = RealVect(0.0, -0.5, 0.0);
  RefCountedPtr<BaseIF> middle = RefCountedPtr<BaseIF> (new new_sphere_if(center, 0.2, false));
  m_dielectrics.push_back(dielectric(middle, eps1));
  }
  {
    const RealVect center = RealVect(-0.5, 0.0, 0.0);
    RefCountedPtr<BaseIF> middle = RefCountedPtr<BaseIF> (new new_sphere_if(center, 0.2, false));
    m_dielectrics.push_back(dielectric(middle, eps1));
  }
#endif

#if 1 // Debug
  {
    const RealVect center = RealVect(0.5/sqrt(2), 0.5/sqrt(2), 0.5);
    RefCountedPtr<BaseIF> middle = RefCountedPtr<BaseIF> (new new_sphere_if(center, 0.2, false));
    m_dielectrics.push_back(dielectric(middle, eps1));
  }

  {
    const RealVect center = RealVect(-0.5/sqrt(2), 0.5/sqrt(2), 0.5);
    RefCountedPtr<BaseIF> middle = RefCountedPtr<BaseIF> (new new_sphere_if(center, 0.2, false));
    m_dielectrics.push_back(dielectric(middle, eps1));
  }
  {
    const RealVect center = RealVect(0.5/sqrt(2), -0.5/sqrt(2), 0.5);
    RefCountedPtr<BaseIF> middle = RefCountedPtr<BaseIF> (new new_sphere_if(center, 0.2, false));
    m_dielectrics.push_back(dielectric(middle, eps1));
  }
  {
    const RealVect center = RealVect(-0.5/sqrt(2), -0.5/sqrt(2), 0.5);
    RefCountedPtr<BaseIF> middle = RefCountedPtr<BaseIF> (new new_sphere_if(center, 0.2, false));
    m_dielectrics.push_back(dielectric(middle, eps1));
  }
#endif

#if 1 // Debug
  {
    const RealVect center = RealVect(0.5/sqrt(2), 0.5/sqrt(2), -0.5);
    RefCountedPtr<BaseIF> middle = RefCountedPtr<BaseIF> (new new_sphere_if(center, 0.2, false));
    m_dielectrics.push_back(dielectric(middle, eps1));
  }

  {
    const RealVect center = RealVect(-0.5/sqrt(2), 0.5/sqrt(2), -0.5);
    RefCountedPtr<BaseIF> middle = RefCountedPtr<BaseIF> (new new_sphere_if(center, 0.2, false));
    m_dielectrics.push_back(dielectric(middle, eps1));
  }
  {
    const RealVect center = RealVect(0.5/sqrt(2), -0.5/sqrt(2), -0.5);
    RefCountedPtr<BaseIF> middle = RefCountedPtr<BaseIF> (new new_sphere_if(center, 0.2, false));
    m_dielectrics.push_back(dielectric(middle, eps1));
  }
  {
    const RealVect center = RealVect(-0.5/sqrt(2), -0.5/sqrt(2), -0.5);
    RefCountedPtr<BaseIF> middle = RefCountedPtr<BaseIF> (new new_sphere_if(center, 0.2, false));
    m_dielectrics.push_back(dielectric(middle, eps1));
  }
#endif
  
}

coaxial_packed_bed::~coaxial_packed_bed(){
  
}
