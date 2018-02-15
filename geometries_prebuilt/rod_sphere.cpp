/*!
  @file rod_sphere_geometry.cpp
  @brief Implementation of rod_sphere_geometry.H
  @author Robert Marskar
  @date Nov. 2017
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
  
  Real eps0        = 1.0;
  {
    ParmParse pp("rod_sphere");
    pp.query("eps0",              eps0);
  }
  this->set_eps0(eps0);
  
  // Electrode parameters
  bool live          = true;
  bool has_electrode = true;
  Real elec_radius   = 1.E-2;
  RealVect center1   = RealVect::Zero;
#if CH_SPACEDIM == 2
  RealVect center2   = RealVect(0.0, 1.0);
#else
  RealVect center2   = RealVect(0.0, 0.0, 1.0);
#endif

  { // Get parameterse for electrode rod
    ParmParse pp("rod_sphere");
    std::string str;
    Vector<Real> vec(SpaceDim);
    pp.query("electrode_radius", elec_radius);
    if(pp.contains("electrode_center1")){
      pp.getarr("electrode_center1", vec, 0, SpaceDim);
      center1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    if(pp.contains("electrode_center2")){
      pp.getarr("electrode_center2", vec, 0, SpaceDim);
      center2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    pp.query("turn_off_electrode", str);
    if(str == "true"){
      has_electrode = false;
    }
    else if(str == "false"){
      has_electrode = true;
    }
    pp.query("electrode_live", str);
    if(str == "true"){
      live = true;
    }
    else if(str == "false"){
      live = false;
    }
  }
  
  //Dielectric sphere
  bool has_dielectric = true;
  Real dielectric_permittivity = 5.0;
  Real radius = 1.0;
  RealVect center = RealVect::Zero;

  {
    ParmParse pp("rod_sphere");
    Vector<Real> vec(SpaceDim);
    std::string str;
    pp.query("dielectric_permittivity", dielectric_permittivity);

    if(pp.contains("dielectric_center")){
      pp.getarr("dielectric_center", vec, 0, SpaceDim);
      center = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    pp.query("dielectric_radius", radius);
    pp.query("turn_off_dielectric", str);
    if(str == "true"){
      has_dielectric = false;
    }
    else if(str == "false"){
      has_dielectric = true;
    }
  }

  if(has_electrode){
    m_electrodes.resize(1);
    RefCountedPtr<BaseIF> rod  = RefCountedPtr<BaseIF> (new rod_if(center1, center2, elec_radius, false));
    m_electrodes[0].define(rod,   live);
  }
  if(has_dielectric){
    m_dielectrics.resize(1);
    RefCountedPtr<BaseIF> slab = RefCountedPtr<BaseIF> (new new_sphere_if(center, radius, false));
    m_dielectrics[0].define(slab, dielectric_permittivity);
  }
}

rod_sphere::~rod_sphere(){
}

