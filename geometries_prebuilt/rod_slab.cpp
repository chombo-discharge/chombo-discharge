/*!
  @file rod_slab_geometry.cpp
  @brief Implementation of rod_slab_geometry.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "rod_slab.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

#include "rounded_box_if.H"
#include "rod_if.H"
#include "new_sphere_if.H"

rod_slab::rod_slab(){
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  Real eps0        = 1.0;
  Real corner_curv = 1.E-3;
  {
    ParmParse pp("rod_slab");
    pp.get("eps0",              eps0);
    pp.get("corner_curvatures", corner_curv);
  }
  this->set_eps0(eps0);
  
  // Electrode parameters
  bool live          = true;
  bool has_electrode = true;
  Real radius        = 1.E-2;
  RealVect center1   = RealVect::Zero;
#if CH_SPACEDIM == 2
  RealVect center2   = RealVect(0.0, 1.0);
#else
  RealVect center2   = RealVect(0.0, 0.0, 1.0);
#endif

  { // Get parameterse for electrode rod
    ParmParse pp("rod_slab");
    std::string str;
    Vector<Real> vec(SpaceDim);
    pp.get("electrode_radius", radius);
    if(pp.contains("electrode_center1")){
      pp.getarr("electrode_center1", vec, 0, SpaceDim);
      center1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    if(pp.contains("electrode_center2")){
      pp.getarr("electrode_center2", vec, 0, SpaceDim);
      center2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    pp.get("turn_off_electrode", str);
    if(str == "true"){
      has_electrode = false;
    }
    else if(str == "false"){
      has_electrode = true;
    }
    pp.get("electrode_live", str);
    if(str == "true"){
      live = true;
    }
    else if(str == "false"){
      live = false;
    }
  }
  //Slab
  bool has_dielectric = true;
  Real dielectric_permittivity = 5.0;
#if CH_SPACEDIM == 2
  RealVect slab_lo = RealVect(-2.0123E-2, -2.0123E-2);
  RealVect slab_hi = RealVect(2.0123E-2,  -1.0123E-2);
#else
  RealVect slab_lo = RealVect(-1.0123E-2, -2.0123E-2,  -2.0123E-2);
  RealVect slab_hi = RealVect( 1E0,        2.0123E-2,  -1.5123E-2);
#endif

  {
    ParmParse pp("rod_slab");
    Vector<Real> vec(SpaceDim);
    std::string str;
    pp.get("dielectric_permittivity", dielectric_permittivity);
    if(pp.contains("dielectric_corner_lo")){
      pp.getarr("dielectric_corner_lo", vec, 0, SpaceDim);
      slab_lo = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    if(pp.contains("dielectric_corner_hi")){
      pp.getarr("dielectric_corner_hi", vec, 0, SpaceDim);
      slab_hi = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    pp.get("turn_off_dielectric", str);
    if(str == "true"){
      has_dielectric = false;
    }
    else if(str == "false"){
      has_dielectric = true;
    }
  }

  if(has_electrode){
    m_electrodes.resize(1);
    RefCountedPtr<BaseIF> rod  = RefCountedPtr<BaseIF> (new rod_if(center1, center2, radius, false));
    m_electrodes[0].define(rod,   live);
  }
  if(has_dielectric){
    m_dielectrics.resize(1);
    RefCountedPtr<BaseIF> slab = RefCountedPtr<BaseIF> (new rounded_box_if(slab_lo, slab_hi, corner_curv, false));
    m_dielectrics[0].define(slab, dielectric_permittivity);
  }
}

rod_slab::~rod_slab(){
}

