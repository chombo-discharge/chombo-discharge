/*!
  @file rough_sphere.cpp
  @brief Implementation of rough_sphere.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "rough_sphere.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>

#include "perlin_sphere_if.H"

rough_sphere::rough_sphere(){

  Real eps0       = 1.0;
  Real R0         = 1.0;
  Real noise_pers = 0.5;
  Real noise_amp  = 0.1;
  RealVect freque = RealVect::Unit;
  RealVect c0     = RealVect::Zero;
  int noise_octav = 1;
  bool reseed     = false;

  bool turn_off_sphere = false;
  bool live            = true;
  
  
  { // Get parameters from input script
    ParmParse pp("rough_sphere");

    std::string str;
    Vector<Real> v0(SpaceDim);
    
    pp.query("eps0",          eps0);
    pp.query("sphere_radius", R0);
    pp.query("noise_amplitude", noise_amp);
    pp.query("noise_octaves", noise_octav);
    pp.query("noise_persistence", noise_pers);

    // Sphere center
    pp.queryarr("sphere_center",   v0, 0, SpaceDim);
    c0 = RealVect(D_DECL(v0[0], v0[1], v0[2]));
    

    pp.queryarr("noise_frequency", v0, 0, SpaceDim);
    freque = RealVect(D_DECL(v0[0], v0[1], v0[2]));

    if(pp.contains("turn_off_sphere")){
      pp.get("turn_off_sphere", str);
      turn_off_sphere = (str == "true") ? true : false;
    }
    if(pp.contains("live")){
      pp.get("live", str);
      live = (str == "true") ? true : false;
    }
    if(pp.contains("noise_reseed")){
      pp.get("noise_reseed", str);
      live = (str == "true") ? true : false;
    }
  }

  // Create geometry
  this->set_eps0(eps0);
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  if(!turn_off_sphere){
    RefCountedPtr<BaseIF> outer = RefCountedPtr<BaseIF> (new perlin_sphere_if(R0,
									      c0,
									      true,
									      noise_amp,
									      freque,
									      noise_pers,
									      noise_octav,
									      reseed));
    m_electrodes.push_back(electrode(outer, true, 1.0));
  }
}

rough_sphere::~rough_sphere(){
  
}
