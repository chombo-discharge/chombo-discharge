/*!
  @file rough_rod.cpp
  @brief Implementation of rough_rod.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "rough_rod.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>

#include "perlin_rod_if.H"

rough_rod::rough_rod(){

  Real eps0       = 1.0;
  Real R0         = 1.0;
  Real noise_pers = 0.5;
  Real noise_amp  = 0.1;
  RealVect f      = RealVect::Unit;
  RealVect c0     = RealVect::Zero;
  RealVect c1     = RealVect::Unit;
  
  int noise_octav   = 1;
  bool reseed       = false;
  bool turn_off_rod = false;
  bool live         = true;
  
  
  // Get parameters from input script
  ParmParse pp("rough_rod");

  std::string str;
  Vector<Real> v0(SpaceDim);
  Vector<Real> v1(SpaceDim);
  Vector<Real> v2(SpaceDim);
    
  pp.get("eps0", eps0);
  pp.get("radius", R0);
  pp.get("noise_octaves", noise_octav);
  pp.get("noise_amplitude", noise_amp);
  pp.get("noise_persistence", noise_pers);

  // Rod center and frequency
  pp.getarr("center1",  v0, 0, SpaceDim);
  pp.getarr("center2",  v1, 0, SpaceDim);
  pp.getarr("noise_frequency", v2, 0, SpaceDim);
    
  c0 = RealVect(D_DECL(v0[0], v0[1], v0[2]));
  c1 = RealVect(D_DECL(v1[0], v1[1], v1[2]));
  f  = RealVect(D_DECL(v2[0], v2[1], v2[2]));

  pp.get("turn_off_rod", str);
  turn_off_rod = (str == "true") ? true : false;
    
  pp.get("live", str);
  live = (str == "true") ? true : false;
    
  pp.get("noise_reseed", str);
  live = (str == "true") ? true : false;

  // Create geometry
  this->set_eps0(eps0);
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  if(!turn_off_rod){
    RefCountedPtr<BaseIF> rod = RefCountedPtr<BaseIF> (new perlin_rod_if(R0,
									 c0,
									 c1,
									 false,
									 noise_amp,
									 f,
									 noise_pers,
									 noise_octav,
									 reseed));
    m_electrodes.push_back(electrode(rod, true, 1.0));
  }
}

rough_rod::~rough_rod(){
  
}
