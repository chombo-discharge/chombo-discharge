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

  std::string str;
  Vector<Real> v(SpaceDim);
  Real eps0, R0, noise_pers, noise_amp;
  RealVect f, c0, c1;
  int noise_octav   = 1;
  bool reseed, turn_off_rod, live;
  
  ParmParse pp("rough_rod");
  pp.get   ("eps0",             eps0);
  pp.get   ("radius",            R0);
  pp.get   ("noise_octaves",     noise_octav);
  pp.get   ("noise_amplitude",   noise_amp);
  pp.get   ("noise_persistence", noise_pers);
  pp.get   ("turn_off_rod",      str);               turn_off_rod = (str == "true") ? true : false;
  pp.get   ("live",              str);                       live = (str == "true") ? true : false;
  pp.get   ("noise_reseed",      str);               reseed = (str == "true") ? true : false;
  pp.getarr("center1",           v, 0, SpaceDim); c0 = RealVect(D_DECL(v0[0], v0[1], v0[2]));
  pp.getarr("center2",           v, 0, SpaceDim); c1 = RealVect(D_DECL(v1[0], v1[1], v1[2]));
  pp.getarr("noise_frequency",   v, 0, SpaceDim); f  = RealVect(D_DECL(v2[0], v2[1], v2[2]));

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
