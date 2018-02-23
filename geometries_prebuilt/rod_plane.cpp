/*!
  @file rod_plane_geometry.cpp
  @brief Implementation of rod_plane_geometry.H
  @author Robert Marskar
  @date Nov. 2017
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

  Real eps0        = 1.0;
  {
    ParmParse pp("rod_plane");
    pp.get("eps0",              eps0);
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
    ParmParse pp("rod_plane");
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


  
  // Dielectric plane
  bool has_dielectric      = true;
  Real plane_permittivity  = 5.0;
  RealVect plane_normal    = RealVect(BASISV(SpaceDim - 1));
  RealVect plane_point     = RealVect::Zero;

  RealVect noise_frequency = RealVect::Unit;  // Noise frequency
  Real noise_amplitude     = 0.0;             // Noise amplitude
  Real noise_persistence   = 0.5;             // Noise persistence
  int noise_octaves        = 1;               // Number of noise octaves
  bool noise_reseed        = false;           // Reseed
    

  {
    ParmParse pp("rod_plane");
    Vector<Real> vec(SpaceDim);
    std::string str;
    pp.get("turn_off_dielectric", str);
    if(str == "true"){
      has_dielectric = false;
    }
    else if(str == "false"){
      has_dielectric = true;
    }

    pp.query("dielectric_permittivity", plane_permittivity);
    if(pp.contains("plane_normal")){
      pp.getarr("plane_normal", vec, 0, SpaceDim);
      plane_normal = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    if(pp.contains("plane_point")){
      pp.getarr("plane_point", vec, 0, SpaceDim);
      plane_point = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    
    pp.query("noise_amplitude",   noise_amplitude);
    pp.query("noise_persistence", noise_persistence);
    pp.query("noise_octaves", noise_octaves);
    
    if(pp.contains("noise_frequency")){
      pp.getarr("noise_frequency", vec, 0, SpaceDim);
      noise_frequency = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    if(pp.contains("noise_reseed")){
      pp.get("noise_reseed", str);
      if(str == "true"){
	noise_reseed = true;
      }
      else if(str == "false"){
	noise_reseed = false;
      }
    }
  }

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

rod_plane::~rod_plane(){
  
}
