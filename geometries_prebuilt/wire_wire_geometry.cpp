/*!
  @file wire_wire_geometry.cpp
  @brief Implementation of wire_wire_geometry.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "wire_wire_geometry.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

#include "perlin_sphere_if.H"
#include "box_if.H"
#include "new_sphere_if.H"

wire_wire_geometry::wire_wire_geometry(){

  Real eps0            = 1.0;
  Real wire1_radius    = 0.1;
  Real wire2_radius    = 0.2;
  Real wire1_potential = -0.5;
  Real wire2_potential = 0.5;

#if CH_SPACEDIM == 2
  RealVect wire1_center1 = RealVect(-1.0, 0.0);
  RealVect wire2_center1 = RealVect( 1.0, 0.0);
#elif CH_SPACEDIM == 3
  RealVect wire1_center1 = RealVect(-1.0, 0.0, -100);
  RealVect wire2_center1 = RealVect( 1.0, 0.0, -100);
  RealVect wire1_center2 = RealVect(-1.0, 0.0,  100);
  RealVect wire2_center2 = RealVect( 1.0, 0.0,  100);
#endif

  bool turn_off_wire1 = false;
  bool turn_off_wire2 = false;
  bool wire1_live     = true;
  bool wire2_live     = true;
  
  
  { // Get parameters from input script
    ParmParse pp("wire_wire_geometry");
    pp.query("eps0",         eps0);
    pp.query("wire1_radius", wire1_radius);
    pp.query("wire2_radius", wire2_radius);
    pp.query("wire1_potential", wire1_potential);
    pp.query("wire2_potential", wire2_potential);

    Vector<Real> vec(SpaceDim);
    std::string str;
    if(pp.contains("wire1_center1")){
      pp.getarr("wire1_center1", vec, 0, SpaceDim);
      wire1_center1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    if(pp.contains("wire2_center1")){
      pp.getarr("wire2_center1", vec, 0, SpaceDim);
      wire2_center1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
#if CH_SPACEDIM == 3
    if(pp.contains("wire1_center2")){
      pp.getarr("wire1_center2", vec, 0, SpaceDim);
      wire1_center2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    if(pp.contains("wire2_center1")){
      pp.getarr("wire2_center2", vec, 0, SpaceDim);
      wire2_center2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
#endif

    if(pp.contains("turn_off_wire1")){
      pp.get("turn_off_wire1", str);
      if(str == "true"){
	turn_off_wire1 = true;
      }
    }
    if(pp.contains("turn_off_wire2")){
      pp.get("turn_off_wire2", str);
      if(str == "true"){
	turn_off_wire2 = true;
      }
    }

    if(pp.contains("wire1_live")){
      pp.get("wire1_live", str);
      if(str == "true"){
	wire1_live = true;
      }
    }

    if(pp.contains("wire2_live")){
      pp.get("wire2_live", str);
      if(str == "true"){
	wire2_live = true;
      }
    }
  }

  // Create geometry



  // Create geometry
  this->set_eps0(eps0);
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  if(!turn_off_wire1){
#if CH_SPACEDIM == 2
    RefCountedPtr<BaseIF> wire1 = RefCountedPtr<BaseIF> (new new_sphere_if(wire1_center1, wire1_radius, false));
#elif CH_SPACEDIM == 3
    RefCountedPtr<BaseIF> wire1 = RefCountedPtr<BaseIF> (new new_sphere_if(wire1_center1, wire1_center2, wire1_radius, false));
#endif
    m_electrodes.push_back(electrode(wire1, wire1_live, wire1_potential));

  }
  
  if(!turn_off_wire2){
#if CH_SPACEDIM == 2
    RefCountedPtr<BaseIF> wire2 = RefCountedPtr<BaseIF> (new new_sphere_if(wire2_center1, wire2_radius, false));
#elif CH_SPACEDIM == 3
    RefCountedPtr<BaseIF> wire2 = RefCountedPtr<BaseIF> (new new_sphere_if(wire2_center1, wire2_center2, wire2_radius, false));
#endif
    m_electrodes.push_back(electrode(wire2, wire2_live, wire2_potential));
  }
}

wire_wire_geometry::~wire_wire_geometry(){
  
}
