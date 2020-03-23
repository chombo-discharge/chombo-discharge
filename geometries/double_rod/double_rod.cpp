/*!
  @file   double_rod.cpp
  @brief  Implementation of double_rod.H
  @author Robert Marskar
  @date   Nov. 2017
*/

#include "double_rod.H"
#include "rod_if.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>

double_rod::double_rod(){
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  Real eps0        = 1.0;
  {
    ParmParse pp("double_rod");
    pp.get("eps0",              eps0);
  }
  this->set_eps0(eps0);

  
  
  // Electrode parameters
  bool turn_off_rod1    = false;
  bool turn_off_rod2    = false;
  bool rod1_live        = false;
  bool rod2_live        = false;
  
  Real rod1_radius      = 1.E-3;
  Real rod2_radius      = 1.E-3;
  
  RealVect rod1_center1 = RealVect::Zero;
  RealVect rod1_center2 = RealVect::Unit;
  RealVect rod2_center1 = RealVect::Zero;
  RealVect rod2_center2 = RealVect::Unit;

  { // Get parameterse for electrode rod
    ParmParse pp("double_rod");
    std::string str1, str2;
    Vector<Real> vec(SpaceDim);

    pp.get("turn_off_rod1", str1);
    pp.get("turn_off_rod2", str2);
    turn_off_rod1 = (str1 == "true") ? true : false;
    turn_off_rod2 = (str2 == "true") ? true : false;

    pp.get("rod1_live", str1);
    pp.get("rod2_live", str2);
    rod1_live = (str1 == "true") ? true : false;
    rod2_live = (str2 == "true") ? true : false;
    
    pp.get("rod1_radius", rod1_radius);
    pp.get("rod2_radius", rod2_radius);

    
    if(pp.contains("rod1_center1")){
      pp.getarr("rod1_center1", vec, 0, SpaceDim);
      rod1_center1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    if(pp.contains("rod1_center2")){
      pp.getarr("rod1_center2", vec, 0, SpaceDim);
      rod1_center2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }

    if(pp.contains("rod2_center1")){
      pp.getarr("rod2_center1", vec, 0, SpaceDim);
      rod2_center1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    if(pp.contains("rod2_center2")){
      pp.getarr("rod2_center2", vec, 0, SpaceDim);
      rod2_center2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
  }

  
  if(!turn_off_rod1){
    RefCountedPtr<BaseIF> rod = RefCountedPtr<BaseIF> (new rod_if(rod1_center1, rod1_center2, rod1_radius, false));
    electrode e = electrode(rod, rod1_live, 1.0);
    m_electrodes.push_back(e);
  }
  if(!turn_off_rod2){
    RefCountedPtr<BaseIF> rod = RefCountedPtr<BaseIF> (new rod_if(rod2_center1, rod2_center2, rod2_radius, false));
    electrode e = electrode(rod, rod2_live, 1.0);
    m_electrodes.push_back(e);
  }
}

double_rod::~double_rod(){
  
}
