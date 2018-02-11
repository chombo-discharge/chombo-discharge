/*!
  @file mechanical_shaft.cpp
  @brief Implementation of mechanical_shaft.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "mechanical_shaft.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

#include "perlin_sphere_if.H"
#include "box_if.H"
#include "new_sphere_if.H"
#include "polygon_rod_if.H"
#include "hollow_cylinder_if.H"

mechanical_shaft::mechanical_shaft(){
#if CH_SPACEDIM == 2
  MayDay::Abort("mechanical_shaft::mechanical_shaft - this class is for 3D only");
#endif

  m_electrodes.resize(0);
  m_dielectrics.resize(0);

  std::string str;
  Real eps0         = 1.0;
  Real corner_curv  = 250E-6;
  {
    ParmParse pp("mechanical_shaft");
    pp.query("eps0", eps0);
    pp.query("corner_curvaturse", corner_curv);
    this->set_eps0(eps0);
  }

  
  bool has_electrode           = false;
  bool electrode_live         = true;
  Real electrode_inner_rad    = 7.E-3;
  Real electrode_outer_rad    = 1.E-2;
  Real electrode_height       = 7.E-3;
  RealVect electrode_center   = RealVect(D_DECL(0., 0., 0.7));

  { // Get parameters for the electrode
    ParmParse pp("mechanical_shaft");
    Vector<Real> center;
    pp.query("electrode_inner_radius", electrode_inner_rad);
    pp.query("electrode_outer_radius", electrode_outer_rad);
    pp.query("electrode_height",       electrode_height);
    pp.query("turn_off_electrode",     str);
    if(str == "true"){
      has_electrode = false;
    }
    else if(str == "false"){
      has_electrode = true;
    }

    if(pp.contains("electrode_center")){
      pp.getarr("electrode_center", center, 0, SpaceDim);
      electrode_center = RealVect(D_DECL(center[0], center[1], center[2]));
    }

    pp.query("electrode_live", str);
    if(str == "true"){
      electrode_live = true;
    }
    else if(str == "false"){
      electrode_live = false;
    }
      
  }
  
  bool has_dielectric           = false;
  int dielectric_num_sides     = 6;       // Number of sides for the rod rod
  Real dielectric_radius       = 5.E-3;    // Rod radius
  Real dielectric_length       = 1;      // Rod length
  Real dielectric_permittivity = 4.0;     // Rod permittivity

  { // Get parameters for dielectric
    ParmParse pp("mechanical_shaft");
    pp.query("dielectric_rod_sides",    dielectric_num_sides);
    pp.query("dielectric_rod_radius",   dielectric_radius);
    pp.query("dielectric_rod_length",   dielectric_length);
    pp.query("dielectric_permittivity", dielectric_permittivity);
    pp.query("turn_off_dielectric",     str);
    if(str == "false"){
      has_dielectric = true;
    }
    else if(str == "true"){
      has_dielectric = false;
    }
  }
  

  if(has_electrode){ // Define electrode
    m_electrodes.resize(1);
    RefCountedPtr<BaseIF> electrode = RefCountedPtr<BaseIF> (new hollow_cylinder_if(electrode_inner_rad,
										    electrode_outer_rad,
										    electrode_height,
										    corner_curv,
										    electrode_center,
										    0));
    m_electrodes[0].define(electrode, electrode_live);    
  }

  
  if(has_dielectric){ // Define dielectric
    m_dielectrics.resize(1);
    RefCountedPtr<BaseIF> dielectric = RefCountedPtr<BaseIF> (new polygon_rod_if(dielectric_num_sides,
										 dielectric_radius,
										 dielectric_length,
										 corner_curv,
										 0));
    m_dielectrics[0].define(dielectric, dielectric_permittivity);

  }
}

mechanical_shaft::~mechanical_shaft(){
  
}
