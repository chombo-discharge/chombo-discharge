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
#include "cylinder_if.H"
#include "profile_cylinder_if.H"
#include "polygon_rod_if.H"
#include "hollow_cylinder_if.H"
#include "rounded_cylinder_if.H"

mechanical_shaft::mechanical_shaft(){
#if CH_SPACEDIM == 2
  MayDay::Abort("mechanical_shaft::mechanical_shaft - this class is for 3D only");
#endif

  m_electrodes.resize(0);
  m_dielectrics.resize(0);

  std::string str;
  Vector<Real> vec;
  
  Real eps0;
  Real corner_curv;

  ParmParse pp("mechanical_shaft");


  std::string shape;
  int dielectric_num_sides     = 6;       // Number of sides for the rod rod
  bool electrode_live          = true;
  bool has_electrode           = false;
  bool has_dielectric          = false;
  Real dielectric_radius       = 5.E-3;   // Rod radius
  Real dielectric_length       = 1;       // Rod length
  Real dielectric_permittivity = 4.0;     // Rod permittivity
  Real electrode_inner_rad    = 7.E-3;
  Real electrode_outer_rad    = 1.E-2;
  Real electrode_height       = 7.E-3;
  RealVect electrode_center   = RealVect(D_DECL(0., 0., 0.7));

  
  pp.get("eps0", eps0);
  pp.get("shaft_shape", shape);
  pp.get("corner_curvatures", corner_curv);
  pp.get("electrode_inner_radius", electrode_inner_rad);
  pp.get("electrode_outer_radius", electrode_outer_rad);
  pp.get("electrode_height",       electrode_height);
  pp.get("dielectric_rod_sides",    dielectric_num_sides);
  pp.get("dielectric_rod_radius",   dielectric_radius);
  pp.get("dielectric_rod_length",   dielectric_length);
  pp.get("dielectric_permittivity", dielectric_permittivity);

  pp.get("turn_off_dielectric",     str);           has_dielectric   = (str == "true") ? false : true;
  pp.get("turn_off_electrode",  str);               has_electrode    = (str == "true") ? false : true;
  pp.get("electrode_live",      str);               electrode_live   = (str == "true") ? true : false;
  pp.getarr("electrode_center", vec, 0, SpaceDim);  electrode_center = RealVect(D_DECL(vec[0], vec[1], vec[2]));

  

  if(has_electrode){ // Define electrode
    m_electrodes.resize(1);
#if 0 // Original code
    RefCountedPtr<BaseIF> electrode = RefCountedPtr<BaseIF> (new hollow_cylinder_if(electrode_inner_rad,
										    electrode_outer_rad,
										    electrode_height,
										    corner_curv,
										    electrode_center,
										    0));
#else // Test code
    const Real radius = 1E-3;
    const Real curv = 0.0;
    const RealVect center1(-4E-2, -4E-2, 0.);
    const RealVect center2(4E-2, 4E-2, 4E-2);
    RefCountedPtr<BaseIF> electrode = RefCountedPtr<BaseIF> (new rounded_cylinder_if(center1, center2, radius, curv, false));
#endif
    m_electrodes[0].define(electrode, electrode_live);
  }

  
  if(has_dielectric){ // Define dielectric
    m_dielectrics.resize(1);
    RefCountedPtr<BaseIF> shaft;
    
    if(shape == "polygon"){
      shaft = RefCountedPtr<BaseIF> (new polygon_rod_if(dielectric_num_sides,
							dielectric_radius,
							dielectric_length,
							corner_curv,
							0));
    }
    else if(shape == "cylinder"){
      const RealVect zhat = RealVect(BASISV(2));
      shaft = RefCountedPtr<BaseIF> (new cylinder_if(electrode_center - 0.5*dielectric_length*zhat,
						     electrode_center + 0.5*dielectric_length*zhat,
						     dielectric_radius,
						     false));
						     
    }
    else if(shape == "cyl_profile"){
      int nleft, nright;
      Real rad, offset, shift, dist, curv;
      pp.get("cylprofile_nleft", nleft);
      pp.get("cylprofile_nright", nright);
      pp.get("cylprofile_rad", rad);
      pp.get("cylprofile_offset", offset);
      pp.get("cylprofile_shift", shift);
      pp.get("cylprofile_dist", dist);
      pp.get("cylprofile_curv", curv);
      shaft = RefCountedPtr<BaseIF> (new profile_cylinder_if(electrode_center,
							     dielectric_length,
							     dielectric_radius,
							     nleft,
							     nright,
							     rad,
							     offset,
							     shift,
							     dist, 
							     curv,
							     false));
    }


    m_dielectrics[0].define(shaft, dielectric_permittivity);
  }

  set_eps0(eps0);
}

mechanical_shaft::~mechanical_shaft(){
  
}
