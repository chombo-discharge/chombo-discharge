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

namespace ChomboDischarge {

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
    bool electrode_live;
    bool has_electrode;
    bool has_dielectric;
    Real dielectric_radius       = 5.E-3;   // Rod radius
    Real dielectric_length       = 1;       // Rod length
    Real dielectric_permittivity = 4.0;     // Rod permittivity
  
    pp.get("eps0",               eps0);
    pp.get("turn_on_dielectric", has_dielectric);
    pp.get("turn_on_electrode",  has_electrode);

  
    // Define shit. 
    if(has_electrode)  this->define_electrode();
    if(has_dielectric) this->define_dielectric();

    set_eps0(eps0);
  }

  mechanical_shaft::~mechanical_shaft(){
  
  }

  void mechanical_shaft::define_electrode(){
    ParmParse pp("mechanical_shaft.electrode");

    Vector<Real> vec(SpaceDim);
    bool live;
    Real innerRadius, outerRadius, curvature;
    RealVect c1, c2;

    pp.get("live",         live);
    pp.get("outer_radius", outerRadius);
    pp.get("inner_radius", innerRadius);
    pp.get("curvature",    curvature);

    pp.getarr("endpoint1", vec, 0, SpaceDim); c1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    pp.getarr("endpoint2", vec, 0, SpaceDim); c2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));

    RefCountedPtr<BaseIF> elec = RefCountedPtr<BaseIF> (new hollow_cylinder_if(c1, c2, outerRadius, innerRadius, curvature, false));

    m_electrodes.resize(1);
    m_electrodes[0].define(elec, live);
  }

  void mechanical_shaft::define_dielectric(){
    ParmParse pp("mechanical_shaft.dielectric");

    std::string str;
    Real eps;
  
    pp.get("shaft_shape", str);
    pp.get("permittivity", eps);

    // Get the shape
    RefCountedPtr<BaseIF> shaft = RefCountedPtr<BaseIF>(nullptr);
    if(str == "polygon"){
      shaft = this->get_polygon();
    }
    else if(str == "cylinder"){
      shaft = this->get_cylinder();
    }
    else if(str == "cyl_profile"){
      shaft = this->get_cylinder_profile();
    }
    else{
      MayDay::Abort("mechanical_shaft::define_dielectric - unknown argument 'shaft_shape' passed");
    }

    m_dielectrics.resize(1);
    m_dielectrics[0].define(shaft, eps);
  }

  RefCountedPtr<BaseIF> mechanical_shaft::get_polygon(){
    ParmParse pp("mechanical_shaft.dielectric.polygon");

    int numSides;
    RealVect c1, c2;
    Real radius, curv;

    Vector<Real> vec;

    pp.get("num_sides", numSides);
    pp.get("radius",    radius);
    pp.get("curvature", curv);

    pp.getarr("endpoint1", vec, 0, SpaceDim); c1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    pp.getarr("endpoint2", vec, 0, SpaceDim); c2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));

    return RefCountedPtr<BaseIF> (new polygon_rod_if(c1, c2, radius, curv, numSides, false));
  }

  RefCountedPtr<BaseIF> mechanical_shaft::get_cylinder(){
    ParmParse pp("mechanical_shaft.dielectric.polygon");
    
    int numSides;
    RealVect c1, c2;
    Real radius, length, curv;

    Vector<Real> vec;

    pp.get("radius",    radius);
    pp.get("curvature", curv);

    pp.getarr("endpoint1", vec, 0, SpaceDim); c1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    pp.getarr("endpoint2", vec, 0, SpaceDim); c2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));

    return RefCountedPtr<BaseIF> (new rounded_cylinder_if(c1, c2, radius, curv, false));
  }

  RefCountedPtr<BaseIF> mechanical_shaft::get_cylinder_profile(){
    ParmParse pp("mechanical_shaft.dielectric.cyl_profile");
    
    RealVect c1, c2;
    Real cylRad, torusMajor, torusMinor, ccDist, shift, curv, nLeft, nRight;
    Vector<Real> vec;

    pp.getarr("endpoint1", vec, 0, SpaceDim); c1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    pp.getarr("endpoint2", vec, 0, SpaceDim); c2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  
    pp.get("cylinder_radius", cylRad);
    pp.get("torus_major",     torusMajor);
    pp.get("torus_minor",     torusMinor);
    pp.get("torus_distance",  ccDist);
    pp.get("shift",           shift);
    pp.get("curvature",       curv);
    pp.get("nleft",           nLeft);
    pp.get("nright",          nRight);

    return RefCountedPtr<BaseIF> (new profile_cylinder_if(c1, c2, cylRad, torusMajor, torusMinor, ccDist, shift, curv, nLeft, nRight, false));
  }
}
