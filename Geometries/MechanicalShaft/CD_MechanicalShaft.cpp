/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_MechanicalShaft.cpp
  @brief  Implementation of CD_MechanicalShaft.H
  @author Robert Marskar
*/

// Std includes
#include <string>
#include <iostream>
#include <fstream>

// Chombo includes
#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

// Our includes
n#include <CD_MechanicalShaft.H>
#include <CD_PerlinSphereSdf.H>
#include <CD_BoxSdf.H>
#include <CD_CylinderSdf.H>
#include <CD_ProfileCylinderIF.H>
#include <CD_PolygonRodIF.H>
#include <CD_HollowCylinderIF.H>
#include <CD_RoundedCylinderIF.H>
#include <CD_NamespaceHeader.H>

MechanicalShaft::MechanicalShaft(){
#if CH_SPACEDIM == 2
  MayDay::Abort("MechanicalShaft::MechanicalShaft - this class is for 3D only");
#endif

  m_electrodes.resize(0);
  m_dielectrics.resize(0);

  std::string str;
  Vector<Real> vec;
  
  Real eps0;
  Real corner_curv;

  ParmParse pp("MechanicalShaft");


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
  if(has_electrode)  this->defineElectrode();
  if(has_dielectric) this->defineDielectric();

  setGasPermittivity(eps0);
}

MechanicalShaft::~MechanicalShaft(){
  
}

void MechanicalShaft::defineElectrode(){
  ParmParse pp("MechanicalShaft.electrode");

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

  RefCountedPtr<BaseIF> elec = RefCountedPtr<BaseIF> (new HollowCylinderIF(c1, c2, outerRadius, innerRadius, curvature, false));

  m_electrodes.resize(1);
  m_electrodes[0].define(elec, live);
}

void MechanicalShaft::defineDielectric(){
  ParmParse pp("MechanicalShaft.dielectric");

  std::string str;
  Real eps;
  
  pp.get("shaft_shape", str);
  pp.get("permittivity", eps);

  // Get the shape
  RefCountedPtr<BaseIF> shaft = RefCountedPtr<BaseIF>(nullptr);
  if(str == "polygon"){
    shaft = this->getPolygon();
  }
  else if(str == "cylinder"){
    shaft = this->getCylinder();
  }
  else if(str == "cyl_profile"){
    shaft = this->getCylinderProfile();
  }
  else{
    MayDay::Abort("MechanicalShaft::defineDielectric - unknown argument 'shaft_shape' passed");
  }

  m_dielectrics.resize(1);
  m_dielectrics[0].define(shaft, eps);
}

RefCountedPtr<BaseIF> MechanicalShaft::getPolygon(){
  ParmParse pp("MechanicalShaft.dielectric.polygon");

  int numSides;
  RealVect c1, c2;
  Real radius, curv;

  Vector<Real> vec;

  pp.get("num_sides", numSides);
  pp.get("radius",    radius);
  pp.get("curvature", curv);

  pp.getarr("endpoint1", vec, 0, SpaceDim); c1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("endpoint2", vec, 0, SpaceDim); c2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));

  return RefCountedPtr<BaseIF> (new PolygonRodIF(c1, c2, radius, curv, numSides, false));
}

RefCountedPtr<BaseIF> MechanicalShaft::getCylinder(){
  ParmParse pp("MechanicalShaft.dielectric.polygon");
    
  int numSides;
  RealVect c1, c2;
  Real radius, length, curv;

  Vector<Real> vec;

  pp.get("radius",    radius);
  pp.get("curvature", curv);

  pp.getarr("endpoint1", vec, 0, SpaceDim); c1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("endpoint2", vec, 0, SpaceDim); c2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));

  return RefCountedPtr<BaseIF> (new RoundedCylinderIF(c1, c2, radius, curv, false));
}

RefCountedPtr<BaseIF> MechanicalShaft::getCylinderProfile(){
  ParmParse pp("MechanicalShaft.dielectric.cylProfile");
    
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

  return RefCountedPtr<BaseIF> (new ProfileCylinderIF(c1, c2, cylRad, torusMajor, torusMinor, ccDist, shift, curv, nLeft, nRight, false));
}

#include <CD_NamespaceFooter.H>
