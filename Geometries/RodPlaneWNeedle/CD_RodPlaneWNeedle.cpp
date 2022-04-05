/* chombo-discharge
 * Copyright © 2022 NTNU.
 * Copyright © 2022 Fanny Skirbekk. 
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RodPlaneWNeedle.cpp
  @brief  Implementation of CD_RodPlaneWNeedle.H 
  @author Fanny Skirbekk
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_RodPlaneWNeedle.H>
#include <CD_RodIF.H>
#include <CD_NamespaceHeader.H>

RodPlaneWNeedle::RodPlaneWNeedle(){
  this->setGasPermittivity(1.0);
  
  ParmParse ppRod("RodPlaneWNeedle.rod");
  ParmParse ppNeedle("RodPlaneWNeedle.needle");
  ParmParse ppPlane("RodPlaneWNeedle.plane");

  bool useRod;
  bool useNeedle;
  bool usePlane;

  ppRod.get("on", useRod);
  ppNeedle.get("on", useNeedle);
  ppPlane.get("on", usePlane);

  if(useRod && useNeedle) this->defineRodWNeedle();
  else if(useRod)         this ->defineRod();
  else if(useNeedle)      this->defineNeedle();
  if(usePlane)            this->definePlane();
}

RodPlaneWNeedle::~RodPlaneWNeedle(){

}

void RodPlaneWNeedle::defineRodWNeedle(){
  Vector<Real> v(SpaceDim); 
  Real rodRadius;
  RealVect rodE1, rodE2;
  bool rodLive;

  Real needleRadius;
  Real tipRadius;
  RealVect needleE1, needleE2;
  bool needleLive;

  Real cornerCurve;

  ParmParse ppRod("RodPlaneWNeedle.rod");
  ParmParse ppNeedle("RodPlaneWNeedle.needle");

  ppRod.get("radius", rodRadius);
  ppRod.get("live", rodLive);

  ppNeedle.get("radius", needleRadius);
  ppNeedle.get("tipRadius", tipRadius);
  ppNeedle.get("live", needleLive);

  ppNeedle.get("cornerCurve", cornerCurve);

  ppRod.getarr("endpoint1", v, 0, SpaceDim); rodE1 = RealVect(D_DECL(v[0], v[1], v[2])); 
  ppRod.getarr("endpoint2", v, 0, SpaceDim); rodE2 = RealVect(D_DECL(v[0], v[1], v[2]));

  ppNeedle.getarr("endpointtip", v, 0, SpaceDim); needleE1 = RealVect(D_DECL(v[0], v[1], v[2])); 
  ppNeedle.getarr("endpointback", v, 0, SpaceDim); needleE2 = RealVect(D_DECL(v[0], v[1], v[2]));

  Vector<BaseIF*> electrodeParts;
  electrodeParts.push_back((BaseIF*) new RodIF(e1, e2, r, false));
  //add needle


  RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF> (new SmoothUnion(electrodeParts, cornerCurve)); 

  m_electrodes.push_back(Electrode(bif, (rodLive&&needleLive)));

  for(int i = 0; i < electrodeParts.size(); ++i){
    delete electrodeParts[i];
  }
}

void RodPlaneWNeedle::defineRod(){
  Vector<Real> v(SpaceDim); 
  Real radius;
  RealVect e1, e2;
  bool live;

  ParmParse pp("RodPlaneWNeedle.rod");

  pp.get("radius", radius);
  pp.get("live", live);

  pp.getarr("endpoint1", v, 0, SpaceDim); e1 = RealVect(D_DECL(v[0], v[1], v[2])); 
  pp.getarr("endpoint2", v, 0, SpaceDim); e2 = RealVect(D_DECL(v[0], v[1], v[2]));

  RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF> (new RodIF(e1, e2, r, false)); 

  m_electrodes.push_back(Electrode(bif, live));
}

void RodPlaneWNeedle::defineNeedle(){
  Vector<Real> v(SpaceDim); 
  Real radius;
  Real tipRadius;
  RealVect e1, e2;
  bool live;

  ParmParse pp("RodPlaneWNeedle.needle");

  pp.get("radius", radius);
  pp.get("tipRadius", tipRadius);
  pp.get("live", live);

  pp.getarr("endpointtip", v, 0, SpaceDim); e1 = RealVect(D_DECL(v[0], v[1], v[2])); 
  pp.getarr("endpointback", v, 0, SpaceDim); e2 = RealVect(D_DECL(v[0], v[1], v[2]));

  // create needle

  // RefCountedPtr<BaseIF> bif =

  // m_electrodes.push_back(Electrode(bif, live));
}

void RodPlaneWNeedle::definePlane(){
  Vector<Real> v(SpaceDim);
  Real radius;
  Real thickness;
  Real curve;
  RealVect normal, point;
  bool live;
  
  ParmParse pp("RodPlaneWNeedle.plane");

  pp.get("radius", radius);
  pp.get("thickness", thickness);
  pp.get("curve", curve);
  pp.get("live", live);
  pp.getarr("normal", v, 0, SpaceDim); normal = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("point", v, 0, SpaceDim); point = RealVect(D_DECL(v[0], v[1], v[2]));

  RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF> (new RoundedCylinderIF(point, point - thickness*normal, radius, curve, false)); 

  m_electrodes.push_back(Electrode(bif, live));
}

#include <CD_NamespaceFooter.H>
