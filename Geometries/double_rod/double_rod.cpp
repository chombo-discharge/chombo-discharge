/*!
  @file   double_rod.cpp
  @brief  Implementation of double_rod.H
  @author Robert Marskar
  @date   Nov. 2017
*/

#include "double_rod.H"
#include "rod_if.H"

#include <ParmParse.H>

#include "CD_NamespaceHeader.H"

double_rod::double_rod(){
  this->set_eps0(1.0);

  ParmParse pp1("double_rod.rod1");
  ParmParse pp2("double_rod.rod2");

  bool use_rod1;
  bool use_rod2;

  pp1.get("on", use_rod1);
  pp2.get("on", use_rod2);

  Vector<Real> v(SpaceDim);
  Real radius;
  RealVect e1, e2;
  bool live;
  
  if(use_rod1){
    pp1.get("radius", radius);
    pp1.get("live",   live);
    pp1.getarr("endpoint1", v, 0, SpaceDim); e1 = RealVect(D_DECL(v[0], v[1], v[2]));
    pp1.getarr("endpoint2", v, 0, SpaceDim); e2 = RealVect(D_DECL(v[0], v[1], v[2]));

    RefCountedPtr<BaseIF> rod1 = RefCountedPtr<BaseIF> (new rod_if(e1, e2, radius, false));

    m_electrodes.push_back(electrode(rod1, live));
  }

  if(use_rod2){
    pp2.get("radius", radius);
    pp2.get("live",   live);
    pp2.getarr("endpoint1", v, 0, SpaceDim); e1 = RealVect(D_DECL(v[0], v[1], v[2]));
    pp2.getarr("endpoint2", v, 0, SpaceDim); e2 = RealVect(D_DECL(v[0], v[1], v[2]));

    RefCountedPtr<BaseIF> rod1 = RefCountedPtr<BaseIF> (new rod_if(e1, e2, radius, false));

    m_electrodes.push_back(electrode(rod1, live));
  }
}

double_rod::~double_rod(){
  
}
#include "CD_NamespaceFooter.H"
