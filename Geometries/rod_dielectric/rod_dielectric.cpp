/*!
  @file   rod_dielectric.cpp
  @brief  Implementation of rod_dielectric.H
  @author Robert Marskar
  @date   Nov. 2017
*/

#include "rod_dielectric.H"
#include "rod_if.H"
#include "perlin_slab_if.H"
#include "new_sphere_if.H"
#include "wedge_if.H"
#include <CD_RoundedBoxIF.H>

#include <PlaneIF.H>
#include <ParmParse.H>

#include "CD_NamespaceHeader.H"

rod_dielectric::rod_dielectric(){
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  this->setGasPermittivity(1.0);

  ParmParse ppRod("rod_dielectric.electrode");
  ParmParse ppIns("rod_dielectric.dielectric");

  bool useRod;
  bool useIns;

  ppRod.get("on", useRod);
  ppIns.get("on", useIns);

  if(useRod) this->define_electrode();
  if(useIns) this->define_insulator();
}

rod_dielectric::~rod_dielectric(){
  
}

void rod_dielectric::define_electrode(){
  Vector<Real> v(SpaceDim);
  RealVect e1, e2;
  Real r;
  bool live;

  ParmParse pp("rod_dielectric.electrode");

  pp.get("radius", r);
  pp.get("live",   live);
  pp.getarr("endpoint1", v, 0, SpaceDim); e1 = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("endpoint2", v, 0, SpaceDim); e2 = RealVect(D_DECL(v[0], v[1], v[2]));

  RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF> (new rod_if(e1, e2, r, false));

  m_electrodes.push_back(electrode(bif, live));
}

void rod_dielectric::define_insulator(){
  ParmParse pp("rod_dielectric.dielectric");

  std::string str;
  Real eps;
  
  pp.get("shape", str);
  pp.get("permittivity", eps);

  RefCountedPtr<BaseIF> bif;
  if(str == "plane"){
    bif = this->get_plane();
  }
  else if(str == "box"){
    bif = this->get_box();
  }
  else if(str == "perlin_box"){
    bif = this->get_perlin_box();
  }
  else if(str == "sphere"){
    bif = this->get_sphere();
  }
  else{
    MayDay::Abort("rod_dielectric.:define_insulator - unsupported shape");
  }
  
  m_dielectrics.push_back(dielectric(bif, eps));
}



RefCountedPtr<BaseIF> rod_dielectric::get_box(){
  ParmParse pp("rod_dielectric.box");

  Vector<Real> v(SpaceDim);
  RealVect lo, hi;
  Real curv;

  pp.get("curvature", curv);
  pp.getarr("lo_corner", v, 0, SpaceDim); lo = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("hi_corner", v, 0, SpaceDim); hi = RealVect(D_DECL(v[0], v[1], v[2]));


  return RefCountedPtr<BaseIF> (new RoundedBoxIF(lo, hi, curv, false));
}

RefCountedPtr<BaseIF> rod_dielectric::get_plane(){
  ParmParse pp("rod_dielectric.plane");

  Vector<Real> v(SpaceDim);
  RealVect p, n;

  pp.getarr("point",  v, 0, SpaceDim); p    = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("normal", v, 0, SpaceDim); n    = RealVect(D_DECL(v[0], v[1], v[2]));

  return RefCountedPtr<BaseIF> (new PlaneIF(n, p, true));
}

RefCountedPtr<BaseIF> rod_dielectric::get_sphere(){
  ParmParse pp("rod_dielectric.sphere");

  Vector<Real> v(SpaceDim);
  RealVect p;
  Real r;

  pp.get("radius", r);
  pp.getarr("center", v, 0, SpaceDim); p = RealVect(D_DECL(v[0], v[1], v[2]));

  return RefCountedPtr<BaseIF> (new new_sphere_if(p, r, false));
}

RefCountedPtr<BaseIF> rod_dielectric::get_perlin_box(){
  ParmParse pp("rod_dielectric.perlin_box");
  
  Vector<Real> v(SpaceDim);
  RealVect p, n, xyz, freq;
  Real amp, persist, curv;
  int octaves;
  bool reseed;

  pp.get("curvature",     curv);
  pp.get("noise_octaves", octaves);
  pp.get("noise_amp",     amp);
  pp.get("noise_reseed",  reseed);
  pp.get("noise_persist", persist);

  pp.getarr("point",      v, 0, SpaceDim); p    = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("normal",     v, 0, SpaceDim); n    = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("dimensions", v, 0, SpaceDim); xyz  = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("noise_freq", v, 0, SpaceDim); freq = RealVect(D_DECL(v[0], v[1], v[2]));

  return RefCountedPtr<BaseIF> (new perlin_slab_if(p,n, xyz, freq, octaves, amp, persist, curv, reseed, false));
}
#include "CD_NamespaceFooter.H"
