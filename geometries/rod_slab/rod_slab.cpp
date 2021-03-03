/*!
  @file   rod_slab.cpp
  @brief  Implementation of rod_slab.H
  @author Robert Marskar
  @date   Nov. 2017
*/

#include "rod_slab.H"
#include "perlin_slab_if.H"
#include "rod_if.H"

#include <ParmParse.H>

rod_slab::rod_slab(){
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  this->set_eps0(1.0);

  ParmParse ppRod("rod_slab.electrode");
  ParmParse ppSlab("rod_slab.dielectric");

  bool useRod;
  bool useSlab;

  ppRod.get ("on", useRod);
  ppSlab.get("on", useSlab);

  if(useRod){ // Set up electrode
    Vector<Real> v(SpaceDim);
    RealVect e1, e2;
    Real r;
    bool live;

    ppRod.get("radius", r);
    ppRod.get("live",   live);
    ppRod.getarr("endpoint1", v, 0, SpaceDim); e1 = RealVect(D_DECL(v[0], v[1], v[2]));
    ppRod.getarr("endpoint2", v, 0, SpaceDim); e2 = RealVect(D_DECL(v[0], v[1], v[2]));

    RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF> (new rod_if(e1, e2, r, false));

    m_electrodes.push_back(electrode(bif, live));
  }

  if(useSlab){ // Set up dielectric
    Vector<Real> v(SpaceDim);
    RealVect p, n, xyz, freq;
    Real eps, amp, persist, curv;
    int octaves;
    bool reseed;

    ppSlab.get("eps",           eps);
    ppSlab.get("curvature",     curv);
    ppSlab.get("noise_octaves", octaves);
    ppSlab.get("noise_amp",     amp);
    ppSlab.get("noise_reseed",  reseed);
    ppSlab.get("noise_persist", persist);

    ppSlab.getarr("point",      v, 0, SpaceDim); p = RealVect(D_DECL(v[0], v[1], v[2]));
    ppSlab.getarr("normal",     v, 0, SpaceDim); n = RealVect(D_DECL(v[0], v[1], v[2]));
    ppSlab.getarr("dimensions", v, 0, SpaceDim); xyz = RealVect(D_DECL(v[0], v[1], v[2]));
    ppSlab.getarr("noise_freq", v, 0, SpaceDim); freq = RealVect(D_DECL(v[0], v[1], v[2]));

    RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF> (new perlin_slab_if(p,n, xyz, freq, octaves, amp, persist, curv, reseed, false));

    m_dielectrics.push_back(dielectric(bif, eps));
  }
}

rod_slab::~rod_slab(){
  
}
