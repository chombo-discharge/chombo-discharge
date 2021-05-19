/*!
  @file   rough_rod.cpp
  @brief  Implementation of rough_rod.H
  @author Robert Marskar
  @date   Nov. 2017
*/

#include "rough_rod.H"
#include "perlin_rod_if.H"

#include <ParmParse.H>

#include "CD_NamespaceHeader.H"

rough_rod::rough_rod(){
  this->setGasPermittivity(1.0);

  ParmParse pp("rough_rod");

  bool useRod;
  
  pp.get("on", useRod);

  if(useRod){
    bool live, reseed;
    RealVect f, e1, e2;
    Real r, persist, amp;
    int octaves;
    Vector<Real> v;

    pp.get("radius", r);
    pp.get("live", live);
    pp.get("noise_reseed", reseed);
    pp.get("noise_amplitude", amp);
    pp.get("noise_octaves", octaves);
    pp.get("noise_persistence", persist);

    pp.getarr("noise_frequency", v, 0 ,SpaceDim); f  = RealVect(D_DECL(v[0], v[1], v[2]));
    pp.getarr("endpoint1",       v, 0 ,SpaceDim); e1 = RealVect(D_DECL(v[0], v[1], v[2]));
    pp.getarr("endpoint2",       v, 0 ,SpaceDim); e2 = RealVect(D_DECL(v[0], v[1], v[2]));

    RefCountedPtr<BaseIF> rod = RefCountedPtr<BaseIF> (new perlin_rod_if(r, e1, e2, false, amp, f, persist, octaves, reseed));

    m_electrodes.push_back(electrode(rod, live));
  }
}

rough_rod::~rough_rod(){
  
}
#include "CD_NamespaceFooter.H"
