/*!
  @file rough_sphere.cpp
  @brief Implementation of rough_sphere.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "rough_sphere.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>

#include "perlin_sphere_if.H"

#include "CD_NamespaceHeader.H"

rough_sphere::rough_sphere(){
  this->setGasPermittivity(1.0);

  ParmParse pp("rough_sphere");

  bool useSphere;
  std::string whichMaterial;

  pp.get("on",       useSphere);
  pp.get("material", whichMaterial);

  if(useSphere){
    bool reseed, live;
    RealVect f, c;
    Real r, persist, amp, eps;
    int octaves;
    Vector<Real> v;
    

    pp.get("radius", r);
    pp.get("live", live);
    pp.get("eps",  eps);
    pp.get("noise_reseed", reseed);
    pp.get("noise_amplitude", amp);
    pp.get("noise_octaves", octaves);
    pp.get("noise_persistence", persist);

    pp.getarr("noise_frequency", v, 0 ,SpaceDim); f  = RealVect(D_DECL(v[0], v[1], v[2]));
    pp.getarr("center",          v, 0 ,SpaceDim); c = RealVect(D_DECL(v[0], v[1], v[2]));

    RefCountedPtr<BaseIF> sph = RefCountedPtr<BaseIF> (new perlin_sphere_if(r, c, false, amp, f, persist, octaves, reseed));

    if      (whichMaterial == "electrode")  m_electrodes.push_back(Electrode(sph, live));
    else if (whichMaterial == "dielectric") m_dielectrics.push_back(Dielectric(sph, eps));
    else    MayDay::Abort("rough_sphere::rough_sphere - unknown material requested");
  }


}

rough_sphere::~rough_sphere(){
  
}
#include "CD_NamespaceFooter.H"
