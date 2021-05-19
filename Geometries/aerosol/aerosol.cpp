/*!
  @file aerosol.cpp
  @brief Implementation of aerosol.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "aerosol.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>

#include "perlin_sphere_if.H"

#include "CD_NamespaceHeader.H"

aerosol::aerosol(){

  Real eps, eps0, noise_persistence, noise_amplitude;
  RealVect noise_frequency;
  int noise_octaves, num_spheres;
  bool noise_reseed;

  Vector<Real> v(SpaceDim);

  // Get a bunch of parameters from the input script
  ParmParse pp("aerosol");
  
  pp.get   ("eps0",              eps0);
  pp.get   ("permittivity",      eps);
  pp.get   ("num_spheres",       num_spheres);
  pp.get   ("noise_amplitude",   noise_amplitude);
  pp.get   ("noise_octaves",     noise_octaves);
  pp.get   ("noise_persistence", noise_persistence);
  pp.get   ("noise_reseed",      noise_reseed);
  pp.getarr("noise_frequency",   v, 0, SpaceDim);  noise_frequency = RealVect(D_DECL(v[0], v[1], v[2]));


  setGasPermittivity(eps0);
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  // Add spheres
  for (int i = 0; i < num_spheres; i++){
    const int ndigits = (int) log10((double) num_spheres) + 1;

    RealVect center;
    Real radius;

    char* cstr = new char[ndigits];
    sprintf(cstr, "%d", 1+i);
    const std::string str = "aerosol.sphere" + std::string(cstr);

    ParmParse pp2(str);
    pp2.get   ("radius", radius);
    pp2.getarr("center", v, 0, SpaceDim); center = RealVect(D_DECL(v[0], v[1], v[2]));

    // Get the perlin sphere. 
    RefCountedPtr<BaseIF> sph = RefCountedPtr<BaseIF> (new perlin_sphere_if(radius,
									    center,
									    false,
									    noise_amplitude,
									    noise_frequency,
									    noise_persistence,
									    noise_octaves,
									    noise_reseed));

    // Add dielectric sphere
    m_dielectrics.push_back(Dielectric(sph, eps));
  }
}

aerosol::~aerosol(){
  
}
#include "CD_NamespaceFooter.H"
