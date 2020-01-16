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

aerosol::aerosol(){

  Real eps, eps0, sphere_radius, noise_persistence, noise_amplitude;
  RealVect noise_frequency, sphere_center;
  int noise_octaves;
  int num_spheres;
  bool noise_reseed, turn_off_sphere, live, use_new_gshop, inside;

  std::string str, which_material;
  Vector<Real> v0(SpaceDim);

  // Get a bunch of parameters from the input script
  ParmParse pp("aerosol");
  pp.get   ("eps0",              eps0);
  pp.get   ("permittivity",      eps);
  pp.get   ("num_spheres",       num_spheres);

  pp.get   ("noise_amplitude",   noise_amplitude);
  pp.get   ("noise_octaves",     noise_octaves);
  pp.get   ("noise_persistence", noise_persistence);
  pp.get   ("noise_reseed",      str);              noise_reseed    = (str == "true") ? true : false;
  pp.get   ("use_new_gshop",     str);              use_new_gshop   = (str == "true") ? true : false;
  pp.getarr("noise_frequency",   v0, 0, SpaceDim);  noise_frequency = RealVect(D_DECL(v0[0], v0[1], v0[2]));


  //
  set_eps0(eps0);
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  

  // Create geometry
  if(num_spheres > 0){
    const real_box full_box(-1E10*RealVect::Unit, 1.E10*RealVect::Unit);
    const int ndigits = (int) log10((double) num_spheres) + 1;

    if(use_new_gshop == true){
      computational_geometry::s_use_new_gshop = true;
      m_regular_voxels_gas.push_back(full_box);
      m_covered_voxels_sol.push_back(full_box);
    }

    // Add spheres as we go
    for (int i = 0; i < num_spheres; i++){

      char* cstr = new char[ndigits];
      sprintf(cstr, "%d", 1+i);

      const std::string str1 = "sphere" + std::string(cstr) + "_center";
      const std::string str2 = "sphere" + std::string(cstr) + "_radius";

      pp.get   (str2.c_str(), sphere_radius);
      pp.getarr(str1.c_str(), v0, 0, SpaceDim);   sphere_center   = RealVect(D_DECL(v0[0], v0[1], v0[2]));

      // sphere box for voxelized load balancing
      const real_box sph_box(sphere_center - 1.25*(sphere_radius+noise_amplitude)*RealVect::Unit,
			     sphere_center + 1.25*(sphere_radius+noise_amplitude)*RealVect::Unit);

      // actual IF object
      RefCountedPtr<BaseIF> sph = RefCountedPtr<BaseIF> (new perlin_sphere_if(sphere_radius,
									      sphere_center,
									      false,
									      noise_amplitude,
									      noise_frequency,
									      noise_persistence,
									      noise_octaves,
									      noise_reseed));

      // Add dielectric sphere
      m_dielectrics.push_back(dielectric(sph, eps));

      if(use_new_gshop){
	m_bounded_voxels_gas.push_back(sph_box);
	m_bounded_voxels_sol.push_back(sph_box);
      }
    }
  }
}

aerosol::~aerosol(){
  
}
