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

rough_sphere::rough_sphere(){

  Real eps, eps0, sphere_radius, noise_persistence, noise_amplitude;
  RealVect noise_frequency, sphere_center;
  int noise_octaves;
  bool noise_reseed, turn_off_sphere, live, use_new_gshop, inside;

  std::string str, which_material;
  Vector<Real> v0(SpaceDim);

  // Get a bunch of parameters from the input script
  ParmParse pp("rough_sphere");
  pp.get   ("eps0",              eps0);
  pp.get   ("permittivity",      eps);
  pp.get   ("sphere_radius",     sphere_radius);
  pp.get   ("noise_amplitude",   noise_amplitude);
  pp.get   ("noise_octaves",     noise_octaves);
  pp.get   ("noise_persistence", noise_persistence);
  pp.get   ("which_material",    which_material);
  pp.get   ("turn_off_sphere",   str);             turn_off_sphere = (str == "true") ? true : false;
  pp.get   ("live",              str);             live            = (str == "true") ? true : false;
  pp.get   ("noise_reseed",      str);             noise_reseed    = (str == "true") ? true : false;
  pp.get   ("use_new_gshop",     str);             use_new_gshop   = (str == "true") ? true : false;
  pp.getarr("sphere_center",     v0, 0, SpaceDim); sphere_center   = RealVect(D_DECL(v0[0], v0[1], v0[2]));
  pp.getarr("noise_frequency",   v0, 0, SpaceDim); noise_frequency = RealVect(D_DECL(v0[0], v0[1], v0[2]));

  set_eps0(eps0);
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  // Create geometry
  if(!turn_off_sphere){

    bool inside;
    if(which_material == "electrode"){
      inside = false;
    }
    else if(which_material == "dielectric"){
      inside = false;
    }
    else{
      MayDay::Abort("rough_sphere::rough_sphere - unknown parameter 'which_material'");
    }

    const real_box full_box(-1E10*RealVect::Unit, 1.E10*RealVect::Unit);
    const real_box sph_box(sphere_center - 1.25*(sphere_radius+noise_amplitude)*RealVect::Unit,
			   sphere_center + 1.25*(sphere_radius+noise_amplitude)*RealVect::Unit);

    RefCountedPtr<BaseIF> sph = RefCountedPtr<BaseIF> (new perlin_sphere_if(sphere_radius,
									    sphere_center,
									    inside,
									    noise_amplitude,
									    noise_frequency,
									    noise_persistence,
									    noise_octaves,
									    noise_reseed));

    if(which_material == "electrode"){
      m_electrodes.resize(1);
      m_electrodes[0].define(sph, live, 1.0);

      if(use_new_gshop){
	computational_geometry::s_use_new_gshop = true;
	
	m_regular_voxels_gas.push_back(full_box);
	m_bounded_voxels_gas.push_back(sph_box);
      }
    }
    else if(which_material == "dielectric"){
      m_dielectrics.resize(1);
      m_dielectrics[0].define(sph, eps);

      if(use_new_gshop){
	computational_geometry::s_use_new_gshop = true;

	m_regular_voxels_gas.push_back(full_box);
	m_bounded_voxels_gas.push_back(sph_box);
	
	m_covered_voxels_sol.push_back(full_box);
	m_bounded_voxels_sol.push_back(sph_box);
      }
    }
  }
}

rough_sphere::~rough_sphere(){
  
}
