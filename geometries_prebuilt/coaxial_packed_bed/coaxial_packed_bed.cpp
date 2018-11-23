/*!
  @file coaxial_packed_bed.cpp
  @brief Implementation of coaxial_packed_bed.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "coaxial_packed_bed.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>
#include <MultiSphereIF.H>

#include "cylinder_if.H"
#include "perlin_sphere_if.H"
#include "box_if.H"
#include "new_sphere_if.H"
#include "cylinder_if.H"

#define hard_spheres 1

coaxial_packed_bed::coaxial_packed_bed(){
  this->build_conductors();
  this->build_packed_bed();
}

coaxial_packed_bed::~coaxial_packed_bed(){
  
}

void coaxial_packed_bed::reseed(){
  int seed = time(NULL);   // Reseed the RNG

  // Haven't had trouble with this, but I'm putting this here for safety. 
#ifdef CH_MPI
  int result = MPI_Bcast(&seed, 1, MPI_INT, 0, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("coaxial_packed_bed::reseed() - broadcast failed");
  }
#endif
  std::srand(seed);
}

Real coaxial_packed_bed::random(){
  return  static_cast <Real> (rand()) / static_cast <Real> (RAND_MAX);
}

void coaxial_packed_bed::build_conductors(){
  m_rad_inner = 0.1;
  m_rad_outer = 0.9;
  
  const RealVect c0 = RealVect(D_DECL(0.0, 0.0, -2.0));
  const RealVect c1 = RealVect(D_DECL(0.0, 0.0,  2.0));

  std::string str          = "";
  bool turn_off_outer      = false;
  bool turn_off_inner      = false;
  bool inner_live          = true;
  bool outer_live          = false;
  
  ParmParse pp("coaxial_packed_bed");

  pp.get("outer_radius", m_rad_outer);
  pp.get("inner_radius", m_rad_inner);
  
  if(pp.contains("turn_off_outer")){
    pp.get("turn_off_outer", str);
    turn_off_outer = (str == "true") ? true : false;
  }
  if(pp.contains("turn_off_inner")){
    pp.get("turn_off_inner", str);
    turn_off_inner = (str == "true") ? true : false;
  }
  if(pp.contains("outer_live")){
    pp.get("outer_live", str);
    outer_live = (str == "true") ? true : false;
  }
  if(pp.contains("inner_live")){
    pp.get("inner_live", str);
    inner_live = (str == "true") ? true : false;
  }

  // Create geometry
  const Real eps0 = 1.0;
  this->set_eps0(eps0);

  m_electrodes.resize(0);
  if(!turn_off_outer){
    RefCountedPtr<BaseIF> outer = RefCountedPtr<BaseIF> (new cylinder_if(c0, c1, m_rad_outer, true));
    m_electrodes.push_back(electrode(outer, outer_live, 1.0));
  }
  if(!turn_off_inner){
    RefCountedPtr<BaseIF> inner = RefCountedPtr<BaseIF> (new cylinder_if(c0, c1, m_rad_inner, false));
    m_electrodes.push_back(electrode(inner, inner_live, 1.0));
  }

#if 0 // debug
  m_electrodes.resize(0);
  RefCountedPtr<BaseIF> sph = RefCountedPtr<BaseIF> (new new_sphere_if(RealVect::Zero, 0.2, false));
  m_electrodes.push_back(electrode(sph, true, 1.0));
#endif
}

void coaxial_packed_bed::build_packed_bed(){
  ParmParse pp("coaxial_packed_bed");
  
  // Init RNG
  std::string str;
  pp.get("reseed", str);
  if(str == "true"){
    this->reseed();
  }
  else{
    std::srand(0);
  }

  bool turn_off_spheres = false;
  if(pp.contains("turn_off_spheres")){
    pp.get("turn_off_spheres", str);
    turn_off_spheres = (str == "true") ? true : false;
  }

  if(turn_off_spheres){
    m_dielectrics.resize(0);
  }
  else{
    Real eps1 = 1.0;
    Real zlo, zhi;
    pp.get("num_spheres",         m_num_spheres);
    pp.get("sphere_radius",       m_rad_sphere);
    pp.get("sphere_min_dist",     m_min_dist);
    pp.get("sphere_permittivity", eps1);
    pp.get("z_low",               zlo);
    pp.get("z_high",              zhi);

    // Draw sphere centers
    int ndrawn = 0;
    Vector<RealVect> centers;
    Vector<Real> radii;
    while(ndrawn < m_num_spheres){
      RealVect c;
      for (int dir = 0; dir < SpaceDim; dir++){
	c[dir] = 2.0*m_rad_outer*(random()-0.5);
      }
#if CH_SPACEDIM==3
      c[2] = zlo + random()*(zhi-zlo);
#endif
      const RealVect rdist = RealVect(D_DECL(c[0], c[1], 0.0));
      const bool test_dist_inner = (rdist.vectorLength() - (m_rad_sphere+m_rad_inner) >= m_min_dist);     // True = OK!
      const bool test_dist_outer = (rdist.vectorLength() + m_rad_sphere <=  (m_rad_outer - m_min_dist));  // True = OK!
#if CH_SPACEDIM==3
      const bool test_dist_low   = (c[2] - zlo) >= m_min_dist + m_rad_sphere;
      const bool test_dist_high  = (zhi - c[2]) >= m_min_dist + m_rad_sphere;
#else
      const bool test_dist_low  = true;
      const bool test_dist_high = true;
#endif
      if(test_dist_inner && test_dist_outer && test_dist_low && test_dist_high){ // Does not intersect coaxial geometry

	bool isect_others = false;
	for (int i = 0; i < centers.size(); i++){
	  const Real dist = (c-centers[i]).vectorLength() - 2*m_rad_sphere;
	  if(dist < m_min_dist) isect_others = true;

	}

	if(!isect_others){
	  centers.push_back(c);
	  ndrawn++;
	}
      }
    }
    radii.resize(centers.size(), m_rad_sphere);

    RefCountedPtr<BaseIF> spheres = RefCountedPtr<BaseIF> (new MultiSphereIF(radii, centers, false));
    m_dielectrics.push_back(dielectric(spheres, eps1));
  }
}
