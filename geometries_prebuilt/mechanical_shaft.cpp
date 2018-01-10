/*!
  @file mechanical_shaft.cpp
  @brief Implementation of mechanical_shaft.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "mechanical_shaft.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

#include "perlin_sphere_if.H"
#include "box_if.H"
#include "new_sphere_if.H"
#include "polygon_rod_if.H"
#include "hollow_cylinder_if.H"

mechanical_shaft::mechanical_shaft(){
  
  Real eps0         = 1.0;
  Real eps_mat      = 5.0;

  Real inner_rad    = 0.2;
  Real outer_rad    = 0.5;
  Real height       = 0.3;
  Real curv         = 0.05;
  RealVect center   = RealVect(D_DECL(0., 0., -1));

  int nsides  = 6;
  Real radius = 0.15;
  Real length = 10;
  
  bool live   = true;

  // Create geometry
  this->set_eps0(eps0);
  m_dielectrics.resize(1);
  m_electrodes.resize(1);

  RefCountedPtr<BaseIF> dielectric = RefCountedPtr<BaseIF> (new polygon_rod_if(nsides,
									       radius,
									       length,
									       curv,
									       0));

  RefCountedPtr<BaseIF> electrode = RefCountedPtr<BaseIF> (new hollow_cylinder_if(inner_rad,
										  outer_rad,
										  height,
										  curv,
										  center,
										  0));
  m_dielectrics[0].define(dielectric, eps_mat);
  m_electrodes[0].define(electrode, live);
}

mechanical_shaft::~mechanical_shaft(){
}
