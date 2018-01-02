/*!
  @file sphere_sphere_geometry.cpp
  @brief Implementation of sphere_sphere_geometry.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "sphere_sphere_geometry.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

#include "perlin_sphere_if.H"
#include "box_if.H"
#include "new_sphere_if.H"

sphere_sphere_geometry::sphere_sphere_geometry(){

  Real eps0         = 1.0;
  Real eps_mat      = 10.0;
  Real elec_rad     = 0.15;
  Real diel_rad     = 0.15;
  Real elec_noise   = 3.E-3;
  Real diel_noise   = 3.E-3;
  Real elec_persist = 0.5;
  Real diel_persist = 0.5;
  
  RealVect elec_center     = -0.25*RealVect::Unit;
  RealVect diel_center     = 0.25*RealVect::Unit;
  RealVect elec_noise_freq = 400.*RealVect::Unit;
  RealVect diel_noise_freq = 400.*RealVect::Unit;
  
  bool live   = true;
  bool reseed = true;

  int elec_octaves = 1;
  int diel_octaves = 1;


  // ParmParse pp("sphere_sphere_geometry");
  // pp.get("electrode_radius",  elec_rad);
  // pp.getarr("electrode_center",  center, 0, SpaceDim); elec_center = RealVect(D_DECL(center[0], center[1], center[2]));
  // pp.get("electrode_live",    live); 
  // pp.get("eps0",              eps0);
  // pp.get("dielectric_radius", diel_rad);
  // pp.getarr("dielectric_center", center, 0, SpaceDim); diel_center = RealVect(D_DECL(center[0], center[1], center[2]));
  // pp.get("dielectric_eps",    eps_mat);


  // Create geometry
  this->set_eps0(eps0);
  m_dielectrics.resize(1);
  m_electrodes.resize(1);

  //  RefCountedPtr<BaseIF> diel_sphere = RefCountedPtr<BaseIF> (new SphereIF(diel_rad, diel_center, 0));
  RefCountedPtr<BaseIF> diel_sphere = RefCountedPtr<BaseIF> (new perlin_sphere_if(diel_rad,
  										  diel_center,
  										  0,
  										  diel_noise,
  										  diel_noise_freq,
  										  diel_persist,
  										  diel_octaves,
  										  false));
  RefCountedPtr<BaseIF> elec_sphere = RefCountedPtr<BaseIF> (new perlin_sphere_if(elec_rad,
  										  elec_center,
  										  0,
  										  elec_noise,
  										  elec_noise_freq,
  										  elec_persist,
  										  elec_octaves,
  										  false));
  m_dielectrics[0].define(diel_sphere, eps_mat);
  m_electrodes[0].define(elec_sphere, live);
}

//
sphere_sphere_geometry::~sphere_sphere_geometry(){
}

//
void sphere_sphere_geometry::dumpScript(std::string a_filename){

#ifdef CH_MPI
  if(procID() == 0){
#endif
  std::ofstream outfile(a_filename);

  outfile << "# SPHERE_SPHERE_GEOMETRY (sphere_sphere_geometry.H)\n";
  outfile << "# -------------------------------------------------\n";
  outfile << "# Template script for defining a prebuilt geometry\n";
  outfile << "# consisting of a dielectric sphere and an electrode\n";
  outfile << "# sphere\n\n";

  outfile << "sphere_sphere_geometry.electrode_radius  = \t\t # Electrode radius (float)\n";
  outfile << "sphere_sphere_geometry.electrode_center  = \t\t # Electrode sphere center (spacedim array of floats)\n";
  outfile << "sphere_sphere_geometry.electrode_live    = \t\t # Electrode is live or not (boolean) \n";
  outfile << "sphere_sphere_geometry.eps0              = \t\t # Gas permittivity\n";
  outfile << "sphere_sphere_geometry.dielectric_radius = \t\t # Dielectric radius (float) \n";
  outfile << "sphere_sphere_geometry.dielectric_center = \t\t # Dielectric center (spacedim array of floats)\n";
  outfile << "sphere_sphere_geometry.dielectric_eps    = \t\t # Dielectric permittivity (float) \n";

  outfile.close();

#ifdef CH_MPI
  }
#endif
}
