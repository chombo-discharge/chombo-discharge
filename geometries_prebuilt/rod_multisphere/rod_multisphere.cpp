/*!
  @file   rod_multisphere.cpp
  @brief  Implementation of rod_multisphere.H
  @author Robert Marskar
  @date   Aug. 2018
*/

#include "rod_multisphere.H"

#include <ParmParse.H>
#include <new_sphere_if.H>
#include <rod_if.H>
#include <math.h>


rod_multisphere::rod_multisphere(){
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  { // Get gas phase permittivity
    Real eps0 = 1.0;
    ParmParse pp("rod_multisphere");
    pp.query("eps0", eps0);
    this->set_eps0(eps0);
  }


  { // Define electrode
    bool live = true;
    bool has_electrode = true;
    Real elec_radius   = 1.0;
    RealVect center1 = RealVect(D_DECL(0.0, 0.0, 0.0));
    RealVect center2 = RealVect(D_DECL(0.0, 0.0, 0.0));

    ParmParse pp("rod_multisphere");

    std::string str;
    Vector<Real> vec(SpaceDim);
    pp.query("electrode_radius", elec_radius);
    if(pp.contains("electrode_center1")){
      pp.getarr("electrode_center1", vec, 0, SpaceDim);
      center1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    if(pp.contains("electrode_center2")){
      pp.getarr("electrode_center2", vec, 0, SpaceDim);
      center2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
    }
    pp.query("turn_off_electrode", str);
    if(str == "true"){
      has_electrode = false;
    }
    else if(str == "false"){
      has_electrode = true;
    }
    pp.query("electrode_live", str);
    if(str == "true"){
      live = true;
    }
    else if(str == "false"){
      live = false;
    }

    // Define electrode
    if(has_electrode){
      m_electrodes.resize(1);
      RefCountedPtr<BaseIF> rod = RefCountedPtr<BaseIF> (new rod_if(center1, center2, elec_radius, false));
      m_electrodes[0].define(rod, live);
    }
  }

  { // Define dielectrics
    int num_spheres = 0;
    ParmParse pp("rod_multisphere");
    pp.get("num_spheres", num_spheres);

    if(num_spheres > 0){

      m_dielectrics.resize(num_spheres);
      for (int i = 0; i < num_spheres; i++){
	Real radius;
	Real permittivity;
	Vector<Real> center;

	// Create a sphereX string
	const int ndigits = (int) log10((double) num_spheres) + 1;
	
	char* str_rad = new char[ndigits];
	char* str_eps = new char[ndigits];
	char* str_cen = new char[ndigits];
	
	sprintf(str_rad, "%s%d%s", "sphere", i, "_radius");
	sprintf(str_cen, "%s%d%s", "sphere", i, "_center");
	sprintf(str_eps, "%s%d%s", "sphere", i, "_eps");

	pp.get(str_rad, radius);
	pp.get(str_eps, permittivity);
	pp.getarr(str_cen, center, 0, SpaceDim);

	const RealVect cen = RealVect(D_DECL(center[0], center[1], center[2]));

	RefCountedPtr<BaseIF> sphere = RefCountedPtr<BaseIF> (new new_sphere_if(cen, radius, false));
	m_dielectrics[i].define(sphere, permittivity);
      }
    }
  }
}

rod_multisphere::~rod_multisphere(){

}
