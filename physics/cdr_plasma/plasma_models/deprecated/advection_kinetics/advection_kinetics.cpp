/*!
  @file   advection_kinetics.cpp
  @brief  Implementation of advection_kinetics.H
  @author Robert Marskar
  @date   July 2018
*/

#include "advection_kinetics.H"

#include <ParmParse.H>

advection_kinetics::advection_kinetics(){
  m_num_species = 1;
  m_num_photons = 0;

  m_species.resize(m_num_species);
  m_species[0] = RefCountedPtr<species> (new advection_kinetics::phi_advect());

  m_outflow = false;

  { // Get parameters from input script
    ParmParse pp("advection_kinetics");
    std::string str;

    if(pp.contains("which_ebbc")){
      pp.get("which_ebbc", str);
      if(str == "wall"){
	m_outflow = false;
      }
      else if(str == "outflow"){
	m_outflow = true;
      }
      else {
	MayDay::Abort("advection_kinetics::advection_kinetics - unknown ebbc type requested");
      }
    }
  }
}

advection_kinetics::~advection_kinetics(){

}

Vector<Real> advection_kinetics::compute_cdr_diffusion_coefficients(const Real         a_time,
								    const RealVect     a_pos,
								    const RealVect     a_E,
								    const Vector<Real> a_cdr_densities) const {
  Vector<Real> diffco(m_num_species, 0.0);
  return diffco;
}

Vector<RealVect> advection_kinetics::compute_cdr_velocities(const Real         a_time,
							    const RealVect     a_pos,
							    const RealVect     a_E,
							    const Vector<Real> a_cdr_densities) const {
  Vector<RealVect> velo(m_num_species);
  velo[0] = a_E;
  return velo;
}

Vector<Real> advection_kinetics::compute_cdr_source_terms(const Real             a_time,
							  const Real             a_kappa,
							  const Real             a_dx,
							  const RealVect         a_pos,
							  const RealVect         a_E,
							  const RealVect         a_gradE,
							  const Vector<Real>     a_cdr_densities,
							  const Vector<Real>     a_rte_densities,
							  const Vector<RealVect> a_grad_cdr) const {
  Vector<Real> source(m_num_species, 0.0);
  return source;
}

Vector<Real> advection_kinetics::compute_cdr_electrode_fluxes(const Real         a_time,
							      const RealVect     a_pos,
							      const RealVect     a_normal,
							      const RealVect     a_E,
							      const Vector<Real> a_cdr_densities,
							      const Vector<Real> a_cdr_velocities,
							      const Vector<Real> a_cdr_gradients,
							      const Vector<Real> a_rte_fluxes,
							      const Vector<Real> a_extrap_cdr_fluxes) const {
  Vector<Real> flux(m_num_species, 0.0);
  flux[0] = m_outflow ? a_extrap_cdr_fluxes[0] : 0.0;
  return flux;
}

Vector<Real> advection_kinetics::compute_cdr_dielectric_fluxes(const Real         a_time,
							       const RealVect     a_pos,
							       const RealVect     a_normal,
							       const RealVect     a_E,
							       const Vector<Real> a_cdr_densities,
							       const Vector<Real> a_cdr_velocities,
							       const Vector<Real> a_cdr_gradients,
							       const Vector<Real> a_rte_fluxes,
							       const Vector<Real> a_extrap_cdr_fluxes) const {
  Vector<Real> flux(m_num_species, 0.0);
  flux[0] = m_outflow ? a_extrap_cdr_fluxes[0] : 0.0;
  return flux;
}

Vector<Real> advection_kinetics::compute_rte_source_terms(const Real         a_time,
							  const Real         a_kappa,
							  const Real         a_dx,
							  const RealVect     a_pos,
							  const RealVect     a_E,
							  const Vector<Real> a_cdr_densities) const {
  Vector<Real> source(m_num_species, 0.0);
  return source;
}

Real advection_kinetics::initial_sigma(const Real a_time, const RealVect a_pos) const {
  return 0.0;
}

advection_kinetics::phi_advect::phi_advect() {

  m_name      = "phi";
  m_unit      = "m-3";
  m_diffusive = false;
  m_mobile    = true;
  m_charge    = 1;
  m_pulse     = "square";
  
  // Get parameters from input script
  {
    ParmParse pp("advection_kinetics");
    std::string str = "true";

    Vector<Real> arr1;
    Vector<Real> arr2;
    pp.getarr("square_center", arr1, 0, SpaceDim);
    pp.getarr("square_width",  arr2, 0, SpaceDim);

    m_center = RealVect(D_DECL(arr1[0], arr1[1], arr1[2]));
    m_width  = RealVect(D_DECL(arr2[0], arr2[1], arr2[2]));

    // Check if we should turn off mobility
    if(pp.contains("mobile")){
      pp.get("mobile", str);
      if(str == "false"){
	m_mobile = false;
      }
    }
    if(pp.contains("pulse_shape")){
      pp.get("pulse_shape", str);
      if(str == "square"){
	m_pulse = "square";
      }
      if(str == "gaussian"){
	m_pulse = "gaussian";
      }
    }
  }
}

advection_kinetics::phi_advect::~phi_advect() {

}

Real advection_kinetics::phi_advect::initial_data(const RealVect a_pos, const Real a_time) const {

  const RealVect new_pos = a_pos - m_center;

  Real ret;
  if(m_pulse == "square"){
    ret = 1.0;
    for (int dir = 0; dir < SpaceDim; dir++){
      if(new_pos[dir] > 0.5*m_width[dir] || new_pos[dir] < -0.5*m_width[dir]){
	ret = 0.0;
      }
    }
  }
  else if(m_pulse == "gaussian"){
    const RealVect v1 = new_pos;
    const RealVect v2 = m_width;
    const RealVect factor = RealVect(D_DECL(v1[0]/v2[0], v1[1]/v2[1], v1[2]/v2[2]));
    const Real flen = factor.vectorLength();
    ret = exp(-0.5*flen*flen);
  }

  return ret;
#include "CD_NamespaceFooter.H"
