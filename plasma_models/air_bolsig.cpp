/*!
  @file   air_bolsig.cpp
  @brief  Implementation of air_bolsig.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "air_bolsig.H"

std::string air_bolsig::s_bolsig_path  = "./plasma_models/bolsig";
std::string air_bolsig::s_lxcat_path   = "./plasma_models/bolsig";
std::string air_bolsig::s_lxcat_file   = "LXCat-June2013.txt";
std::string air_bolsig::s_script_file  = "bolsig_in.";
std::string air_bolsig::s_data_file    = "bolsig_out.";
std::string air_bolsig::s_bolsig_alpha = "E/N (Td)	Townsend ioniz. coef. alpha/N (m2)";
std::string air_bolsig::s_bolsig_mob   = "E/N (Td)	Mobility *N (1/m/V/s)";
std::string air_bolsig::s_bolsig_diff  = "E/N (Td)	Diffusion coefficient *N (1/m/s)";
std::string air_bolsig::s_bolsig_eta   = "E/N (Td)	Townsend attach. coef. eta/N (m2)";
std::string air_bolsig::s_bolsig_kex23 = "C23   N2    Excitation    12.25 eV";                              
std::string air_bolsig::s_bolsig_kex24 = "C24   N2    Excitation    13.00 eV";

air_bolsig::air_bolsig(){

}

air_bolsig::~air_bolsig(){

}

Vector<RealVect> air_bolsig::compute_velocities(const RealVect& a_E) const{

}


Vector<Real> air_bolsig::compute_source_terms(const Vector<Real>& a_cdr_densities, 
					      const Vector<Real>& a_rte_densities,
					      const RealVect&     a_E) const{

}


Vector<Real> air_bolsig::compute_rte_source_terms(const Vector<Real>& a_cdr_densities, const RealVect& a_E) const{

}

Vector<Real> air_bolsig::compute_diffusion_coefficients(const RealVect& a_E) const{

}


Vector<Real> air_bolsig::compute_conductor_fluxes(const Vector<Real>& a_extrapolated_fluxes,
						  const Vector<Real>& a_cdr_densities,
						  const Vector<Real>& a_cdr_velocities,
						  const Vector<Real>& a_rte_fluxes,
						  const RealVect&     a_E,
						  const RealVect&     a_pos,
						  const RealVect&     a_normal,
						  const Real&         a_time) const{

}


Vector<Real> air_bolsig::compute_dielectric_fluxes(const Vector<Real>& a_extrapolated_fluxes,
						   const Vector<Real>& a_cdr_densities,
						   const Vector<Real>& a_cdr_velocities,
						   const Vector<Real>& a_rte_fluxes,
						   const RealVect&     a_E,
						   const RealVect&     a_pos,
						   const RealVect&     a_normal,
						   const Real&         a_time) const{

}

Real air_bolsig::initial_sigma(const RealVect& a_pos) const {

}

void air_bolsig::call_bolsig(){

}

void air_bolsig::build_script(){

}

void air_bolsig::build_data(){

}

void air_bolsig::extract_data(){

}

void air_bolsig::delete_data(){

}

void air_bolsig::delete_script(){

}
