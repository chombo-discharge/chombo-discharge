/*!
  @file   air2.cpp
  @brief  Implementation of air2.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air2.H"
#include "air2_species.H"
#include "units.H"
#include "data_ops.H"

#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include <PolyGeom.H>
#include <ParmParse.H>

air2::air2(){
  m_num_species  = 2;  
  m_num_photons  = 3;
  
  m_electron_idx = 0;
  m_positive_idx = 1;
  
  m_photon1_idx  = 0;
  m_photon2_idx  = 1;
  m_photon3_idx  = 2;

  ParmParse pp("air2");
  std::string str;
  pp.get("electrode_townsend2"       ,    m_townsend2_electrode);
  pp.get("dielectric_townsend2"       ,   m_townsend2_dielectric);
  pp.get("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
  pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
  pp.get("photoionization_efficiency",    m_photoionization_efficiency);
  pp.get("excitation_efficiency",         m_excitation_efficiency);
  pp.get("rng_seed", m_seed);
  pp.get("poiss_exp_swap", m_poiss_exp_swap);
  pp.get("quenching_pressure", m_pq);
  pp.get("use_fhd",str);

  m_fhd  = (str == "true") ? true : false;
  m_seed = (m_seed < 0) ? std::chrono::system_clock::now().time_since_epoch().count() : m_seed;
  m_pq   = m_pq*units::s_atm2pascal;
  m_p    = 1.0*units::s_atm2pascal;

  // Instantiate cdr species
  m_species.resize(m_num_species);
  m_photons.resize(m_num_photons);
  
  m_species[m_electron_idx] = RefCountedPtr<species> (new air2::electron());
  m_species[m_positive_idx] = RefCountedPtr<species> (new air2::positive_ion());
  
  m_photons[m_photon1_idx]  = RefCountedPtr<photon_group> (new air2::photon_one());
  m_photons[m_photon2_idx]  = RefCountedPtr<photon_group> (new air2::photon_two());
  m_photons[m_photon3_idx]  = RefCountedPtr<photon_group> (new air2::photon_three());

  // Init rng
  m_rng = new std::mt19937_64(m_seed);
}

air2::~air2(){

}

Vector<Real> air2::compute_cdr_diffusion_coefficients(const Real&         a_time,
						      const RealVect&     a_pos,
						      const RealVect&     a_E,
						      const Vector<Real>& a_cdr_densities) const {

  Vector<Real> diffco(m_num_species, 0.0);
  diffco[0] = 4.3628E-3*pow(a_E.vectorLength(), 0.22);

  return diffco;
}

Vector<RealVect> air2::compute_cdr_velocities(const Real&         a_time,
					      const RealVect&     a_pos,
					      const RealVect&     a_E,
					      const Vector<Real>& a_cdr_densities) const {
  Vector<RealVect> velocities(m_num_species, RealVect::Zero);
  velocities[m_electron_idx] = -a_E*get_mobility(a_E.vectorLength());
  return velocities;
}

Real air2::get_alpha(const Real a_E) const {
  return (1.1944E6 + 4.3666E26/(a_E*a_E*a_E))*exp(-2.73E7/a_E);
}

Real air2::get_mobility(const Real a_E) const {
  return 2.3987*pow(a_E, -0.26);
}
  
Vector<Real> air2::compute_cdr_source_terms(const Real              a_time,
					    const Real              a_kappa,
					    const Real              a_dx,
					    const RealVect&         a_pos,
					    const RealVect&         a_E,
					    const RealVect&         a_gradE,
					    const Vector<Real>&     a_cdr_densities,
					    const Vector<Real>&     a_rte_densities,
					    const Vector<RealVect>& a_grad_cdr) const {
  Vector<Real> source(m_num_species, 0.0);

  const Real E      = a_E.vectorLength();
  const Real alpha  = air2::get_alpha(E);
  const Real eta    = 340.75;
  const Real mu     = air2::get_mobility(E);
  const Real vol    = pow(a_dx, SpaceDim);
  const Real pO2    = 0.2*m_p;;


  Real products, p;
  
  // Stochastic or deterministic impact ionization
  p = mu*(alpha-eta)*E*a_cdr_densities[m_electron_idx];
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_electron_idx] += products;
  source[m_positive_idx] += products;

  // Stochastic or deterministic photoionization
  const air2::photon_one*   photon1 = static_cast<air2::photon_one*>   (&(*m_photons[m_photon1_idx]));
  const air2::photon_two*   photon2 = static_cast<air2::photon_two*>   (&(*m_photons[m_photon2_idx]));
  const air2::photon_three* photon3 = static_cast<air2::photon_three*> (&(*m_photons[m_photon3_idx]));
  p = m_photoionization_efficiency*units::s_c0*pO2*(photon1->get_A()*a_rte_densities[m_photon1_idx]
						    + photon2->get_A()*a_rte_densities[m_photon2_idx]
						    + photon3->get_A()*a_rte_densities[m_photon3_idx]);
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;

  source[m_electron_idx] += products;
  source[m_positive_idx] += products;
  
		      
  return source;
}

Vector<Real> air2::compute_cdr_fluxes(const Real&         a_time,
				      const RealVect&     a_pos,
				      const RealVect&     a_normal,
				      const RealVect&     a_E,
				      const Vector<Real>& a_cdr_densities,
				      const Vector<Real>& a_cdr_velocities,
				      const Vector<Real>& a_cdr_gradients,
				      const Vector<Real>& a_rte_fluxes,
				      const Vector<Real>& a_extrap_cdr_fluxes,
				      const Real&         a_townsend2,
				      const Real&         a_quantum_efficiency) const {

  Vector<Real> fluxes(m_num_species, 0.0);  
  const bool cathode = PolyGeom::dot(a_E, a_normal) < 0.0;
  const bool anode   = PolyGeom::dot(a_E, a_normal) > 0.0;

  // Switch for setting drift flux to zero for charge species
  Vector<Real> aj(m_num_species, 0.0);
  for (int i = 0; i < m_num_species; i++){
    if(data_ops::sgn(m_species[i]->get_charge())*PolyGeom::dot(a_E, a_normal) < 0){
      aj[i] = 1.0;
    }
    else {
      aj[i] = 0.0;
    }

    fluxes[i] = aj[i]*a_extrap_cdr_fluxes[i];
  }

  // Add secondary emission
  if(cathode){
    fluxes[m_electron_idx] += -a_rte_fluxes[m_photon1_idx]*a_quantum_efficiency;
  }

  return fluxes;
}

Vector<Real> air2::compute_cdr_electrode_fluxes(const Real&         a_time,
						const RealVect&     a_pos,
						const RealVect&     a_normal,
						const RealVect&     a_E,
						const Vector<Real>& a_cdr_densities,
						const Vector<Real>& a_cdr_velocities,
						const Vector<Real>& a_cdr_gradients,
						const Vector<Real>& a_rte_fluxes,
						const Vector<Real>& a_extrap_cdr_fluxes) const {

  return this->compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
				  a_extrap_cdr_fluxes, m_townsend2_electrode, m_electrode_quantum_efficiency);
}

Vector<Real> air2::compute_cdr_dielectric_fluxes(const Real&         a_time,
						 const RealVect&     a_pos,
						 const RealVect&     a_normal,
						 const RealVect&     a_E,
						 const Vector<Real>& a_cdr_densities,
						 const Vector<Real>& a_cdr_velocities,
						 const Vector<Real>& a_cdr_gradients,
						 const Vector<Real>& a_rte_fluxes,
						 const Vector<Real>& a_extrap_cdr_fluxes) const {

  return this->compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
				  a_extrap_cdr_fluxes, m_townsend2_dielectric, m_dielectric_quantum_efficiency);
}

Vector<Real> air2::compute_rte_source_terms(const Real&         a_time,
					    const RealVect&     a_pos,
					    const RealVect&     a_E,
					    const Vector<Real>& a_cdr_densities) const {

  // We take the source terms as Se = alpha*Ne*ve
  Vector<Real> ret(m_num_photons, 0.0);

  const Real E     = a_E.vectorLength();
  const Real mu    = air2::get_mobility(E);
  const Real alpha = air2::get_alpha(a_E.vectorLength());
  const Real p     = mu*alpha*E*a_cdr_densities[m_electron_idx]*m_excitation_efficiency;
  const Real vol = 1.E-12;
  const Real Se  = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  
  ret[m_photon1_idx] = Se*(m_pq/(m_pq + m_p));
  ret[m_photon2_idx] = Se*(m_pq/(m_pq + m_p));
  ret[m_photon3_idx] = Se*(m_pq/(m_pq + m_p));

  return ret;
}

Real air2::initial_sigma(const Real a_time, const RealVect& a_pos) const {
  return 0.0;
}

Real air2::stochastic_reaction(const Real a_S, const Real a_vol, const Real a_dt) const{
  Real value = 0.0;
  const Real mean = a_S*a_vol*a_dt;
  if(mean < m_poiss_exp_swap){
    std::poisson_distribution<int> dist(mean);
    value = Max(0.0, 1.0*dist(*m_rng)/(a_vol * a_dt));
  }
  else{
    std::normal_distribution<double> dist(mean, sqrt(mean));
    value = Max(0.0, 1.0*dist(*m_rng)/(a_vol * a_dt));
  }

  return value;
}
