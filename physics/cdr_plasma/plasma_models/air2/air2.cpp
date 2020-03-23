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
  pp.get("quenching_pressure",            m_pq);

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
}

air2::~air2(){

}

Vector<Real> air2::compute_cdr_diffusion_coefficients(const Real         a_time,
						      const RealVect     a_pos,
						      const RealVect     a_E,
						      const Vector<Real> a_cdr_densities) const {

  Vector<Real> diffco(m_num_species, 0.0);
  diffco[0] = 4.3628E-3*pow(a_E.vectorLength(), 0.22);

  return diffco;
}

Vector<RealVect> air2::compute_cdr_velocities(const Real         a_time,
					      const RealVect     a_pos,
					      const RealVect     a_E,
					      const Vector<Real> a_cdr_densities) const {
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

void air2::advance_reaction_network(Vector<Real>&          a_particle_sources,
				    Vector<Real>&          a_photon_sources,
				    const Vector<Real>     a_particle_densities,
				    const Vector<RealVect> a_particle_gradients,
				    const Vector<Real>     a_photon_densities,
				    const RealVect         a_E,
				    const RealVect         a_pos,
				    const Real             a_dx,
				    const Real             a_dt,
				    const Real             a_time,
				    const Real             a_kappa) const{

  const Real E      = a_E.vectorLength();
  const Real alpha  = air2::get_alpha(E);
  const Real eta    = 340.75;
  const Real mu     = air2::get_mobility(E);
  const Real vol    = pow(a_dx, SpaceDim);
  const Real pO2    = 0.2*m_p;;

  const air2::photon_one*   photon1 = static_cast<air2::photon_one*>   (&(*m_photons[m_photon1_idx]));
  const air2::photon_two*   photon2 = static_cast<air2::photon_two*>   (&(*m_photons[m_photon2_idx]));
  const air2::photon_three* photon3 = static_cast<air2::photon_three*> (&(*m_photons[m_photon3_idx]));
  
  // Impact ionization and recomb
  const Real p1 = mu*alpha*E*a_particle_densities[m_electron_idx];
  const Real p2 = mu*eta*E*a_particle_densities[m_electron_idx];


  a_particle_sources[m_electron_idx] = p1;
  a_particle_sources[m_positive_idx] = p1;

  a_particle_sources[m_electron_idx] -= p2;
  a_particle_sources[m_positive_idx] -= p2;

  // Photoionization
  const Real Sph = m_photoionization_efficiency*units::s_c0*pO2*(photon1->get_A()*a_photon_densities[m_photon1_idx]
								 + photon2->get_A()*a_photon_densities[m_photon2_idx]
								 + photon3->get_A()*a_photon_densities[m_photon3_idx]);
  a_particle_sources[m_electron_idx] += Sph;
  a_particle_sources[m_positive_idx] += Sph;


  // Emission of photons
  const Real p     = mu*alpha*E*a_particle_densities[m_electron_idx]*m_excitation_efficiency*(m_pq/(m_pq + m_p));;
  a_photon_sources[m_photon1_idx] = p;
  a_photon_sources[m_photon2_idx] = p;
  a_photon_sources[m_photon3_idx] = p;
}

Vector<Real> air2::compute_cdr_fluxes(const Real         a_time,
				      const RealVect     a_pos,
				      const RealVect     a_normal,
				      const RealVect     a_E,
				      const Vector<Real> a_cdr_densities,
				      const Vector<Real> a_cdr_velocities,
				      const Vector<Real> a_cdr_gradients,
				      const Vector<Real> a_rte_fluxes,
				      const Vector<Real> a_extrap_cdr_fluxes,
				      const Real         a_townsend2,
				      const Real         a_quantum_efficiency) const {

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

Vector<Real> air2::compute_cdr_electrode_fluxes(const Real         a_time,
						const RealVect     a_pos,
						const RealVect     a_normal,
						const RealVect     a_E,
						const Vector<Real> a_cdr_densities,
						const Vector<Real> a_cdr_velocities,
						const Vector<Real> a_cdr_gradients,
						const Vector<Real> a_rte_fluxes,
						const Vector<Real> a_extrap_cdr_fluxes) const {

  return this->compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
				  a_extrap_cdr_fluxes, m_townsend2_electrode, m_electrode_quantum_efficiency);
}

Vector<Real> air2::compute_cdr_dielectric_fluxes(const Real         a_time,
						 const RealVect     a_pos,
						 const RealVect     a_normal,
						 const RealVect     a_E,
						 const Vector<Real> a_cdr_densities,
						 const Vector<Real> a_cdr_velocities,
						 const Vector<Real> a_cdr_gradients,
						 const Vector<Real> a_rte_fluxes,
						 const Vector<Real> a_extrap_cdr_fluxes) const {

  return this->compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients, a_rte_fluxes,
				  a_extrap_cdr_fluxes, m_townsend2_dielectric, m_dielectric_quantum_efficiency);
}

Real air2::initial_sigma(const Real a_time, const RealVect a_pos) const {
  return 0.0;
}
