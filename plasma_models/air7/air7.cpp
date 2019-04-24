#include <random>

/*!
  @file   air7.cpp
  @brief  Implementation of air7.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air7.H"
#include "air7_species.H"
#include "units.H"
#include "data_ops.H"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include <PolyGeom.H>
#include <ParmParse.H>

std::string air7::s_bolsig_N2_alpha       = "C25   N2    Ionization    15.60 eV";
std::string air7::s_bolsig_O2_alpha       = "C42   O2    Ionization    12.06 eV";
std::string air7::s_bolsig_mobility       = "E/N (Td)	Mobility *N (1/m/V/s)";
std::string air7::s_bolsig_diffusion      = "E/N (Td)	Diffusion coefficient *N (1/m/s)";
std::string air7::s_bolsig_mean_energy    = "E/N (Td)	Mean energy (eV)";
std::string air7::s_skip_fit              = "Fit coefficients y=exp(A+B*ln(x)+C/x+D/x^2+E/x^3)";

air7::air7(){

  // Comment: I've taken out the dissociative excitation losses and all excitation losses

  MayDay::Warning("air7::air7 - this class is really not done...");

  m_num_species = 7;  // seven species
  m_num_photons = 3;  // Bourdon model for photons

  m_seed = 0;
  m_poiss_exp_swap = 500.;

  air7::get_gas_parameters(m_Tg, m_p, m_N, m_O2frac, m_N2frac); // Get gas parameters

  { // Emission coefficients at boundaries. Can be overridden from input script.
    m_townsend2_electrode           = 1.E-3;
    m_townsend2_dielectric          = 1.E-6;
    m_electrode_quantum_efficiency  = 1.E-2;
    m_dielectric_quantum_efficiency = 1.E-4;
    m_photoionization_efficiency    = 0.1;
    m_excitation_efficiency         = 0.6;

    ParmParse pp("air7");
    pp.get("transport_file",                  m_transport_file);
    pp.get("electrode_townsend2"       ,    m_townsend2_electrode);
    pp.get("dielectric_townsend2"       ,   m_townsend2_dielectric);
    pp.get("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
    pp.get("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
    pp.get("photoionization_efficiency",    m_photoionization_efficiency);
    pp.get("excitation_efficiency",         m_excitation_efficiency);

    pp.get("rng_seed", m_seed);
    pp.get("poiss_exp_swap", m_poiss_exp_swap);
  }

  if(m_seed < 0){ // Use time for seed
    m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }

  { // Quenching pressure
    m_pq        = 0.03947; 
    ParmParse pp("air7");
    pp.get("quenching_pressure", m_pq);
    
    m_pq *= units::s_atm2pascal;
  }

  { // Mobile ions
    m_ion_mobility   = 2.E-4;
    m_mobile_ions    = true;
    m_diffusive_ions = true;
    std::string str;

    ParmParse pp("air7");
    pp.get("ion_mobility", m_ion_mobility);
    if(pp.contains("mobile_ions")){
      pp.get("mobile_ions", str);
      if(str == "true"){
	m_mobile_ions = true;
      }
      else if(str == "false"){
	m_mobile_ions = false;
      }
    }
    if(pp.contains("diffusive_ions")){
      pp.get("diffusive_ions", str);
      if(str == "true"){
	m_diffusive_ions = true;
      }
      else if(str == "false"){
	m_diffusive_ions = false;
      }
    }
    if(pp.contains("use_fhd")){
      pp.get("use_fhd", str);
      if(str == "true"){
	m_fhd = true;
      }
      else if(str == "false"){
	m_fhd = false;
      }
    }
    if(pp.contains("alpha_corr")){
      pp.get("alpha_corr", str);
      if(str == "true"){
	m_alpha_corr = true;
      }
      else if(str == "false"){
	m_alpha_corr = false;
      }
    }
  }

  // Instantiate cdr species
  m_species.resize(m_num_species);
  m_electron_idx = 0;
  m_N2plus_idx   = 1;
  m_N4plus_idx   = 2;
  m_O2plus_idx   = 3;
  m_O4plus_idx   = 4;
  m_O2plusN2_idx = 5;
  m_O2minus_idx  = 6;
  
  m_species[m_electron_idx] = RefCountedPtr<species> (new air7::electron());
  m_species[m_N2plus_idx]   = RefCountedPtr<species> (new air7::N2plus());
  m_species[m_N4plus_idx]   = RefCountedPtr<species> (new air7::N4plus());
  m_species[m_O2plus_idx]   = RefCountedPtr<species> (new air7::O2plus());
  m_species[m_O4plus_idx]   = RefCountedPtr<species> (new air7::O4plus());
  m_species[m_O2plusN2_idx] = RefCountedPtr<species> (new air7::O2plusN2());
  m_species[m_O2minus_idx]  = RefCountedPtr<species> (new air7::O2minus());

  // Instantiate photon solvers
  m_photons.resize(m_num_photons);
  m_photon1_idx = 0;
  m_photon2_idx = 1;
  m_photon3_idx = 2;
  m_photons[m_photon1_idx] = RefCountedPtr<photon_group> (new air7::photon_one());
  m_photons[m_photon2_idx] = RefCountedPtr<photon_group> (new air7::photon_two());
  m_photons[m_photon3_idx] = RefCountedPtr<photon_group> (new air7::photon_three());

  // Compute transport coefficients
  this->compute_transport_coefficients();
  m_e_mobility.scale_y(1./m_N); // Need to scale
  m_e_diffco.scale_y(1./m_N); // Need to scale

  // Init rng
  m_rng = new std::mt19937_64(m_seed);
}

air7::~air7(){

}

void air7::get_gas_parameters(Real& a_Tg, Real& a_p, Real& a_N, Real& a_O2frac, Real& a_N2frac){
  ParmParse pp("air7");
  pp.get("gas_temperature", a_Tg);
  pp.get("gas_pressure", a_p);

  a_Tg = 300.;
  a_O2frac = 0.20;
  a_N2frac = 0.80; 

  const Real tot_frac = a_O2frac + a_N2frac; 
  a_p      = a_p*units::s_atm2pascal;
  a_O2frac = a_O2frac/tot_frac; // Normalize to one
  a_N2frac = a_N2frac/tot_frac;
  a_N      = a_p*units::s_Na/(a_Tg*units::s_R);
}

void air7::compute_transport_coefficients(){
  Real energy;
  Real entry;
  bool readLine = false;
  lookup_table* which_table = NULL;
  std::ifstream infile(m_transport_file);
  std::string line;
  while (std::getline(infile, line)){

    // Right trim string
    line.erase(line.find_last_not_of(" \n\r\t")+1);

    if(line == air7::s_bolsig_O2_alpha){
      which_table = &m_e_O2_alpha;
      readLine = true;
      continue;
    }
    else if(line == air7::s_bolsig_N2_alpha){
      which_table = &m_e_N2_alpha;
      readLine = true;
      continue;
    }
    else if(line == air7::s_bolsig_mobility){
      which_table = &m_e_mobility;
      readLine   = true;
      continue;
    }
    else if(line == air7::s_bolsig_diffusion){
      which_table = &m_e_diffco;
      readLine   = true;
      continue;
    }
    else if(line == air7::s_bolsig_mean_energy){
      which_table = &m_e_energy;
      readLine   = true;
      continue;
    }

    // Stop when we encounter an empty line
    if((line == "" || line == s_skip_fit) && readLine){
      readLine = false;
      continue;
    }

    if(readLine){
      std::istringstream iss(line);
      if (!(iss >> energy >> entry)) {
    	continue;
      }
      which_table->add_entry(energy, entry);
    }
  }
  infile.close();
}

Vector<Real> air7::compute_cdr_diffusion_coefficients(const Real&         a_time,
						      const RealVect&     a_pos,
						      const RealVect&     a_E,
						      const Vector<Real>& a_cdr_densities) const {

  Vector<Real> diffco(m_num_species, 0.0);

  const Real EbyN        = (a_E/(m_N*units::s_Td)).vectorLength();
  diffco[m_electron_idx] = this->compute_electron_diffco(EbyN);
  if(m_diffusive_ions){
    diffco[m_N2plus_idx]   = this->compute_N2plus_diffco(m_Tg);
    diffco[m_N4plus_idx]   = this->compute_N4plus_diffco(m_Tg);
    diffco[m_O2plus_idx]   = this->compute_O2plus_diffco(m_Tg);
    diffco[m_O4plus_idx]   = this->compute_O4plus_diffco(m_Tg);
    diffco[m_O2plusN2_idx] = this->compute_O2plusN2_diffco(m_Tg);
    diffco[m_O2minus_idx]  = this->compute_O2minus_diffco(m_Tg);
  }

  return diffco;
}

Vector<RealVect> air7::compute_cdr_velocities(const Real&         a_time,
					      const RealVect&     a_pos,
					      const RealVect&     a_E,
					      const Vector<Real>& a_cdr_densities) const {
  Vector<RealVect> velocities(m_num_species, RealVect::Zero);

  const Real EbyN             = (a_E/(m_N*units::s_Td)).vectorLength();
  velocities[m_electron_idx]   = -a_E*this->compute_electron_mobility(EbyN);
  if(m_mobile_ions){
    velocities[m_N2plus_idx]   =  a_E*this->compute_N2plus_mobility(EbyN);
    velocities[m_N4plus_idx]   =  a_E*this->compute_N4plus_mobility(EbyN);
    velocities[m_O2plus_idx]   =  a_E*this->compute_O2plus_mobility(EbyN);
    velocities[m_O4plus_idx]   =  a_E*this->compute_O4plus_mobility(EbyN);
    velocities[m_O2plusN2_idx] =  a_E*this->compute_O2plusN2_mobility(EbyN);
    velocities[m_O2minus_idx]  = -a_E*this->compute_O2minus_mobility(EbyN);
  }

  return velocities;
}
  
Vector<Real> air7::compute_cdr_source_terms(const Real              a_time,
					    const Real              a_kappa,
					    const Real              a_dx,
					    const RealVect&         a_pos,
					    const RealVect&         a_E,
					    const RealVect&         a_gradE,
					    const Vector<Real>&     a_cdr_densities,
					    const Vector<Real>&     a_rte_densities,
					    const Vector<RealVect>& a_grad_cdr) const {
  Vector<Real> source(m_num_species, 0.0);

#if 0 // Debug
  return source;
#endif

  // Reduced field and electron temperature
  const Real EbyN = (a_E/(m_N*units::s_Td)).vectorLength();
  const Real Te   = this->compute_Te(EbyN);
  
  // Rate constants
  const Real k1  = this->compute_electron_N2_impact_ionization(EbyN);
  const Real k2  = this->compute_electron_O2_impact_ionization(EbyN);
  const Real k3  = this->compute_N2plus_N2_M_to_N4plus_M();                         
  const Real k4  = this->compute_N4plus_O2_to_O2_2N2();
  const Real k5  = this->compute_N2plus_O2_to_O2plus_N2(m_Tg);
  const Real k6  = this->compute_O2plus_2N2_to_O2plusN2_N2(m_Tg);
  const Real k7  = this->compute_O2plusN2_N2_to_O2plus_2N2(m_Tg);
  const Real k8  = this->compute_O2plusN2_O2_to_O4plus_N2();
  const Real k9  = this->compute_O2plus_O2_M_to_O4plus_M(m_Tg);
  const Real k10 = this->compute_e_O4plus_to_2O2(Te);
  const Real k11 = this->compute_e_O2plus_to_O2(Te);
  const Real k12 = this->compute_e_2O2_to_O2minus_O2(Te);
  const Real k13 = this->compute_O2minus_O4plus_to_3O2();
  const Real k14 = this->compute_O2minus_O4plus_M_to_3O2_M(m_Tg);
  const Real k15 = this->compute_O2minus_O2plus_M_to_2O2_M(m_Tg);
  const Real k16 = this->compute_Oplus_O2_to_O_O2(m_Tg);

  const Real n_N2    = m_N*m_N2frac;
  const Real n_O2    = m_N*m_O2frac;
  
  const Real n_e     = a_cdr_densities[m_electron_idx];
  const Real n_N2p   = a_cdr_densities[m_N2plus_idx];
  const Real n_N4p   = a_cdr_densities[m_N4plus_idx];
  const Real n_O2p   = a_cdr_densities[m_O2plus_idx];
  const Real n_O4p   = a_cdr_densities[m_O4plus_idx];
  const Real n_O2pN2 = a_cdr_densities[m_O2plusN2_idx];
  const Real n_O2m   = a_cdr_densities[m_O2minus_idx];

  Real products, p, mean;

  const Real vol     = pow(a_dx, SpaceDim);
  const RealVect ve  = -a_E*this->compute_electron_mobility(EbyN);
  const Real De      =      this->compute_electron_diffco(EbyN);

  Real k1_corr;
  Real k2_corr;
  if(m_alpha_corr){
    const Real factor  = PolyGeom::dot(a_E,De*a_grad_cdr[m_electron_idx])/((1.0 + n_e)*PolyGeom::dot(ve, a_E));
    k1_corr = k1*Max(0.0, (1.0 - Max(factor, 0.0)));
    k2_corr = k2*Max(0.0, (1.0 - Max(factor, 0.0)));
  }
  else{
    k1_corr = k1;
    k2_corr = k2;
  }

  
  // k1 reaction
  p = k1_corr * n_e * n_N2;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_electron_idx]  += products;
  source[m_N2plus_idx]    += products;

  // k2 reaction
  p = k2_corr * n_e * n_O2;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_electron_idx] += products;
  source[m_O2plus_idx]   += products;

  // k3 reaction. 
  p = k3 * n_N2p * n_N2 * m_N;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_N2plus_idx] -= products;
  source[m_N4plus_idx] += products;

  // k4 reaction
  p = k4 * n_N4p * n_O2;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_N4plus_idx]  -= products;
  source[m_O2plus_idx]  += products;

  // k5 reaction
  p = k5 * n_N2p * n_O2;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_N2plus_idx] -= products;
  source[m_O2plus_idx] += products;

  // k6 reaction
  p = k6 * n_O2p * n_N2 * n_N2;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_O2plus_idx]   -= products;
  source[m_O2plusN2_idx] += products;

  // k7 reaction
  p = k7 * n_O2pN2 * n_N2;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_O2plusN2_idx] -= products;
  source[m_O2plus_idx]   += products;

  // k8 reaction
  p = k8 * n_O2pN2 * n_O2;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_O2plusN2_idx] -= products;
  source[m_O4plus_idx]   += products;

  // k9 reaction
  p = k9 * n_O2p * n_O2 * m_N;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_O2plus_idx] -= products;
  source[m_O4plus_idx] += products;

  // k10 reaction
  p = k10 * n_e * n_O4p;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_electron_idx] -= products;
  source[m_O4plus_idx]   -= products;

  // k11 reaction
  p = k11 * n_e * n_O2p;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_electron_idx] -= products;
  source[m_O2plus_idx]   -= products;

  // k12 reaction
  p = k12 * n_e * n_O2 * n_O2;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_electron_idx] -= products;
  source[m_O2minus_idx]  += products;

  // k13 reaction
  p = k13 * n_O2m * n_O4p;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_O2minus_idx] -= products;
  source[m_O4plus_idx]  -= products;

  // k14 reaction
  p = k14 * n_O2m * n_O4p * m_N;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_O2minus_idx] -= products;
  source[m_O4plus_idx]  -= products;

  // k15 reaction
  p = k15 * n_O2m * n_O2p * m_N;
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  source[m_O2minus_idx] -= products;
  source[m_O2plus_idx]  -= products;



  // Photoionization gamma + O2 -> e + O2+
  const air7::photon_one*   photon1 = static_cast<air7::photon_one*>   (&(*m_photons[m_photon1_idx]));
  const air7::photon_two*   photon2 = static_cast<air7::photon_two*>   (&(*m_photons[m_photon2_idx]));
  const air7::photon_three* photon3 = static_cast<air7::photon_three*> (&(*m_photons[m_photon3_idx]));
  p = m_photoionization_efficiency*units::s_c0*m_O2frac*m_p*(photon1->get_A()*a_rte_densities[m_photon1_idx]
								    + photon2->get_A()*a_rte_densities[m_photon2_idx]
								    + photon3->get_A()*a_rte_densities[m_photon3_idx]);
  products = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;

  source[m_electron_idx] += products;
  source[m_O2plus_idx]   += products;

  return source;
}

Vector<Real> air7::compute_cdr_fluxes(const Real&         a_time,
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
  }

  // Drift outflow for now
  for (int i = 0; i < m_num_species; i++){
    fluxes[i] = aj[i]*a_extrap_cdr_fluxes[i];
  }

  // Add secondary emission
  if(cathode){
    fluxes[m_electron_idx] += -a_rte_fluxes[m_photon1_idx]*a_quantum_efficiency;
    fluxes[m_electron_idx] += -a_rte_fluxes[m_photon2_idx]*a_quantum_efficiency;
    fluxes[m_electron_idx] += -a_rte_fluxes[m_photon3_idx]*a_quantum_efficiency;

    fluxes[m_electron_idx] += -fluxes[m_O2plus_idx]*a_townsend2;
    fluxes[m_electron_idx] += -fluxes[m_N2plus_idx]*a_townsend2;
    fluxes[m_electron_idx] += -fluxes[m_N4plus_idx]*a_townsend2;
    fluxes[m_electron_idx] += -fluxes[m_O4plus_idx]*a_townsend2;
    fluxes[m_electron_idx] += -fluxes[m_O2plusN2_idx]*a_townsend2;
  }

  return fluxes;
}

Vector<Real> air7::compute_cdr_electrode_fluxes(const Real&         a_time,
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

Vector<Real> air7::compute_cdr_dielectric_fluxes(const Real&         a_time,
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

Vector<Real> air7::compute_rte_source_terms(const Real&         a_time,
					    const RealVect&     a_pos,
					    const RealVect&     a_E,
					    const Vector<Real>& a_cdr_densities) const {

  // We take the source terms as Se = alpha*Ne*ve
  Vector<Real> ret(m_num_photons, 0.0);

  const Real EbyN  = (a_E/(m_N*units::s_Td)).vectorLength();
  const Real k1    = this->compute_electron_N2_impact_ionization(EbyN);
  const Real p    = k1*a_cdr_densities[m_electron_idx]*m_N*m_N2frac;

  const Real vol = 1.E-12;
  const Real Se  = m_fhd ? stochastic_reaction(p, vol, m_dt) : p;
  ret[m_photon1_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));
  ret[m_photon2_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));
  ret[m_photon3_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));

  return ret;
}

Real air7::initial_sigma(const Real a_time, const RealVect& a_pos) const {
  return 0.0;
}

Real air7::compute_Te(const Real a_EbyN) const{
  const Real electron_energy = m_e_energy.get_entry(a_EbyN);
  return 2.0*electron_energy*units::s_Qe/(3.0*units::s_kb);
}

Real air7::compute_electron_mobility(const Real a_EbyN) const {return m_e_mobility.get_entry(a_EbyN);}
Real air7::compute_N2plus_mobility(const Real a_EbyN)   const {return m_ion_mobility;}
Real air7::compute_N4plus_mobility(const Real a_EbyN)   const {return m_ion_mobility;}
Real air7::compute_O2plus_mobility(const Real a_EbyN)   const {return m_ion_mobility;}
Real air7::compute_O4plus_mobility(const Real a_EbyN)   const {return m_ion_mobility;}
Real air7::compute_O2plusN2_mobility(const Real a_EbyN) const {return m_ion_mobility;}
Real air7::compute_O2minus_mobility(const Real a_EbyN)  const {return m_ion_mobility;}

Real air7::compute_electron_diffco(const Real a_EbyN) const {//
  return m_e_diffco.get_entry(a_EbyN);
}
Real air7::compute_N2plus_diffco(const Real a_Tg)     const {return units::s_kb*a_Tg*m_ion_mobility/(units::s_Qe);}
Real air7::compute_N4plus_diffco(const Real a_Tg)     const {return units::s_kb*a_Tg*m_ion_mobility/(units::s_Qe);}
Real air7::compute_O2plus_diffco(const Real a_Tg)     const {return units::s_kb*a_Tg*m_ion_mobility/(units::s_Qe);}
Real air7::compute_O4plus_diffco(const Real a_Tg)     const {return units::s_kb*a_Tg*m_ion_mobility/(units::s_Qe);}
Real air7::compute_O2plusN2_diffco(const Real a_Tg)   const {return units::s_kb*a_Tg*m_ion_mobility/(units::s_Qe);}
Real air7::compute_O2minus_diffco(const Real a_Tg)    const {return units::s_kb*a_Tg*m_ion_mobility/(units::s_Qe);}

Real air7::compute_electron_N2_impact_ionization(const Real a_EbyN) const {return m_e_N2_alpha.get_entry(a_EbyN);}
Real air7::compute_electron_O2_impact_ionization(const Real a_EbyN) const {return m_e_O2_alpha.get_entry(a_EbyN);}
Real air7::compute_N2plus_N2_M_to_N4plus_M()                        const {return 5.E-41;}
Real air7::compute_N4plus_O2_to_O2_2N2()                            const {return 2.5E-16;}
Real air7::compute_N2plus_O2_to_O2plus_N2(const Real a_Tg)          const {return 1.05E-15/sqrt(a_Tg);}
Real air7::compute_O2plus_2N2_to_O2plusN2_N2(const Real a_Tg)       const {return 8.1E-38/(a_Tg*a_Tg);}
Real air7::compute_O2plusN2_N2_to_O2plus_2N2(const Real a_Tg)       const {return 14.8*pow(a_Tg, -5.3)*exp(-2357.0/a_Tg);}
Real air7::compute_O2plusN2_O2_to_O4plus_N2()                       const {return 1.E-15;}
Real air7::compute_O2plus_O2_M_to_O4plus_M(const Real a_Tg)         const {return 2.03E-34*pow(a_Tg, -3.2);}
Real air7::compute_e_O4plus_to_2O2(const Real a_Te)                 const {return 2.42E-11/(sqrt(a_Te));}
Real air7::compute_e_O2plus_to_O2(const Real a_Te)                  const {return 6.E-11/a_Te;}
Real air7::compute_e_2O2_to_O2minus_O2(const Real a_Te)             const {return 6E-39/a_Te;}
Real air7::compute_O2minus_O4plus_to_3O2()                          const {return 1.E-13;}
Real air7::compute_O2minus_O4plus_M_to_3O2_M(const Real a_Tg)       const {return 3.12E-31*pow(a_Tg, -2.5);}
Real air7::compute_O2minus_O2plus_M_to_2O2_M(const Real a_Tg)       const {return 3.12E-31*pow(a_Tg, -2.5);}
Real air7::compute_Oplus_O2_to_O_O2(const Real a_Tg)                const {return 3.46E-12/sqrt(a_Tg);}

Real air7::stochastic_reaction(const Real a_S, const Real a_vol, const Real a_dt) const{
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
