/*!
  @file   air9eed.cpp
  @brief  Implementation of air9eed.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air9eed.H"
#include "air9eed_species.H"
#include "units.H"
#include "data_ops.H"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <PolyGeom.H>
#include <ParmParse.H>

std::string air9eed::s_bolsig_energy_E = "Energy (eV) 	Electric field / N (Td)";
std::string air9eed::s_bolsig_mobility = "Energy (eV)	Mobility *N (1/m/V/s)";
std::string air9eed::s_bolsig_N2_alpha = "C25   N2    Ionization    15.60 eV";
std::string air9eed::s_bolsig_O2_alpha = "C42   O2    Ionization    12.06 eV";


air9eed::air9eed(){

  MayDay::Warning("air9eed::air9eed - this class is really not done...");

  m_num_species = 9;  // 8 reactive ones plus the eed
  m_num_photons = 3;  // Bourdon model for photons
  m_eed_solve = true;
  m_eed_index = 0;

  air9eed::get_gas_parameters(m_Tg, m_p, m_N, m_O2frac, m_N2frac); // Get gas parameters

  { // Emission coefficients at boundaries. Can be overridden from input script.
    m_townsend2_electrode           = 1.E-3;
    m_townsend2_dielectric          = 1.E-6;
    m_electrode_quantum_efficiency  = 1.E-2;
    m_dielectric_quantum_efficiency = 1.E-4;
    m_photoionization_efficiency    = 0.1;
    m_excitation_efficiency         = 0.6;

    ParmParse pp("air9eed");
    pp.get("transport_file",                  m_transport_file);
    pp.query("electrode_townsend2"       ,    m_townsend2_electrode);
    pp.query("dielectric_townsend2"       ,   m_townsend2_dielectric);
    pp.query("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
    pp.query("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
    pp.query("photoionization_efficiency",    m_photoionization_efficiency);
    pp.query("excitation_efficiency",         m_excitation_efficiency);
  }

  { // Quenching pressure
    m_pq        = 0.03947; 
    ParmParse pp("air9eed");
    pp.query("quenching_pressure", m_pq);
    
    m_pq *= units::s_atm2pascal;
  }

  // Instantiate cdr species
  m_species.resize(m_num_species);
  m_eed_idx      = 0;
  m_electron_idx = 1;
  m_N2plus_idx   = 2;
  m_N4plus_idx   = 3;
  m_O2plus_idx   = 4;
  m_O4plus_idx   = 5;
  m_O2plusN2_idx = 6;
  m_O2minus_idx  = 7;
  m_Ominus_idx   = 8;
  
  m_species[m_eed_idx]      = RefCountedPtr<species> (new air9eed::eed());
  m_species[m_electron_idx] = RefCountedPtr<species> (new air9eed::electron());
  m_species[m_N2plus_idx]   = RefCountedPtr<species> (new air9eed::N2plus());
  m_species[m_N4plus_idx]   = RefCountedPtr<species> (new air9eed::N4plus());
  m_species[m_O2plus_idx]   = RefCountedPtr<species> (new air9eed::O2plus());
  m_species[m_O4plus_idx]   = RefCountedPtr<species> (new air9eed::O4plus());
  m_species[m_O2plusN2_idx] = RefCountedPtr<species> (new air9eed::O2plusN2());
  m_species[m_O2minus_idx]  = RefCountedPtr<species> (new air9eed::O2minus());
  m_species[m_Ominus_idx]   = RefCountedPtr<species> (new air9eed::Ominus());

  // Instantiate photon solvers
  m_photons.resize(m_num_photons);
  m_photon1_idx = 0;
  m_photon2_idx = 1;
  m_photon3_idx = 2;
  m_photons[m_photon1_idx] = RefCountedPtr<photon_group> (new air9eed::photon_one());
  m_photons[m_photon2_idx] = RefCountedPtr<photon_group> (new air9eed::photon_two());
  m_photons[m_photon3_idx] = RefCountedPtr<photon_group> (new air9eed::photon_three());

  // Compute transport coefficients
  this->compute_transport_coefficients();
  m_e_mobility.scale_y(1./m_N); // Need to scale
  m_init_eed.swap_xy();         // Input table is in reverse order


  //  m_e_mobility.dump_table();
  m_init_eed.dump_table();
}

air9eed::~air9eed(){

}

void air9eed::get_gas_parameters(Real& a_Tg, Real& a_p, Real& a_N, Real& a_O2frac, Real& a_N2frac){
  ParmParse pp("air9eed");
  pp.query("gas_temperature", a_Tg);
  pp.query("gas_pressure", a_p);

  a_Tg = 300.;
  a_O2frac = 0.21;
  a_N2frac = 0.79; 

  const Real tot_frac = a_O2frac + a_N2frac; 
  a_p      = a_p*units::s_atm2pascal;
  a_O2frac = a_O2frac/tot_frac; // Normalize to one
  a_N2frac = a_N2frac/tot_frac;
  a_N      = a_p*units::s_Na/(a_Tg*units::s_R);
}

void air9eed::compute_transport_coefficients(){
  Real energy;
  Real entry;
  bool readLine = false;
  lookup_table* which_table = NULL;
  std::ifstream infile(m_transport_file);
  std::string line;
  while (std::getline(infile, line)){

    // Right trim string
    line.erase(line.find_last_not_of(" \n\r\t")+1);

    if(line == air9eed::s_bolsig_energy_E){
      which_table = &m_init_eed;
      readLine = true;
      continue;
    }
    else if(line == air9eed::s_bolsig_mobility){
      which_table = &m_e_mobility;
      readLine   = true;
      continue;
    }
    else if(line == air9eed::s_bolsig_N2_alpha){
      which_table = &m_e_N2_alpha;
      readLine = true;
      continue;
    }
    else if(line == air9eed::s_bolsig_O2_alpha){
      which_table = &m_e_O2_alpha;
      readLine = true;
      continue;
    }

    // Stop when we encounter an empty line
    if(line == "" && readLine){
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

Vector<Real> air9eed::compute_cdr_diffusion_coefficients(const Real&         a_time,
							 const RealVect&     a_pos,
							 const RealVect&     a_E,
							 const Vector<Real>& a_cdr_densities) const {

  Vector<Real> diffco(m_num_species, 0.0);

  const Real electron_energy = a_cdr_densities[m_eed_idx]/(1.0 + a_cdr_densities[m_electron_idx]);
  const Real EbyN            = (a_E/(m_N*units::s_Td)).vectorLength();
  
  diffco[m_eed_idx]      = this->compute_eed_diffco(electron_energy, m_N);
  diffco[m_electron_idx] = this->compute_electron_diffco(electron_energy, m_N);
  diffco[m_N2plus_idx]   = this->compute_N2plus_diffco();
  diffco[m_N4plus_idx]   = this->compute_N4plus_diffco();
  diffco[m_O2plus_idx]   = this->compute_O2plus_diffco();
  diffco[m_O4plus_idx]   = this->compute_O4plus_diffco();
  diffco[m_O2plusN2_idx] = this->compute_O2plusN2_diffco();
  diffco[m_O2minus_idx]  = this->compute_O2minus_diffco();

  return diffco;
}

Vector<RealVect> air9eed::compute_cdr_velocities(const Real&         a_time,
						 const RealVect&     a_pos,
						 const RealVect&     a_E,
						 const Vector<Real>& a_cdr_densities) const {
  Vector<RealVect> velocities(m_num_species, RealVect::Zero);


  const Real electron_energy = a_cdr_densities[m_eed_idx]/(1.0 + a_cdr_densities[m_electron_idx]);
  const Real EbyN            = (a_E/(m_N*units::s_Td)).vectorLength();

  velocities[m_eed_idx]      = this->compute_eed_mobility(electron_energy, m_N)*(-a_E);
  velocities[m_electron_idx] = this->compute_electron_mobility(electron_energy, m_N)*(-a_E);
  velocities[m_N2plus_idx]   = this->compute_N2plus_mobility(EbyN)*a_E;
  velocities[m_N4plus_idx]   = this->compute_N4plus_mobility(EbyN)*a_E;
  velocities[m_O2plus_idx]   = this->compute_O2plus_mobility(EbyN)*a_E;
  velocities[m_O4plus_idx]   = this->compute_O4plus_mobility(EbyN)*a_E;
  velocities[m_O2plusN2_idx] = this->compute_O2plusN2_mobility(EbyN)*a_E;
  velocities[m_O2minus_idx]  = this->compute_O2minus_mobility(EbyN)*(-a_E);
  velocities[m_Ominus_idx]   = this->compute_Ominus_mobility(EbyN)*(-a_E);

  return velocities;
}
  
Vector<Real> air9eed::compute_cdr_source_terms(const Real              a_time,
					       const RealVect&         a_pos,
					       const RealVect&         a_E,
					       const RealVect&         a_gradE,
					       const Vector<Real>&     a_cdr_densities,
					       const Vector<Real>&     a_rte_densities,
					       const Vector<RealVect>& a_grad_cdr) const {
  Vector<Real> source(m_num_species, 0.0); 

  const Real electron_energy = a_cdr_densities[m_eed_idx]/(1.E0 + a_cdr_densities[m_electron_idx]); // eV
  const Real Te              = 2.0*(electron_energy*units::s_Qe)/(3.0*units::s_kb);  // Kelvin
  const Real EbyN            = (a_E/(m_N*units::s_Td)).vectorLength();

  // Room for improvement: The best thing would be to store the rate coefficients as matrices and then do S = K*n
  
  // Get all rate constant
  const Real k1  = this->compute_electron_N2_impact_ionization(electron_energy, m_N); 
  const Real k2  = this->compute_electron_O2_impact_ionization(electron_energy, m_N); 
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
  const Real k16 = this->compute_e_O2_to_e_2O_c1(electron_energy, m_N);
  const Real k17 = this->compute_e_O2_to_e_2O_c2(electron_energy, m_N);
  const Real k18 = this->compute_e_O2_to_Ominus_O(electron_energy, m_N);
  const Real k19 = this->compute_Oplus_O2_to_O_O2(m_Tg);
  const Real k20 = this->compute_e_N2_to_e_N2(electron_energy, m_N);
  const Real k21 = this->compute_e_O2_to_e_O2(electron_energy, m_N);

  // Electron energy losses for electron collisions
  const Real dE_k1  = this->compute_e_N2_ionization_loss();
  const Real dE_k2  = this->compute_e_O2_ionization_loss();
  const Real dE_k16 = this->compute_e_O2_dissociation_loss_c1();
  const Real dE_k17 = this->compute_e_O2_dissociation_loss_c2();
  const Real dE_k18 = this->compute_e_O2_dissociative_attachment_loss();
  const Real dE_k20 = this->compute_e_N2_scattering_loss();
  const Real dE_k21 = this->compute_e_O2_scattering_loss();

  const Real n_e     = a_cdr_densities[m_electron_idx];
  const Real n_N2    = m_N*m_N2frac;
  const Real n_O2    = m_N*m_O2frac;
  const Real n_N2p   = a_cdr_densities[m_N2plus_idx];
  const Real n_N4p   = a_cdr_densities[m_N4plus_idx];
  const Real n_O2p   = a_cdr_densities[m_O2plus_idx];
  const Real n_O4p   = a_cdr_densities[m_O4plus_idx];
  const Real n_O2pN2 = a_cdr_densities[m_O2plusN2_idx];
  const Real n_O2m   = a_cdr_densities[m_O2minus_idx];
  const Real n_Om    = a_cdr_densities[m_Ominus_idx];

  Real loss;
  Real products;

  // Compute electron velocity, both drift and diffusion
  const RealVect ve = this->compute_electron_mobility(electron_energy, m_N)*(-a_E);
  const Real     De = this->compute_electron_diffco(electron_energy, m_N);
  const RealVect je = ve*n_e - De*a_grad_cdr[m_electron_idx];

  // Joule heating
  loss = PolyGeom::dot(-je, a_E); 
  source[m_eed_idx] += loss;
  
  // k1 reaction
  loss     = dE_k1;
  products = k1 * n_e * n_N2;
  source[m_eed_idx]       -= products*loss; 
  source[m_electron_idx]  += products;
  source[m_N2plus_idx]    += products;

  // k2 reaction
  loss     = dE_k2; 
  products = k2 * n_e * n_O2;
  source[m_eed_idx]      -= products*loss;
  source[m_electron_idx] += products;
  source[m_O2plus_idx]   += products;

#if 1 // Debug
  return source;
#endif

  // k3 reaction. 
  products = k3 * n_N2p * n_N2 * (n_N2 + n_O2);
  source[m_N2plus_idx] -= products;
  source[m_N4plus_idx] += products;

  // k4 reaction
  products = k4 * n_N4p * n_O2;
  source[m_N4plus_idx]  -= products;
  source[m_O2plus_idx]  += products;
  source[m_N2plus_idx]  += 2*products;

  // k5 reaction
  products = k5 * n_N2p * n_O2;
  source[m_N2plus_idx] -= products;
  source[m_O2plus_idx] += products;

  // k6 reaction
  products = k6 * n_O2p * 2.0*n_N2;
  source[m_O2plus_idx]   -= products;
  source[m_O2plusN2_idx] += products;

  // k7 reaction
  products = k7 * n_O2pN2 * n_N2;
  source[m_O2plusN2_idx] -= products;
  source[m_O2plus_idx]   += products;

  // k8 reaction
  products = k8 * n_O2pN2 * n_O2;
  source[m_O2plusN2_idx] -= products;
  source[m_O4plus_idx]   += products;

  // k9 reaction
  products = k9 * n_O2p * n_O2 * (n_O2 + n_N2);
  source[m_O2plus_idx] -= products;
  source[m_O4plus_idx] += products;

  // k10 reaction
  products = k10 * n_e * n_O4p;
  source[m_electron_idx] -= products;
  source[m_O4plus_idx]   -= products;

  // k11 reaction
  products = k11 * n_e * n_O2p;
  source[m_electron_idx] -= products;
  source[m_O2plus_idx]   -= products;

  // k12 reaction
  products = k12 * n_e * n_O2 * n_O2;
  source[m_electron_idx] -= products;
  source[m_O2minus_idx]  += products;

  // k13 reaction
  products = k13 * n_O2m * n_O4p;
  source[m_O2minus_idx] -= products;
  source[m_O4plus_idx]  -= products;

  // k14 reaction
  products = k14 * n_O2m * n_O4p * (n_N2 + n_O2);
  source[m_O2minus_idx] -= products;
  source[m_O4plus_idx]  -= products;

  // k15 reaction
  products = k15 * n_O2m * n_O2p * (n_O2 + n_N2);
  source[m_O2minus_idx] -= products;
  source[m_O2plus_idx]  -= products;

  // k16 reaction
  loss     = dE_k16;
  products = k16 * n_e * n_O2;
  source[m_eed_idx] -= products*loss;

  // k17 reaction
  loss     = dE_k17;
  products = k17 * n_e * n_O2;
  source[m_eed_idx] -= products*loss;

  // k18 reaction
  loss     = dE_k18;
  products = k18 * n_e * n_O2;
  source[m_eed_idx]      -= products*loss;
  source[m_electron_idx] -= products;
  source[m_Ominus_idx]   += products;

  // k19 reaction
  products = k19 * n_Om * n_O2p;
  source[m_Ominus_idx] -= products;
  source[m_O2plus_idx] -= products;

  // k20 reaction
  loss     = dE_k20;
  products = k20 * n_e * n_N2;
  source[m_eed_idx] -= products*loss;

  // k21 reaction
  loss     = dE_k21;
  products = k21 * n_e * n_O2;
  source[m_eed_idx] -= products*loss;

  // Photoionization gamma + O2 -> e + O2+
  const air9eed::photon_one*   photon1 = static_cast<air9eed::photon_one*>   (&(*m_photons[m_photon1_idx]));
  const air9eed::photon_two*   photon2 = static_cast<air9eed::photon_two*>   (&(*m_photons[m_photon2_idx]));
  const air9eed::photon_three* photon3 = static_cast<air9eed::photon_three*> (&(*m_photons[m_photon3_idx]));
  products = m_photoionization_efficiency*units::s_c0*m_O2frac*m_p*(photon1->get_A()*a_rte_densities[m_photon1_idx]
								    + photon2->get_A()*a_rte_densities[m_photon2_idx]
								    + photon3->get_A()*a_rte_densities[m_photon3_idx]);

  source[m_electron_idx] += products;
  source[m_O2plus_idx]   += products;


  return source;
}

Vector<Real> air9eed::compute_cdr_fluxes(const Real&         a_time,
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
  
  const Real electron_energy = a_cdr_densities[m_eed_idx]/(1.0 + a_cdr_densities[m_electron_idx]);
  const Real Te              = 2.0*electron_energy*units::s_Qe/(3.0*units::s_kb);
  const Real EbyN            = (a_E/(m_N*units::s_Td)).vectorLength();
  const Real ion_mass        = 2.65E-26; // kg
  const Real vth_g           = sqrt(units::s_kb*m_Tg/(units::s_pi*ion_mass));  // Ion thermal velocity
  const Real vth_e           = sqrt(units::s_kb*Te/(units::s_pi*units::s_me)); // Electron thermal velocity


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

  // Drift outflow 
  for (int i = 0; i < m_num_species; i++){
    fluxes[i] = aj[i]*a_extrap_cdr_fluxes[i];
    fluxes[i] = Max(0.0, a_extrap_cdr_fluxes[i]);
  }

  return fluxes;
}

Vector<Real> air9eed::compute_cdr_electrode_fluxes(const Real&         a_time,
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

Vector<Real> air9eed::compute_cdr_dielectric_fluxes(const Real&         a_time,
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

Vector<Real> air9eed::compute_rte_source_terms(const Real&         a_time,
					       const RealVect&     a_pos,
					       const RealVect&     a_E,
					       const Vector<Real>& a_cdr_densities) const {

  // We take the source terms as Se = alpha*Ne*ve

  Vector<Real> ret(m_num_photons, 0.0);

  const Real electron_energy = a_cdr_densities[m_eed_idx]/(1.E0 + a_cdr_densities[m_electron_idx]); // eV
  const Real Te              = 2.0*(electron_energy*units::s_Qe)/(3.0*units::s_kb);  // Kelvin
  const Real EbyN            = (a_E/(m_N*units::s_Td)).vectorLength();
  const Real k1              = this->compute_electron_N2_impact_ionization(electron_energy, m_N); 
  const Real Se              = k1*a_cdr_densities[m_electron_idx]*m_N*m_N2frac;

  ret[m_photon1_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));
  ret[m_photon2_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));
  ret[m_photon3_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));

  return ret;
}

Real air9eed::initial_sigma(const Real a_time, const RealVect& a_pos) const {
  return 0.0;
}

Real air9eed::compute_eed_mobility(const Real a_energy, const Real a_N) const {
  return (5.0/3.0)*this->compute_electron_mobility(a_energy, a_N);
}
Real air9eed::compute_electron_mobility(const Real a_energy, const Real a_N) const {
  return m_e_mobility.get_entry(a_energy);
}
Real air9eed::compute_N2plus_mobility(const Real a_EbyN) const   {return 2.E-4;}
Real air9eed::compute_N4plus_mobility(const Real a_EbyN) const   {return 2.E-4;}
Real air9eed::compute_O2plus_mobility(const Real a_EbyN) const   {return 2.E-4;}
Real air9eed::compute_O4plus_mobility(const Real a_EbyN) const   {return 2.E-4;}
Real air9eed::compute_O2plusN2_mobility(const Real a_EbyN) const {return 2.E-4;}
Real air9eed::compute_O2minus_mobility(const Real a_EbyN) const  {return 2.E-4;}
Real air9eed::compute_Ominus_mobility(const Real a_EbyN) const   {return 2.E-4;}

Real air9eed::compute_eed_diffco(const Real a_energy, const Real a_N) const{
  return (5.0/3.0)*this->compute_electron_diffco(a_energy, a_N);
}
Real air9eed::compute_electron_diffco(const Real a_energy, const Real a_N) const{
  return (2.0/3.0)*a_energy*this->compute_electron_mobility(a_energy, a_N);
}
Real air9eed::compute_N2plus_diffco() const {return 0.0;}
Real air9eed::compute_N4plus_diffco() const {return 0.0;}
Real air9eed::compute_O2plus_diffco() const {return 0.0;}
Real air9eed::compute_O4plus_diffco() const {return 0.0;}
Real air9eed::compute_O2plusN2_diffco() const {return 0.0;}
Real air9eed::compute_O2minus_diffco() const {return 0.0;}
Real air9eed::compute_Ominus_diffco() const {return 0.0;}

Real air9eed::compute_electron_N2_impact_ionization(const Real a_energy, const Real a_N) const {
  return m_e_N2_alpha.get_entry(a_energy);
}
Real air9eed::compute_electron_O2_impact_ionization(const Real a_energy, const Real a_N) const {
  return m_e_O2_alpha.get_entry(a_energy);
}
Real air9eed::compute_N2plus_N2_M_to_N4plus_M() const {return 5.E-41;}
Real air9eed::compute_N4plus_O2_to_O2_2N2() const {return 2.5E-16;}
Real air9eed::compute_N2plus_O2_to_O2plus_N2(const Real a_Tg) const {return 1.05E-15/sqrt(a_Tg);}
Real air9eed::compute_O2plus_2N2_to_O2plusN2_N2(const Real a_Tg) const {return 8.1E-38/(a_Tg*a_Tg);}
Real air9eed::compute_O2plusN2_N2_to_O2plus_2N2(const Real a_Tg) const {return 14.8*pow(a_Tg, -5.3)*exp(-2357.0/a_Tg);}
Real air9eed::compute_O2plusN2_O2_to_O4plus_N2() const {return 1.E-15;}
Real air9eed::compute_O2plus_O2_M_to_O4plus_M(const Real a_Tg) const {return 2.03E-34*pow(a_Tg, -3.2);}
Real air9eed::compute_e_O4plus_to_2O2(const Real a_Te) const {return 2.42E-11/(sqrt(a_Te));}
Real air9eed::compute_e_O2plus_to_O2(const Real a_Te) const {return 6.E-11/a_Te;}
Real air9eed::compute_e_2O2_to_O2minus_O2(const Real a_Te) const {return 0.0;}//6E-39/a_Te;}
Real air9eed::compute_O2minus_O4plus_to_3O2() const {return 1.E-13;}
Real air9eed::compute_O2minus_O4plus_M_to_3O2_M(const Real a_Tg) const {return 3.12E-31*pow(a_Tg, -2.5);}
Real air9eed::compute_O2minus_O2plus_M_to_2O2_M(const Real a_Tg) const {return 3.12E-31*pow(a_Tg, -2.5);}
Real air9eed::compute_e_O2_to_e_2O_c1(const Real a_energy, const Real a_N) const {
  Real k = 0.0;
  return 0.0;
    
  const Real min_energy       = 1.0;
  const Real max_energy       = 30.;
  const Real min_energy_coeff = 0.6937E-18;
  const Real max_energy_coeff = 0.8357E-14;

  if(a_energy < min_energy) { // Outside lower end
    k = min_energy_coeff;
  }
  else if(a_energy > max_energy){
    k = max_energy_coeff;
  }
  else {
    const Real A = -29.04;
    const Real B =  0.8249;
    const Real C = -17.40;
    const Real D =  7.278;
    const Real E = -2.656;

    const Real x = a_energy;
    k = exp(A + B*log(x) + C/x + D/(x*x) + E/(x*x*x));
  }

  return k;
} // TABLE
Real air9eed::compute_e_O2_to_e_2O_c2(const Real a_energy, const Real a_N) const {
  Real k = 0.0; return k;
    
  const Real min_energy       = 1.0;
  const Real max_energy       = 30.;
  const Real min_energy_coeff = 0.1149E-17;
  const Real max_energy_coeff = 0.2596E-12;

  if(a_energy < min_energy) { // Outside lower end
    k = min_energy_coeff;
  }
  else if(a_energy > max_energy){
    k = max_energy_coeff;
  }
  else {
    const Real A = -29.54;
    const Real B =  0.2660;
    const Real C = -10.13;
    const Real D = -2.733;
    const Real E =  1.101;

    const Real x = a_energy;
    k = exp(A + B*log(x) + C/x + D/(x*x) + E/(x*x*x));
  }

  return k;
} // TABLE
Real air9eed::compute_e_O2_to_Ominus_O(const Real a_energy, const Real a_N) const {
  Real k = 0.0;return k;
  
  const Real min_energy       = 1.0;
  const Real max_energy       = 30.;
  const Real min_energy_coeff = 0.1456E-17;
  const Real max_energy_coeff = 0.3192E-15;

  if(a_energy < min_energy) { // Outside lower end
    k = min_energy_coeff;
  }
  else if(a_energy > max_energy){
    k = max_energy_coeff;
  }
  else {
    const Real A = -37.26;
    const Real B =  0.4130;
    const Real C =  5.960;
    const Real D = -19.30;
    const Real E =  9.528;

    const Real x = a_energy;
    k = exp(A + B*log(x) + C/x + D/(x*x) + E/(x*x*x));
  }

  return k;
} // TABLE
Real air9eed::compute_Oplus_O2_to_O_O2(const Real a_Tg) const {return 3.46E-12/sqrt(a_Tg);}
Real air9eed::compute_e_N2_to_e_N2(const Real a_energy, const Real a_N) const {

  Real k = 0.0;return k;
  
  const Real min_energy       = 1.0;
  const Real max_energy       = 30.;
  const Real min_energy_coeff = 0.4298e-17;
  const Real max_energy_coeff = 0.3746e-15;

  if(a_energy < min_energy) { // Outside lower end
    k = min_energy_coeff;
  }
  else if(a_energy > max_energy){
    k = max_energy_coeff;
  }
  else {
    const Real A = -39.56;
    const Real B =  1.204;
    const Real C = -1.689;
    const Real D =  3.688;
    const Real E = -2.424;

    const Real x = a_energy;
    k = exp(A + B*log(x) + C/x + D/(x*x) + E/(x*x*x));
  }

  return k;
} // TABLE
Real air9eed::compute_e_O2_to_e_O2(const Real a_energy, const Real a_N) const {
  Real k = 0.0;return k;
  
  const Real min_energy       = 1.0;
  const Real max_energy       = 30.;
  const Real min_energy_coeff = 0.1561E-17;
  const Real max_energy_coeff = 0.3645E-15;

  if(a_energy < min_energy) { // Outside lower end
    k = min_energy_coeff;
  }
  else if(a_energy > max_energy){
    k = max_energy_coeff;
  }
  else {
    const Real A = -40.17;
    const Real B =  1.390;
    const Real C = -3.345;
    const Real D =  4.603;
    const Real E = -2.090;

    const Real x = a_energy;
    k = exp(A + B*log(x) + C/x + D/(x*x) + E/(x*x*x));
  }

  return k;
} // TABLE

Real air9eed::compute_e_N2_ionization_loss() const {
  return 15.6;
}
Real air9eed::compute_e_O2_ionization_loss() const {
  return 12.07;
}
Real air9eed::compute_e_O2_dissociation_loss_c1() const {
  return 5.58;
}
Real air9eed::compute_e_O2_dissociation_loss_c2() const {
  return 8.4;
}
Real air9eed::compute_e_O2_dissociative_attachment_loss() const {
  return 3.6;
}
Real air9eed::compute_e_O2_scattering_loss() const {
  return 1;
}
Real air9eed::compute_e_N2_scattering_loss() const {
  return 1;
}

Real air9eed::init_eed(const RealVect a_pos, const Real a_time, const RealVect a_E){
  const Real EbyN = (a_E/(m_N*units::s_Td)).vectorLength();
  return m_init_eed.direct_lookup(EbyN)*(m_species[m_electron_idx]->initial_data(a_pos, a_time));
}
