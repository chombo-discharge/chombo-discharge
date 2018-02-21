/*!
  @file   air_11eed.cpp
  @brief  Implementation of air_11eed.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air_11eed.H"
#include "air_11eed_species.H"
#include "units.H"

#include <ParmParse.H>

air_11eed::air_11eed(){

  m_num_species = 12; // 11 reactive ones plus the eed
  m_num_photons = 3;  // Bourdon model for photons

 
  air_11eed::get_gas_parameters(m_Tg, m_p, m_N, m_O2frac, m_N2frac); // Get gas parameters


  { // Emission coefficients at boundaries. Can be overridden from input script.
    m_townsend2_electrode           = 1.E-3;
    m_townsend2_dielectric          = 1.E-6;
    m_electrode_quantum_efficiency  = 1.E-2;
    m_dielectric_quantum_efficiency = 1.E-4;

    ParmParse pp("air_11eed");
    pp.query("electrode_townsend2"       ,    m_townsend2_electrode);
    pp.query("dielectric_townsend2"       ,   m_townsend2_dielectric);
    pp.query("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
    pp.query("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
  }

  
  // Instantiate cdr species
  m_species.resize(m_num_species);
  m_eed_idx      = 0;
  m_electron_idx = 1;
  m_N2_idx       = 2;
  m_O2_idx       = 3;
  m_N2plus_idx   = 4;
  m_N4plus_idx   = 5;
  m_O2plus_idx   = 6;
  m_O4plus_idx   = 7;
  m_O2plusN2_idx = 8;
  m_O2minus_idx  = 9;
  m_Ominus_idx   = 10;
  m_O_idx        = 11;
  
  m_species[m_eed_idx]      = RefCountedPtr<species> (new air_11eed::eed());
  m_species[m_electron_idx] = RefCountedPtr<species> (new air_11eed::electron());
  m_species[m_N2_idx]       = RefCountedPtr<species> (new air_11eed::N2());
  m_species[m_O2_idx]       = RefCountedPtr<species> (new air_11eed::O2());
  m_species[m_N2plus_idx]   = RefCountedPtr<species> (new air_11eed::N2plus());
  m_species[m_N4plus_idx]   = RefCountedPtr<species> (new air_11eed::N4plus());
  m_species[m_O2plus_idx]   = RefCountedPtr<species> (new air_11eed::O2plus());
  m_species[m_O4plus_idx]   = RefCountedPtr<species> (new air_11eed::O4plus());
  m_species[m_O2plusN2_idx] = RefCountedPtr<species> (new air_11eed::O2plusN2());
  m_species[m_O2minus_idx]  = RefCountedPtr<species> (new air_11eed::O2minus());
  m_species[m_Ominus_idx]   = RefCountedPtr<species> (new air_11eed::Ominus());
  m_species[m_O_idx]        = RefCountedPtr<species> (new air_11eed::O());


  // Instantiate photon solvers
  m_photons.resize(m_num_photons);
  m_photon1_idx = 0;
  m_photon2_idx = 1;
  m_photon3_idx = 2;
  m_photons[m_photon1_idx] = RefCountedPtr<photon_group> (new air_11eed::photon_one());
  m_photons[m_photon2_idx] = RefCountedPtr<photon_group> (new air_11eed::photon_two());
  m_photons[m_photon3_idx] = RefCountedPtr<photon_group> (new air_11eed::photon_three());
  
}

air_11eed::~air_11eed(){

}

void air_11eed::get_gas_parameters(Real& a_Tg, Real& a_p, Real& a_N, Real& a_O2frac, Real& a_N2frac){
    ParmParse pp("air_11eed");
    pp.get("gas_temperature", a_Tg);
    pp.get("gas_pressure", a_p);
    pp.get("gas_O2_frac", a_O2frac);
    pp.get("gas_N2_frac", a_N2frac);

    const Real tot_frac = a_O2frac + a_N2frac; 
    a_p      = a_p*units::s_atm2pascal;
    a_O2frac = a_O2frac/tot_frac; // Normalize to one
    a_N2frac = a_N2frac/tot_frac;
    a_N      = a_p*units::s_Na/(a_Tg*units::s_R);
}

Vector<Real> air_11eed::compute_cdr_diffusion_coefficients(const Real&         a_time,
							   const RealVect&     a_pos,
							   const RealVect&     a_E,
							   const Vector<Real>& a_cdr_densities) const {

  return Vector<Real>(m_num_species, 0.0);
}

Vector<RealVect> air_11eed::compute_cdr_velocities(const Real&         a_time,
						   const RealVect&     a_pos,
						   const RealVect&     a_E,
						   const Vector<Real>& a_cdr_densities) const {
  Vector<RealVect> velocities(m_num_species, RealVect::Zero);


  const Real electron_energy = a_cdr_densities[m_eed_idx]/(1.E-20 + a_cdr_densities[m_electron_idx]);
  const Real N               = a_cdr_densities[m_O2_idx] + a_cdr_densities[m_N2_idx];
  const Real EbyN            = (a_E/N*units::s_Td).vectorLength();

  velocities[m_eed_idx]      = (5.0/3.0)*this->compute_electron_mobility(electron_energy, N)*a_E;
  velocities[m_electron_idx] = -1.0*this->compute_electron_mobility(electron_energy, N)*a_E;
  velocities[m_N2_idx]       = RealVect::Zero;
  velocities[m_O2_idx]       = RealVect::Zero;
  velocities[m_N2plus_idx]   = this->compute_N2plus_mobility(EbyN)*a_E;
  velocities[m_N4plus_idx]   = this->compute_N4plus_mobility(EbyN)*a_E;
  velocities[m_O2plus_idx]   = this->compute_O2plus_mobility(EbyN)*a_E;
  velocities[m_O4plus_idx]   = this->compute_O4plus_mobility(EbyN)*a_E;
  velocities[m_O2plusN2_idx] = this->compute_O2plusN2_mobility(EbyN)*a_E;
  velocities[m_O2minus_idx]  = -1.0*this->compute_O2minus_mobility(EbyN)*a_E;
  velocities[m_Ominus_idx]   = -1.0*this->compute_Ominus_mobility(EbyN)*a_E;

  return velocities;
}
  
Vector<Real> air_11eed::compute_cdr_source_terms(const Real              a_time,
						 const RealVect&         a_pos,
						 const RealVect&         a_E,
						 const RealVect&         a_gradE,
						 const Vector<Real>&     a_cdr_densities,
						 const Vector<Real>&     a_rte_densities,
						 const Vector<RealVect>& a_grad_cdr) const {
  return Vector<Real>(m_num_species, 0.0);
}

Vector<Real> air_11eed::compute_cdr_electrode_fluxes(const Real&         a_time,
						     const RealVect&     a_pos,
						     const RealVect&     a_normal,
						     const RealVect&     a_E,
						     const Vector<Real>& a_cdr_densities,
						     const Vector<Real>& a_cdr_velocities,
						     const Vector<Real>& a_rte_fluxes,
						     const Vector<Real>& a_extrap_cdr_fluxes) const {
  return Vector<Real>(m_num_species, 0.0);
}

Vector<Real> air_11eed::compute_cdr_dielectric_fluxes(const Real&         a_time,
						      const RealVect&     a_pos,
						      const RealVect&     a_normal,
						      const RealVect&     a_E,
						      const Vector<Real>& a_cdr_densities,
						      const Vector<Real>& a_cdr_velocities,
						      const Vector<Real>& a_rte_fluxes,
						      const Vector<Real>& a_extrap_cdr_fluxes) const {
  return Vector<Real>(m_num_species, 0.0);
}

Vector<Real> air_11eed::compute_rte_source_terms(const Real&         a_time,
						 const RealVect&     a_pos,
						 const RealVect&     a_E,
						 const Vector<Real>& a_cdr_densities) const {
  return Vector<Real>(m_num_photons, 0.0);
}

Real air_11eed::initial_sigma(const Real a_time, const RealVect& a_pos) const {
  return 0.0;
}

Real air_11eed::compute_electron_mobility(const Real a_energy, const Real a_N) const {

  Real mobility;

  // These are the valid ranges from the BOLSIG call
  const Real min_energy     = 1.E-2;
  const Real max_energy     = 20.;
  const Real min_energy_mob = 0.148E27;
  const Real max_energy_mob = 0.537E24;
  
  if(a_energy < min_energy){ // Outside lower end
    mobility = min_energy_mob;
  }
  else if(a_energy > max_energy){
    mobility = max_energy_mob;
  }
  else {
    const Real A =  55.93;
    const Real B = -0.3830;
    const Real C =  0.2688;
    const Real D = -0.1298E-1;
    const Real E = -0.1056E-3;

    const Real x = a_energy;
    mobility = exp(A + B*log(x) + C/x + D/(x*x) + E/(x*x*x));
  }

  mobility *= 1./a_N;

  return mobility;
}

Real air_11eed::compute_N2plus_mobility(const Real a_EbyN) const {
  return 2.E-4;
}

Real air_11eed::compute_N4plus_mobility(const Real a_EbyN) const {
  return 2.E-4;
}

Real air_11eed::compute_O2plus_mobility(const Real a_EbyN) const {
  return 2.E-4;
}

Real air_11eed::compute_O4plus_mobility(const Real a_EbyN) const {
  return 2.E-4;
}

Real air_11eed::compute_O2plusN2_mobility(const Real a_EbyN) const {
  return 2.E-4;
}

Real air_11eed::compute_O2minus_mobility(const Real a_EbyN) const {
  return 2.E-4;
}

Real air_11eed::compute_Ominus_mobility(const Real a_EbyN) const {
  return 2.E-4;
}
