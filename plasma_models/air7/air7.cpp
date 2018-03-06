/*!
  @file   air7.cpp
  @brief  Implementation of air7.H
  @author Robert Marskar
  @date   Feb. 2018
  @todo   Implement corrections for the k1/k2 ionization source terms due to diffusion into high-field regions. 
*/

#include "air7.H"
#include "air7_species.H"
#include "units.H"

#include <ParmParse.H>
#include <PolyGeom.H>

air7::air7(){
  
  m_num_species = 7;
  m_num_photons = 3;

  air7::get_gas_parameters(m_Tg, m_p, m_N, m_O2frac, m_N2frac);

  { // Emission coefficients at boundaries. Can be overridden from input script.
    m_townsend2_electrode           = 1.E-2;
    m_townsend2_dielectric          = 1.E-2;
    m_electrode_quantum_efficiency  = 1.E-1;
    m_dielectric_quantum_efficiency = 1.E-1;
    m_photoionization_efficiency    = 0.1;
    m_excitation_efficiency         = 0.6;

    ParmParse pp("air7");
    pp.query("electrode_townsend2"       ,    m_townsend2_electrode);
    pp.query("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
    pp.query("dielectric_townsend2"       ,   m_townsend2_dielectric);
    pp.query("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
    pp.query("photoionization_efficiency",    m_photoionization_efficiency);
    pp.query("excitation_efficiency",         m_excitation_efficiency);
  }

  { // Quenching pressure
    m_pq        = 0.03947;
    ParmParse pp("air7");
    pp.query("quenching_pressure", m_pq);

    m_pq *= units::s_atm2pascal;
  }

  { // Ion mobility
    m_ion_mobility = 2.E-4;
    ParmParse pp("air7");
    pp.query("ion_mobility", m_ion_mobility);
  }

  // Instantiate species
  m_species.resize(m_num_species);
  m_electron_idx = 0;
  m_N2plus_idx   = 1;
  m_O2plus_idx   = 2;
  m_N4plus_idx   = 3;
  m_O4plus_idx   = 4;
  m_O2plusN2_idx = 5;
  m_O2minus_idx  = 6;
  m_species[m_electron_idx] = RefCountedPtr<species> (new air7::electron());
  m_species[m_N2plus_idx]   = RefCountedPtr<species> (new air7::N2plus());
  m_species[m_O2plus_idx]   = RefCountedPtr<species> (new air7::O2plus());
  m_species[m_N4plus_idx]   = RefCountedPtr<species> (new air7::N4plus());
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
}

air7::~air7(){

}

void air7::get_gas_parameters(Real& a_Tg, Real& a_p, Real& a_N, Real& a_O2frac, Real& a_N2frac){
    ParmParse pp("air7");
    a_p = 1.0;
    pp.query("gas_pressure", a_p); // Only get to adjust pressure for now.
    
    //    pp.query("gas_temperature", a_Tg);
    //    pp.query("gas_O2_frac", a_O2frac);
    //    pp.query("gas_N2_frac", a_N2frac);

    a_Tg = 300.;
    a_O2frac = 0.21;
    a_N2frac = 0.79; 

    const Real tot_frac = a_O2frac + a_N2frac; 
    a_p      = a_p*units::s_atm2pascal;
    a_O2frac = a_O2frac/tot_frac; // Normalize to one
    a_N2frac = a_N2frac/tot_frac;
    a_N      = a_p*units::s_Na/(a_Tg*units::s_R);
}

Vector<Real> air7::compute_cdr_diffusion_coefficients(const Real&         a_time,
						      const RealVect&     a_pos,
						      const RealVect&     a_E,
						      const Vector<Real>& a_cdr_densities) const {

  Vector<Real> diffco(m_num_species, 0.0);

  const Real ET = a_E.vectorLength()/(m_N*units::s_Td);
  diffco[m_electron_idx] = this->compute_electron_diffusion(ET);

  return diffco;
}

Vector<RealVect> air7::compute_cdr_velocities(const Real&         a_time,
					      const RealVect&     a_pos,
					      const RealVect&     a_E,
					      const Vector<Real>& a_cdr_densities) const{
  Vector<RealVect> velo(m_num_species);

  const Real ET = a_E.vectorLength()/(m_N*units::s_Td);
  
  velo[m_electron_idx] = -1.0*this->compute_electron_mobility(ET)*a_E;
  velo[m_N2plus_idx]   =  1.0*m_ion_mobility*a_E;
  velo[m_O2plus_idx]   =  1.0*m_ion_mobility*a_E;
  velo[m_N4plus_idx]   =  1.0*m_ion_mobility*a_E;
  velo[m_O4plus_idx]   =  1.0*m_ion_mobility*a_E;
  velo[m_O2plusN2_idx] =  1.0*m_ion_mobility*a_E;
  velo[m_O2minus_idx]  = -1.0*m_ion_mobility*a_E;

  return velo;
}

Vector<Real> air7::compute_cdr_source_terms(const Real              a_time,
					    const RealVect&         a_pos,
					    const RealVect&         a_E,
					    const RealVect&         a_gradE,
					    const Vector<Real>&     a_cdr_densities,
					    const Vector<Real>&     a_rte_densities,
					    const Vector<RealVect>& a_grad_cdr) const {
  Vector<Real> source(m_num_species, 0.0);

  const Real ET       = a_E.vectorLength()/(m_N*units::s_Td);
  const Real Te       = this->compute_electron_temperature(ET);
  const RealVect vele = -1.0*this->compute_electron_mobility(ET)*(a_E);
  const Real De       = this->compute_electron_diffusion(ET);

  const Real k1  = this->compute_townsend_ionization_N2(ET);
  const Real k2  = this->compute_townsend_ionization_O2(ET);
  const Real k3  = this->compute_N2plus_N2_M_to_N4plus_M();
  const Real k4  = this->compute_N4plus_O2_to_O2plus_2N2();
  const Real k5  = this->compute_N2plus_O2_to_O2plus_N2();
  const Real k6  = this->compute_O2plus_2N2_to_O2plusN2_N2();
  const Real k7  = this->compute_O2plusN2_N2_to_O2plus_2N2();
  const Real k8  = this->compute_O2plusN2_O2_to_O4plus_N2();
  const Real k9  = this->compute_O2plus_O2_M_to_O4plus_M();
  const Real k10 = this->compute_e_O4plus_to_2O2(Te);
  const Real k11 = this->compute_e_O2plus_to_2O(Te);
  const Real k12 = this->compute_e_2O2_to_O2minus_O2(Te);
  const Real k13 = this->compute_O2minus_O4plus_to_3O2();
  const Real k14 = this->compute_O2minus_O4plus_M_to_3O2_M();
  const Real k15 = this->compute_O2minus_O2plus_M_to_2O2_M();


  const Real n_N2    = m_N*m_N2frac;
  const Real n_O2    = m_N*m_O2frac;
  const Real n_Ne    = a_cdr_densities[m_electron_idx];
  const Real n_N2p   = a_cdr_densities[m_N2plus_idx];
  const Real n_O2p   = a_cdr_densities[m_O2plus_idx];
  const Real n_N4p   = a_cdr_densities[m_N4plus_idx];
  const Real n_O4p   = a_cdr_densities[m_O4plus_idx];
  const Real n_O2pN2 = a_cdr_densities[m_O2plusN2_idx];
  const Real n_O2m   = a_cdr_densities[m_O2minus_idx];
  const Real n_M     = n_N2 + n_O2;

  Real S = 0.0;

  const Real a1 = PolyGeom::dot(a_E, De*a_grad_cdr[m_electron_idx]);
  const Real a2 = PolyGeom::dot(a_E, (1.0 + a_cdr_densities[m_electron_idx])*vele);
  const Real corr = (1.0 - a1/a2);
  
  // k1 reaction, e + N2 -> e + e + N2+
  S = k1 * n_Ne * n_N2;
  source[m_electron_idx] += S;
  source[m_N2plus_idx]   += S;

  // k2 reaction, e + O2 -> e + e + O2+
  S = k2 * n_Ne * n_O2;
  source[m_electron_idx] += S;
  source[m_O2plus_idx]   += S;

  // k3 reaction
  S = k3 * n_N2p * n_N2 * n_M;
  source[m_N2plus_idx] -= S;
  source[m_N4plus_idx] += S;

  // k4 reaction
  S = k4 * n_N4p * n_O2;
  source[m_N4plus_idx] -= S;
  source[m_O2plus_idx] += S;

#if 0
  // k5 reaction
  S = k5 * n_N2p * n_O2;
  source[m_N2plus_idx] -= S;
  source[m_O2plus_idx] -= S;

  // k6 reaction
  S = k6 * n_O2p * n_N2 * n_N2;
  source[m_O2plus_idx]   -= S;
  source[m_O2plusN2_idx] += S;

  // k7 reaction
  S = k7 * n_O2pN2 * n_N2;
  source[m_O2plusN2_idx] -= S;
  source[m_O2plus_idx]   += S;

  // k8 reaction
  S = k8 * n_O2pN2 * n_O2;
  source[m_O2plusN2_idx] -= S;
  source[m_O4plus_idx]   += S;

  // k9 reaction
  S = k9 * n_O2p * n_O2 * n_M;
  source[m_O2plus_idx] -= S;
  source[m_O4plus_idx] += S;

  // k10 reaction
  S = k10 * n_Ne * n_O4p;
  source[m_electron_idx] -= S;
  source[m_O4plus_idx]   -= S;
#endif
  
  // k11 reaction
  S = k11 * n_Ne * n_O2p;
  source[m_electron_idx] -= S;
  source[m_O2plus_idx]   -= S;

  // k12 reaction
  S = k12 * n_Ne * n_O2 * n_O2;
  source[m_electron_idx] -= S;
  source[m_O2minus_idx]  += S;

  // k13 reaction
  S = n_O2m * n_O4p;
  source[m_O2minus_idx] -= S;
  source[m_O4plus_idx]  -= S;

  // k14 reaction
  S = k14 * n_O2m * n_O4p * n_M;
  source[m_O2minus_idx] -= S;
  source[m_O4plus_idx]  -= S;

  // k15 reaction
  S = k15 * n_O2m * n_O2p * n_M;
  source[m_O2minus_idx] -= S;
  source[m_O2plus_idx]  -= S;
  
  
  // Photoionization source term
  const air7::photon_one*   photon1 = static_cast<air7::photon_one*>   (&(*m_photons[m_photon1_idx]));
  const air7::photon_two*   photon2 = static_cast<air7::photon_two*>   (&(*m_photons[m_photon2_idx]));
  const air7::photon_three* photon3 = static_cast<air7::photon_three*> (&(*m_photons[m_photon3_idx]));

  const Real Sph = m_photoionization_efficiency*units::s_c0*m_O2frac*m_p*(photon1->get_A()*a_rte_densities[m_photon1_idx]
  							         	  + photon2->get_A()*a_rte_densities[m_photon2_idx]
  									  + photon3->get_A()*a_rte_densities[m_photon3_idx]);

  //  source[m_electron_idx] += Sph;
  //  source[m_O2plus_idx]   += Sph;

  return source;
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
  return this->compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients,
				  a_rte_fluxes, a_extrap_cdr_fluxes, m_townsend2_electrode, m_electrode_quantum_efficiency);
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

  return this->compute_cdr_fluxes(a_time, a_pos, a_normal, a_E, a_cdr_densities, a_cdr_velocities, a_cdr_gradients,
				  a_rte_fluxes, a_extrap_cdr_fluxes, m_townsend2_dielectric, m_dielectric_quantum_efficiency);
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
				      const Real&         a_townsend,
				      const Real&         a_quantum_efficiency) const {

  Vector<Real> fluxes(m_num_species, 0.0);

  const Real ET  = a_E.vectorLength()/(m_N*units::s_Td);;
  const Real Te  = this->compute_electron_temperature(ET);
  const Real Tg  = m_Tg;
  const Real mO2 = 2.65E-26;

  const bool anode   = PolyGeom::dot(a_E, a_normal) > 0.0;
  const bool cathode = PolyGeom::dot(a_E, a_normal) < 0.0;
  
  Vector<Real> aj(m_num_species, 0.0);
  if(anode){
    aj[m_electron_idx] = 1.0;
    aj[m_O2minus_idx]  = 1.0;
  }
  else if(cathode){
    aj[m_N2plus_idx]   = 1.0;
    aj[m_O2plus_idx]   = 1.0;
    aj[m_N4plus_idx]   = 1.0;
    aj[m_O4plus_idx]   = 1.0;
    aj[m_O2plusN2_idx] = 1.0;
  }

  // Drift outflow.
  for (int i = 0; i < m_num_species; i++){
    fluxes[i] = aj[i]*a_extrap_cdr_fluxes[i];
  }

  // Thermal outflow
  const Real vth_g = sqrt(8.0*units::s_kb*Te/(units::s_pi*mO2));
  const Real vth_e = sqrt(8.0*units::s_kb*Te/(units::s_pi*units::s_me));
  for (int i = 0; i < m_num_species; i++){
    Real vth = 0.0;
    if(i == m_electron_idx){
      vth = vth_e;
    }
    else{
      vth = vth_g;
    }
    fluxes[i] += 0.25*vth*a_cdr_densities[i];
  }


#if 1 // Secondary emission
  if(cathode){
    Real ion_bombardment    = 0.0;
    Real photon_bombardment = 0.0;
    
    ion_bombardment += fluxes[m_N2plus_idx];
    ion_bombardment += fluxes[m_O2plus_idx];
    ion_bombardment += fluxes[m_N4plus_idx];
    ion_bombardment += fluxes[m_O4plus_idx];
    ion_bombardment += fluxes[m_O2plusN2_idx];

    photon_bombardment += a_rte_fluxes[m_photon1_idx];
    photon_bombardment += a_rte_fluxes[m_photon2_idx];
    photon_bombardment += a_rte_fluxes[m_photon3_idx];

    fluxes[m_electron_idx] -= ion_bombardment*a_townsend;
    fluxes[m_electron_idx] -= photon_bombardment*a_quantum_efficiency;
  }
#endif
  
  return fluxes;
}

Vector<Real> air7::compute_rte_source_terms(const Real&         a_time,
					    const RealVect&     a_pos,
					    const RealVect&     a_E,
					    const Vector<Real>& a_cdr_densities) const {
  Vector<Real> ret(m_num_photons);

  const Real ET      = a_E.vectorLength()/(units::s_Td*m_N);
  const RealVect vel = -1.0*this->compute_electron_mobility(ET)*a_E;
  const Real alpha   = this->compute_townsend_ionization_N2(ET);
  const Real N2      = m_N*m_N2frac;
  const Real Ne      = a_cdr_densities[m_electron_idx];       // Electron density
  const Real ve      = vel.vectorLength();
  const Real Se      = Max(0., alpha*N2*Ne);                  // Excitations = alpha*Ne*ve

  // Photo emissions = electron excitations * efficiency * quenching
  ret[m_photon1_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));
  ret[m_photon2_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));
  ret[m_photon3_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));

  return ret;
}

Real air7::initial_sigma(const Real a_time, const RealVect& a_pos) const {
  return 0.0;
}

Real air7::compute_electron_temperature(const Real a_EbyN) const {
  Real temp = 0.0;

  const Real safety = 1.0; // To avoid division by zero
  const Real minE   = 10;
  const Real maxE   = 3000;
  const Real min_eV = 0.9632;
  const Real max_eV = 42.19;

  if(a_EbyN < minE){
    temp = min_eV;
  }
  else if(a_EbyN > maxE){
    temp = max_eV;
  }
  else {
    const Real A = -2.689;
    const Real B =  0.8036;
    const Real C = -8.165;
    const Real D =  338.6;
    const Real E = -1768;

    const Real x = a_EbyN;
    temp = exp(A + B*log(x) + C/x + D/(x*x) + E/(x*x*x));
  }

  // temp is in energy so far, make it into Kelvin by E = 1.5*k_b*T
  temp *= units::s_Qe;          // eV -> Joule
  temp *= 1./(1.5*units::s_kb); // Mean temperature

  return safety + temp;
}

Real air7::compute_electron_mobility(const Real a_EbyN) const {

  Real mobility = 0.0;
  
  const Real minE = 10.0;
  const Real maxE = 3000.;
  const Real min_mob = 0.311E25;
  const Real max_mob = 0.7172E24;
  
  if(a_EbyN < minE){
    mobility = min_mob;
  }
  else if(a_EbyN > maxE){
    mobility = max_mob;
  }
  else {
    const Real A =  57.39;
    const Real B = -0.3072;
    const Real C = -15.73;
    const Real D =  298.4;
    const Real E = -1701.0;

    const Real x = a_EbyN;
    mobility = exp(A + B*log(x) + C/x + D/(x*x) + E/(x*x*x));
  }

  return mobility/m_N;
}

Real air7::compute_electron_diffusion(const Real a_EbyN) const {
  Real diffCo = 0.0;
  
  const Real minE = 10.0;
  const Real maxE = 3000.;
  const Real min_diffco = 0.1887E25;
  const Real max_diffco = 0.1318E26;
  
  if(a_EbyN < minE){
    diffCo = min_diffco;
  }
  else if(a_EbyN > maxE){
    diffCo = max_diffco;
  }
  else {
    const Real A =  54.05;
    const Real B =  0.4732;
    const Real C = -0.5350;
    const Real D =  246.8;
    const Real E = -1659.0;

    const Real x = a_EbyN;
    diffCo = exp(A + B*log(x) + C/x + D/(x*x) + E/(x*x*x));
  }

  return diffCo/m_N;
}

Real air7::compute_townsend_ionization_N2(const Real a_EbyN) const {
  Real k = 0.0;
  
  const Real minE = 10.0;
  const Real maxE = 3000.;
  const Real min_k = 0.4537E-25;
  const Real max_k = 0.6042E-13;
  
  if(a_EbyN < minE){
    k = min_k;
  }
  else if(a_EbyN > maxE){
    k = max_k;
  }
  else {
    const Real A = -38.34;
    const Real B =  1.023;
    const Real C = -888.1;
    const Real D =  8400;
    const Real E = -0.1866E6;

    const Real x = a_EbyN;
    k = exp(A + B*log(x) + C/x + D/(x*x) + E/(x*x*x));
  }

  return k;
}

Real air7::compute_townsend_ionization_O2(const Real a_EbyN) const {
  Real k = 0.0;
  
  const Real minE = 10.0;
  const Real maxE = 3000.;
  const Real min_k = 0.3627E-47;
  const Real max_k = 0.6268E-13;
  
  if(a_EbyN < minE){
    k = min_k;
  }
  else if(a_EbyN > maxE){
    k = max_k;
  }
  else {
    const Real A = -40.70;
    const Real B =  1.309;
    const Real C = -543.6;
    const Real D =  2091;
    const Real E = -0.3810E6;

    const Real x = a_EbyN;
    k = exp(A + B*log(x) + C/x + D/(x*x) + E/(x*x*x));
  }

  return k;
}

Real air7::compute_N2plus_N2_M_to_N4plus_M() const {
  return 5.0E-41;
}

Real air7::compute_N4plus_O2_to_O2plus_2N2() const {
  return 2.5E-16;
}

Real air7::compute_N2plus_O2_to_O2plus_N2() const {
  return 6.0E-17;
}

Real air7::compute_O2plus_2N2_to_O2plusN2_N2() const {
  return 9.0E-43;
}

Real air7::compute_O2plusN2_N2_to_O2plus_2N2() const {
  return 4.3E-16;
}

Real air7::compute_O2plusN2_O2_to_O4plus_N2() const {
  return 1.E-15;
}

Real air7::compute_O2plus_O2_M_to_O4plus_M() const {
  return 2.4E-42;
}

Real air7::compute_e_O4plus_to_2O2(const Real a_Te) const {
  return 1.4E-12*sqrt(300./a_Te);
}

Real air7::compute_e_O2plus_to_2O(const Real a_Te) const {
  return 6E-11/(a_Te);
}

Real air7::compute_e_2O2_to_O2minus_O2(const Real a_Te) const {
  return 6E-39/(a_Te);
}

Real air7::compute_O2minus_O4plus_to_3O2() const {
  return 1.E-13;
}

Real air7::compute_O2minus_O4plus_M_to_3O2_M() const {
  return 2.E-37;
}

Real air7::compute_O2minus_O2plus_M_to_2O2_M() const {
  return 2.E-37;
}
