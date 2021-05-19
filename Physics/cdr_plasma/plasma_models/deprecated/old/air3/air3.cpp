/*!
  @file   air3.cpp
  @brief  Implementation of air3.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air3.H"
#include "air3_species.H"
#include "units.H"

#include <ParmParse.H>
#include <PolyGeom.H>

air3::air3(){
  
  m_num_species = 3;
  m_num_Photons = 3;

  air3::get_gas_parameters(m_Tg, m_p, m_N, m_O2frac, m_N2frac);

  { // Emission coefficients at boundaries. Can be overridden from input script.
    m_townsend2_electrode           = 1.E-3;
    m_townsend2_dielectric          = 1.E-6;
    m_electrode_quantum_efficiency  = 1.E-2;
    m_dielectric_quantum_efficiency = 1.E-4;
    m_photoionization_efficiency    = 0.1;
    m_excitation_efficiency         = 0.6;

    ParmParse pp("air3");
    pp.query("electrode_townsend2"       ,    m_townsend2_electrode);
    pp.query("electrode_quantum_efficiency",  m_electrode_quantum_efficiency);
    pp.query("dielectric_townsend2"       ,   m_townsend2_dielectric);
    pp.query("dielectric_quantum_efficiency", m_dielectric_quantum_efficiency);
    pp.query("photoionization_efficiency",    m_photoionization_efficiency);
    pp.query("excitation_efficiency",         m_excitation_efficiency);
  }

  { // Quenching pressure
    m_pq        = 0.03947;
    ParmParse pp("air3");
    pp.query("quenching_pressure", m_pq);

    m_pq *= units::s_atm2pascal;
  }

  { // Ion mobility
    m_ion_mobility = 2.E-4;
    ParmParse pp("air3");
    pp.query("ion_mobility", m_ion_mobility);
  }

  { // Recombination coefficients and attachment
    m_electron_recombination = 5.E-14;
    m_ion_recombination      = 2.0E-12;
    m_electron_detachment    = 1.E-18;

    ParmParse pp("air3");

    pp.query("electron_recombination", m_electron_recombination);
    pp.query("ion_recombination", m_ion_recombination);
    pp.query("electron_detachment", m_electron_detachment);
  }

    

  // Instantiate species
  m_species.resize(m_num_species);
  m_electron_idx = 0;
  m_positive_idx = 1;
  m_negative_idx = 2;
  m_species[m_electron_idx] = RefCountedPtr<species> (new air3::electron());
  m_species[m_positive_idx] = RefCountedPtr<species> (new air3::positive_species());
  m_species[m_negative_idx] = RefCountedPtr<species> (new air3::negative_species());


  // Instantiate Photon solvers
  m_Photons.resize(m_num_Photons);
  m_Photon1_idx = 0;
  m_Photon2_idx = 1;
  m_Photon3_idx = 2;
  m_Photons[m_Photon1_idx] = RefCountedPtr<Photon_group> (new air3::Photon_one());
  m_Photons[m_Photon2_idx] = RefCountedPtr<Photon_group> (new air3::Photon_two());
  m_Photons[m_Photon3_idx] = RefCountedPtr<Photon_group> (new air3::Photon_three());
}

air3::~air3(){

}

void air3::get_gas_parameters(Real& a_Tg, Real& a_p, Real& a_N, Real& a_O2frac, Real& a_N2frac){
  ParmParse pp("air3");
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

Vector<Real> air3::compute_cdr_diffusion_coefficients(const Real&         a_time,
						      const RealVect&     a_pos,
						      const RealVect&     a_E,
						      const Vector<Real>& a_cdr_densities) const {

  Vector<Real> diffco(m_num_species, 0.0);

  const Real ET = a_E.vectorLength()/(m_N*units::s_Td);
  diffco[m_electron_idx] = this->compute_electron_diffusion(ET);

  return diffco;
}

Vector<RealVect> air3::compute_cdr_velocities(const Real&         a_time,
					      const RealVect&     a_pos,
					      const RealVect&     a_E,
					      const Vector<Real>& a_cdr_densities) const{
  Vector<RealVect> velo(m_num_species);

  const Real ET = a_E.vectorLength()/(m_N*units::s_Td);
  
  velo[m_electron_idx] = -1.0*this->compute_electron_mobility(ET)*a_E;
  velo[m_positive_idx] =  1.0*m_ion_mobility*a_E;
  velo[m_negative_idx] = -1.0*m_ion_mobility*a_E;

  return velo;
}

Vector<Real> air3::compute_cdr_source_terms(const Real              a_time,
					    const RealVect&         a_pos,
					    const RealVect&         a_E,
					    const RealVect&         a_gradE,
					    const Vector<Real>&     a_cdr_densities,
					    const Vector<Real>&     a_rte_densities,
					    const Vector<RealVect>& a_grad_cdr) const {
  Vector<Real> source(m_num_species, 0.0);

  const Real ET = a_E.vectorLength()/(m_N*units::s_Td);

  const RealVect vele = -1.0*this->compute_electron_mobility(ET)*a_E;
  const Real ve       = vele.vectorLength();
  const Real De       = this->compute_electron_diffusion(ET);
  const Real alpha    = this->compute_townsend_ionization(ET);
  const Real eta      = this->compute_townsend_attachment(ET);

  const Real bep  = m_electron_recombination;
  const Real bpn  = m_ion_recombination;
  const Real kdet = m_electron_detachment;

  const Real Ne = a_cdr_densities[m_electron_idx];
  const Real Np = a_cdr_densities[m_positive_idx];
  const Real Nn = a_cdr_densities[m_negative_idx];

  const air3::Photon_one*   Photon1 = static_cast<air3::Photon_one*>   (&(*m_Photons[m_Photon1_idx]));
  const air3::Photon_two*   Photon2 = static_cast<air3::Photon_two*>   (&(*m_Photons[m_Photon2_idx]));
  const air3::Photon_three* Photon3 = static_cast<air3::Photon_three*> (&(*m_Photons[m_Photon3_idx]));

  const Real Sph = m_photoionization_efficiency*units::s_c0*m_O2frac*m_p*(Photon1->get_A()*a_rte_densities[m_Photon1_idx]
							         	  + Photon2->get_A()*a_rte_densities[m_Photon2_idx]
									  + Photon3->get_A()*a_rte_densities[m_Photon3_idx]);

  const Real alpha_corr = alpha*(1 - PolyGeom::dot(a_E,De*a_grad_cdr[m_electron_idx])/((1.0 + Ne)*PolyGeom::dot(vele, a_E)));
  Real& Se = source[m_electron_idx];
  Real& Sp = source[m_positive_idx];
  Real& Sn = source[m_negative_idx];
  
  Se = alpha*Ne*ve - eta*Ne*ve - bep*Ne*Np + Sph + kdet*Nn*m_N;
  Sp = alpha*Ne*ve - bep*Ne*Np - bpn*Np*Nn + Sph;
  Sn = eta*Ne*ve   - bpn*Np*Nn - kdet*Nn*m_N;

  return source;
}

Vector<Real> air3::compute_cdr_electrode_fluxes(const Real&         a_time,
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

Vector<Real> air3::compute_cdr_dielectric_fluxes(const Real&         a_time,
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

Vector<Real> air3::compute_cdr_fluxes(const Real&         a_time,
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

  const bool anode   = PolyGeom::dot(a_E, a_normal) > 0.0;
  const bool cathode = PolyGeom::dot(a_E, a_normal) < 0.0;
  
  Vector<Real> aj(m_num_species, 0.0);
  if(anode){
    aj[m_electron_idx] = 1.0;
    aj[m_negative_idx] = 1.0;
  }
  else if(cathode){
    aj[m_positive_idx] = 1.0;
  }

  // Drift outflow.
  for (int i = 0; i < m_num_species; i++){
    fluxes[i] = aj[i]*a_extrap_cdr_fluxes[i];
  }

  // Secondary emission for cathode surfaces due to ion and Photon bombardment
  if(cathode){
    Real ion_bombardment    = 0.0;
    Real Photon_bombardment = 0.0;
    
    ion_bombardment += fluxes[m_positive_idx];

    Photon_bombardment += a_rte_fluxes[m_Photon1_idx];
    Photon_bombardment += a_rte_fluxes[m_Photon2_idx];
    Photon_bombardment += a_rte_fluxes[m_Photon3_idx];

    fluxes[m_electron_idx] -= ion_bombardment*a_townsend;
    fluxes[m_electron_idx] -= Photon_bombardment*a_quantum_efficiency;
  }

  return fluxes;
}

Vector<Real> air3::compute_rte_source_terms(const Real&         a_time,
					    const RealVect&     a_pos,
					    const RealVect&     a_E,
					    const Vector<Real>& a_cdr_densities) const {
  Vector<Real> ret(m_num_Photons);

  const Real ET      = a_E.vectorLength()/(units::s_Td*m_N);
  const RealVect vel = -1.0*this->compute_electron_mobility(ET)*a_E;
  const Real alpha   = this->compute_townsend_ionization(ET);
  const Real Ne      = a_cdr_densities[m_electron_idx];       // Electron density
  const Real ve      = vel.vectorLength();
  const Real Se      = Max(0., alpha*Ne*ve);               // Excitations = alpha*Ne*ve

  // Photo emissions = electron excitations * efficiency * quenching
  ret[m_Photon1_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));
  ret[m_Photon2_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));
  ret[m_Photon3_idx] = Se*m_excitation_efficiency*(m_pq/(m_pq + m_p));

  return ret;
}

Real air3::initial_sigma(const Real a_time, const RealVect& a_pos) const {
  return 0.0;
}

Real air3::compute_electron_mobility(const Real a_EbyN) const {

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

Real air3::compute_electron_diffusion(const Real a_EbyN) const {
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

Real air3::compute_townsend_ionization(const Real a_EbyN) const {
  Real alpha = 0.0;
  
  const Real minE = 10.0;
  const Real maxE = 3000.;
  const Real min_alpha = 0.13E-63;
  const Real max_alpha = 0.4489E-19;
  
  if(a_EbyN < minE){
    alpha = min_alpha;
  }
  else if(a_EbyN > maxE){
    alpha = max_alpha;
  }
  else {
    const Real A = -48.01;
    const Real B =  0.4634;
    const Real C = -768.3;
    const Real D =  9753.0;
    const Real E = -0.1208E6;

    const Real x = a_EbyN;
    alpha = exp(A + B*log(x) + C/x + D/(x*x) + E/(x*x*x));
  }

  return alpha*m_N;
}

Real air3::compute_townsend_attachment(const Real a_EbyN) const {
  Real eta = 0.0;
  
  const Real minE = 10.0;
  const Real maxE = 3000.;
  const Real min_eta = 0.2684E-25;
  const Real max_eta = 0.4933E-23;
  
  if(a_EbyN < minE){
    eta = min_eta;
  }
  else if(a_EbyN > maxE){
    eta = max_eta;
  }
  else {
    const Real A = -46.47;
    const Real B = -0.8973;
    const Real C = -38.94;
    const Real D = -8975;
    const Real E =  0.8329E5;

    const Real x = a_EbyN;
    eta = exp(A + B*log(x) + C/x + D/(x*x) + E/(x*x*x));
  }

  return eta*m_N;
#include "CD_NamespaceFooter.H"
