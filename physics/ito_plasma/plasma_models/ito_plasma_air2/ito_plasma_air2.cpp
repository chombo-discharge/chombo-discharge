/*!
  @file   ito_plasma_air2.cpp
  @brief  Implementation of ito_plasma_air2.H
  @author Robert Marskar
  @date   June 2020
*/

#include "ito_plasma_air2.H"
#include "units.H"

#include <ParmParse.H>

using namespace physics::ito_plasma;

ito_plasma_air2::ito_plasma_air2(){
  m_num_ito_species = 2;
  m_num_rte_species = 1;

  ParmParse pp("ito_plasma_air2");
  Vector<Real> v;
  
  // Stuff for initial particles
  pp.get   ("seed",            m_seed);
  pp.get   ("blob_radius",     m_blob_radius);
  pp.get   ("num_particles",   m_num_particles);
  pp.get   ("particle_weight", m_particle_weight);
  pp.get   ("react_ppc",       m_ppc);
  pp.getarr("blob_center",     v, 0, SpaceDim); m_blob_center = RealVect(D_DECL(v[0], v[1], v[2]));

  // Reaction stuff
  pp.get("quenching_pressure", m_pq);
  pp.get("photoi_factor",      m_photoi_factor);
  pp.get("tau_switch",         m_tau_switch);
  pp.get("poisson_switch",     m_poisson_switch);

  // Standard air. 
  m_p = 1.0;
  m_T = 300;
  m_N2frac = 0.8;
  m_O2frac = 0.2;

  // Convert to SI units
  m_p  = m_p*units::s_atm2pascal;
  m_pq = m_pq*units::s_atm2pascal;
  m_N  = m_p*units::s_Na/(m_T*units::s_R);

  // Set up species
  m_ito_species.resize(m_num_ito_species);
  m_rte_species.resize(m_num_rte_species);

  m_electron_idx = 0;
  m_positive_idx = 1;
  m_photonZ_idx  = 0;

  m_ito_species[m_electron_idx] = RefCountedPtr<ito_species> (new electron());
  m_ito_species[m_positive_idx] = RefCountedPtr<ito_species> (new positive());
  m_rte_species[m_photonZ_idx]  = RefCountedPtr<rte_species> (new photonZ());

  // To avoid that MPI ranks draw the same particle positions, increment the seed for each rank
  m_seed += procID();
  m_rng   = std::mt19937_64(m_seed);

  List<ito_particle>& electrons = m_ito_species[m_electron_idx]->get_initial_particles();
  List<ito_particle>& positives = m_ito_species[m_positive_idx]->get_initial_particles();
  this->draw_gaussian_particles(electrons, positives, m_num_particles, m_blob_center, m_blob_radius, m_particle_weight);

  m_reactions.clear();
  m_reactions.emplace("impact_ionization", ito_reaction({m_electron_idx}, {m_electron_idx, m_electron_idx, m_positive_idx}));
  m_reactions.emplace("recombination",     ito_reaction({m_electron_idx, m_positive_idx}, {}));
  m_reactions.emplace("photo_excitation",  ito_reaction({m_electron_idx}, {m_electron_idx}, {m_photonZ_idx}));
}

ito_plasma_air2::~ito_plasma_air2(){

}

Real ito_plasma_air2::compute_alpha(const RealVect a_E) const {
  Real E = a_E.vectorLength();

  const Real alpha = (1.1944E6 + 4.3666E26/(E*E*E))*exp(-2.73E7/E);
  
  return alpha;
}

Real ito_plasma_air2::compute_eta(const RealVect a_E) const {
  return 340.75;
}

Vector<RealVect> ito_plasma_air2::compute_ito_velocities(const Real         a_time,
							 const RealVect     a_pos,
							 const RealVect     a_E,
							 const Vector<Real> a_cdr_densities) const {
  
  Vector<RealVect> velo(m_num_ito_species, RealVect::Zero);

  velo[m_electron_idx] = this->compute_electron_velocity(a_E);
  
  return velo;
}

RealVect ito_plasma_air2::compute_electron_velocity(const RealVect a_E) const {
  const Real E = a_E.vectorLength();
  const Real mu = 2.3987*pow(E, -0.26);
  
  return -mu*a_E;
}

Vector<Real> ito_plasma_air2::compute_ito_diffusion(const Real         a_time,
						    const RealVect     a_pos,
						    const RealVect     a_E,
						    const Vector<Real> a_cdr_densities) const {

  Vector<Real> D(m_num_ito_species, 0.0);

  D[m_electron_idx] = 4.3628E-3*pow(a_E.vectorLength(), 0.22);
  
  return D;
}

void ito_plasma_air2::advance_reaction_network(Vector<List<ito_particle>* >& a_particles,
					       Vector<List<photon>* >&       a_photons,
					       Vector<List<photon>* >&       a_newPhotons,
					       const RealVect                a_E,
					       const RealVect                a_cellPos,
					       const RealVect                a_centroid,
					       const RealVect                a_bndryCentroid,
					       const RealVect                a_bndryNormal,
					       const RealVect                a_lo,
					       const RealVect                a_hi,
					       const Real                    a_dx,
					       const Real                    a_kappa, 
					       const Real                    a_dt) const {

  this->advance_reaction_network_tau(a_particles, a_photons, a_newPhotons, a_E, a_cellPos, a_centroid,
				     a_bndryCentroid, a_bndryNormal, a_lo, a_hi, a_dx, a_kappa, a_dt);
}

void ito_plasma_air2::advance_reaction_network_tau(Vector<List<ito_particle>* >& a_particles,
						   Vector<List<photon>* >&       a_photons,
						   Vector<List<photon>* >&       a_newPhotons,
						   const RealVect                a_E,
						   const RealVect                a_pos,
						   const RealVect                a_centroid,
						   const RealVect                a_bndryCentroid,
						   const RealVect                a_bndryNormal,
						   const RealVect                a_lo,
						   const RealVect                a_hi,
						   const Real                    a_dx,
						   const Real                    a_kappa, 
						   const Real                    a_dt) const{



  Vector<int> photon_count(m_num_rte_species, 0);
  Vector<int> particle_count = this->get_particle_count(a_particles);

  const int Xe = particle_count[m_electron_idx];
  const int Xp = particle_count[m_positive_idx];

  // Compute the reaction rates. 
  const Real E       = a_E.vectorLength();
  const Real alpha   = this->compute_alpha(a_E);
  const Real eta     = this->compute_eta(a_E);
  const Real velo    = this->compute_electron_velocity(a_E).vectorLength();
  const Real xfactor = (m_pq/(m_p + m_pq))*excitation_rates(E)*sergey_factor(m_O2frac)*m_photoi_factor;

  m_reactions["impact_ionization"].rate() = alpha*velo;
  m_reactions["photo_excitation"].rate()  = alpha*velo*xfactor;
  m_reactions["recombination"].rate()     = eta*velo;

  // Do a tau-leaping step
  this->tau_leap(particle_count, photon_count, a_dt);

  // ==================================================
  // OLD CODE BELOW HERE
  // ==================================================

  // Compute number of electrons
  Real num_electrons = 0;
  for (ListIterator<ito_particle> lit(*a_particles[m_electron_idx]); lit.ok(); ++lit){
    num_electrons += lit().mass();
  }

  // Propensity functions
  const Real ionizationProp = alpha*velo*num_electrons;
  const Real recombProp     = eta*velo*num_electrons;
  const Real photoexcProp   = alpha*velo*num_electrons*xfactor;
  const Real totalProp      = ionizationProp + recombProp + photoexcProp;

  // Number of firings in each reaction channel
  const int num_ionizations = this->poisson_reaction(ionizationProp, a_dt);
  const int num_recomb      = this->poisson_reaction(recombProp, a_dt);
  const int num_photoexc    = this->poisson_reaction(photoexcProp, a_dt);

  List<ito_particle>& electrons = *a_particles[m_electron_idx];
  List<ito_particle>& posIons   = *a_particles[m_positive_idx];
  
  // Add electron-ion pairs
  this->add_particles(electrons, num_ionizations, a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
  this->add_particles(posIons,   num_ionizations, a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);

  // Photogeneration and photoionization
  this->add_photons(*a_newPhotons[m_photonZ_idx], *m_rte_species[m_photonZ_idx], num_photoexc, a_pos, a_lo, a_hi,
		    a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);

  // Photionization
  this->add_photoionization(*a_particles[m_electron_idx], *a_particles[m_positive_idx], *a_photons[m_photonZ_idx]);
}

Real ito_plasma_air2::excitation_rates(const Real a_E) const{
  const Real Etd = a_E/(m_N*units::s_Td);

  Real y = 1.0;
  if(Etd > 100){
    y = 0.1*exp(233/Etd);
  }

  return y;
}

Real ito_plasma_air2::sergey_factor(const Real a_O2frac) const{
  return 3E-2 + 0.4*pow(a_O2frac, 0.6);
}

Real ito_plasma_air2::photo_rate(const Real a_E) const {
  return excitation_rates(a_E)*sergey_factor(m_O2frac)*m_photoi_factor;
}

ito_plasma_air2::electron::electron(){
  m_mobile    = true;
  m_diffusive = true;
  m_name      = "electron";
  m_charge    = -1;
}

ito_plasma_air2::electron::~electron(){

}

ito_plasma_air2::positive::positive(){
  m_mobile    = false;
  m_diffusive = false;
  m_name      = "positive";
  m_charge    = 1;
}  

ito_plasma_air2::positive::~positive(){

}  

ito_plasma_air2::photonZ::photonZ(){
  m_name   = "photonZ";

  const Real O2_frac  = 0.2;
  const Real pressure = 1.0;
  
  ParmParse pp("ito_plasma_air2");
  
  pp.get("photoi_f1",   m_f1);
  pp.get("photoi_f2",   m_f2);
  pp.get("photoi_K1",   m_K1);
  pp.get("photoi_K2",   m_K2);
  pp.get("photoi_seed", m_seed);

  // Convert units
  m_pO2 = pressure*O2_frac*units::s_atm2pascal;
  m_K1  = m_K1*m_pO2;
  m_K2  = m_K2*m_pO2;

  // Seed the RNG
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  m_rng = new std::mt19937_64(m_seed);
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
}

ito_plasma_air2::photonZ::~photonZ(){

}

Real ito_plasma_air2::photonZ::get_kappa(const RealVect a_pos) const {
  const Real f = m_f1 + (*m_udist01)(*m_rng)*(m_f2 - m_f1);
  return m_K1*pow(m_K2/m_K1, (f-m_f1)/(m_f2-m_f1));
}
