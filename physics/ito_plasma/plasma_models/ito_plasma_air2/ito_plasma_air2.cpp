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
  pp.get   ("seed",           m_seed);
  pp.get   ("blob_radius",    m_blob_radius);
  pp.getarr("blob_center",    v, 0, SpaceDim); m_blob_center = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.get   ("num_particles",  m_num_particles);

  // Photostuff
  pp.get("quenching_pressure", m_pq);
  pp.get("photoi_factor",      m_photoi_factor);

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


  this->draw_initial_particles();
}

ito_plasma_air2::~ito_plasma_air2(){

}

void ito_plasma_air2::draw_initial_particles(){

  // To avoid that MPI ranks draw the same particle positions, increment the seed for each rank
  m_seed += procID();

  // Set up the RNG
  m_rng     = std::mt19937_64(m_seed);
  m_gauss   = std::normal_distribution<Real>(0.0, m_blob_radius);
  m_udist11 = std::uniform_real_distribution<Real>(-1., 1.);

  // Each MPI process draws the desired number of particles from a distribution
  const int quotient  = m_num_particles/numProc();
  const int remainder = m_num_particles % numProc();
  
  Vector<int> particlesPerRank(numProc(), quotient);
  
  for (int i = 0; i < remainder; i++){ 
    particlesPerRank[i] += 1;
  }

  List<ito_particle>& electrons = m_ito_species[m_electron_idx]->get_initial_particles();
  List<ito_particle>& positives = m_ito_species[m_positive_idx]->get_initial_particles();

  electrons.clear();
  positives.clear();

  // Now make the particles
  for (int i = 0; i < particlesPerRank[procID()]; i++){
    const Real weight  = 1.0;
    const RealVect pos = m_blob_center + random_gaussian();
    
    electrons.add(ito_particle(weight, pos));
    positives.add(ito_particle(weight, pos));
  }
}

RealVect ito_plasma_air2::random_gaussian(){

  const Real rad = m_gauss(m_rng);
  return rad*random_direction();
}

RealVect ito_plasma_air2::random_direction(){
  const Real EPS = 1.E-8;
#if CH_SPACEDIM==2
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = m_udist11(m_rng);
    x2 = m_udist11(m_rng);
    r  = x1*x1 + x2*x2;
  }

  return RealVect(x1,x2)/sqrt(r);
#elif CH_SPACEDIM==3
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = m_udist11(m_rng);
    x2 = m_udist11(m_rng);
    r  = x1*x1 + x2*x2;
  }

  const Real x = 2*x1*sqrt(1-r);
  const Real y = 2*x2*sqrt(1-r);
  const Real z = 1 - 2*r;

  return RealVect(x,y,z);
#endif
}

Real ito_plasma_air2::compute_alpha(const RealVect a_E) const {
  return 0.0;
}

Vector<RealVect> ito_plasma_air2::compute_ito_velocities(const Real         a_time,
							 const RealVect     a_pos,
							 const RealVect     a_E,
							 const Vector<Real> a_cdr_densities) const {
  
  Vector<RealVect> velo(m_num_ito_species, RealVect::Zero);
  
  const Real E = a_E.vectorLength();
  const Real mu = 2.3987*pow(E, -0.26);

  velo[m_electron_idx] = -mu*a_E;
  
  return velo;
}

Vector<Real> ito_plasma_air2::compute_ito_diffusion(const Real         a_time,
						    const RealVect     a_pos,
						    const RealVect     a_E,
						    const Vector<Real> a_cdr_densities) const{

  Vector<Real> D(m_num_ito_species, 0.0);

  D[m_electron_idx] = 4.3628E-3*pow(a_E.vectorLength(), 0.22);
  
  return D;
}

void ito_plasma_air2::advance_reaction_network(Vector<List<ito_particle>* >& a_particles,
					       Vector<List<photon>* >&       a_photons,
					       Vector<List<photon>* >&       a_newPhotons,
					       const RealVect                a_E,           
					       const RealVect                a_pos,
					       const Real                    a_dx,
					       const Real                    a_kappa, 
					       const Real                    a_dt) const {

  // Add a photon that propagates along -y. This is development code!
  List<photon>& srcPhotons = *a_newPhotons[m_photonZ_idx];
  const RealVect v = units::s_c0*RealVect(0.0, -1.0);
  srcPhotons.clear();
  srcPhotons.add(photon(a_pos, v, m_rte_species[m_photonZ_idx]->get_kappa(a_pos), 1.0));

  return;
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
