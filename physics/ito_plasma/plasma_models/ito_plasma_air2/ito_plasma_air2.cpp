/*!
  @file   ito_plasma_air2.cpp
  @brief  Implementation of ito_plasma_air2.H
  @author Robert Marskar
  @date   June 2020
*/

#include "ito_plasma_air2.H"

#include <ParmParse.H>

using namespace physics::ito_plasma;

ito_plasma_air2::ito_plasma_air2(){
  m_num_ito_species = 2;

  ParmParse pp("ito_plasma_air2");
  Vector<Real> v;
  
  // Get input parameters
  pp.get   ("seed",           m_seed);
  pp.get   ("blob_radius",    m_blob_radius);
  pp.getarr("blob_center",    v, 0, SpaceDim); m_blob_center = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.get   ("num_particles",  m_num_particles);

  // Set up species
  m_ito_species.resize(m_num_ito_species);
  m_rte_species.resize(0);

  m_electron_idx = 0;
  m_positive_idx = 1;

  m_ito_species[m_electron_idx] = RefCountedPtr<ito_species> (new electron());
  m_ito_species[m_positive_idx] = RefCountedPtr<ito_species> (new positive());

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

