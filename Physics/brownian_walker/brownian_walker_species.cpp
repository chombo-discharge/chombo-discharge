/*!
  @file   brownian_walker_species.cpp
  @brief  Implementation of brownian_walker_species.H
  @author Robert Marskar
  @date   March 2020
*/

#include "brownian_walker_species.H"

#include <ParmParse.H>

#include "CD_NamespaceHeader.H"
using namespace physics::brownian_walker;

brownian_walker_species::brownian_walker_species() : m_rng(), m_gauss(0., 1.) {
  m_name   = "scalar species";
  m_chargeNumber = 0;

  ParmParse pp("brownian_walker");
  Vector<Real> v;

  pp.get   ("seed",           m_seed);
  pp.get   ("diffusion",      m_isDiffusive);
  pp.get   ("advection",      m_isMobile);
  pp.get   ("blob_amplitude", m_blob_amplitude);
  pp.get   ("blob_radius",    m_blob_radius);
  pp.getarr("blob_center",    v, 0, SpaceDim); m_blob_center = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.get   ("num_particles",  m_num_particles);


  draw_initial_particles();
}

brownian_walker_species::~brownian_walker_species(){

}

void brownian_walker_species::draw_initial_particles(){

  // To avoid that MPI ranks draw the same particle positions, increment the seed for each rank
  //  m_seed += procID();

  // Set up the RNG
  m_rng = std::mt19937_64(m_seed);
  m_gauss = std::normal_distribution<Real>(0.0, m_blob_radius);
  m_udist11 = std::uniform_real_distribution<Real>(-1., 1.);

  // Each MPI process draws the desired number of particles from a distribution
  const int quotient  = m_num_particles/numProc();
  const int remainder = m_num_particles % numProc();
  
  Vector<int> particlesPerRank(numProc(), quotient);
  
  for (int i = 0; i < remainder; i++){ 
    particlesPerRank[i] += 1;
  }

  // Now make the particles
  m_initialParticles.clear();
  for (int i = 0; i < particlesPerRank[procID()]; i++){
    const Real weight  = 1.0;
    const RealVect pos = m_blob_center + random_gaussian();
    m_initialParticles.add(ito_particle(weight, pos));
  }
}

RealVect brownian_walker_species::random_gaussian(){

  const Real rad = m_gauss(m_rng);
  return rad*randomDirection();
}

RealVect brownian_walker_species::randomDirection(){
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
#include "CD_NamespaceFooter.H"
