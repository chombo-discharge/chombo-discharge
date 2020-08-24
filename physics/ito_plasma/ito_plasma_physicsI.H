/*!
  @file   ito_plasma_physicsI.H 
  @brief  Implementation of ito_plasma_physics.H
  @author Robert Marskar
  @date   June 2020
*/

#ifndef _ITO_PLASMA_PHYSICSI_H_
#define _ITO_PLASMA_PHYSICSI_H_

#include "ito_plasma_physics.H"
#include "units.H"

#include <PolyGeom.H>

using namespace physics::ito_plasma;

inline RealVect ito_plasma_physics::random_position(const RealVect a_pos,
						    const RealVect a_lo,
						    const RealVect a_hi,
						    const RealVect a_bndryCentroid,
						    const RealVect a_normal,
						    const Real     a_dx,
						    const Real     a_kappa) const {

  RealVect pos;
  if(a_kappa < 1.0){ // Rejection sampling. 
    pos = this->random_position(a_lo, a_hi, a_bndryCentroid, a_normal);
  }
  else{ // Regular cell. Get a position. 
    pos = this->random_position(a_lo, a_hi);
  }

  pos = a_pos + pos*a_dx;

  return pos;
}

inline RealVect ito_plasma_physics::random_position(const RealVect a_lo,
						    const RealVect a_hi,
						    const RealVect a_bndryCentroid,
						    const RealVect a_normal) const {

  RealVect pos = this->random_position(a_lo, a_hi);
  bool valid   = PolyGeom::dot(pos-a_bndryCentroid, a_normal) >= 0.0;

  while(!valid){
    pos    = this->random_position(a_lo, a_hi);
    valid = PolyGeom::dot(pos-a_bndryCentroid, a_normal) > 0.0;
  }

  return pos;
}

inline RealVect ito_plasma_physics::random_position(const RealVect a_lo, const RealVect a_hi) const {

  RealVect pos = RealVect::Unit;

  for (int dir = 0; dir < SpaceDim; dir++){
    pos[dir] = a_lo[dir] + 0.5*(1.0 + m_udist11(m_rng))*(a_hi[dir] - a_lo[dir]);
  }

  return pos;
}

inline RealVect ito_plasma_physics::random_direction() const {
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

inline bool ito_plasma_physics::valid_advance(const Vector<List<ito_particle>* >& a_particles,
					      const Vector<int>&                  a_newParticles) const {
  bool ret = true;
  for (int i = 0; i < a_particles.size(); i++){
    if(a_particles[i]->length() + a_newParticles[i] < 0){
      ret = false;
      break;
    }
  }

  return ret;
}

inline void ito_plasma_physics::generate_particles(List<ito_particle>& a_particles,
						   const int           a_numNewParticles,
						   const RealVect      a_pos,
						   const RealVect      a_lo,
						   const RealVect      a_hi,
						   const RealVect      a_bndryCentroid,
						   const RealVect      a_bndryNormal,
						   const Real          a_dx,
						   const Real          a_kappa) const {

  if(a_numNewParticles > 0){ // We will add at most m_ppc particles
    this->add_particles(a_particles, a_numNewParticles, a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
  }
  else if(a_numNewParticles < 0){ // Need to remove mass
    this->remove_particles(a_particles, -a_numNewParticles, a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);

  }
}

inline void ito_plasma_physics::compute_particle_weights(int& a_weight, int& a_num, int& a_remainder, const int a_numNewParticles) const {
  if(a_numNewParticles <= m_ppc){  
    a_weight    = 1;
    a_num       = a_numNewParticles; 
    a_remainder = 0;
  }
  else{ // Add superparticles
    a_weight    = a_numNewParticles/m_ppc;
    a_remainder = a_numNewParticles%m_ppc;
    a_num       = (a_remainder == 0) ? m_ppc : m_ppc - 1;
  }
}

inline void ito_plasma_physics::add_particles(List<ito_particle>& a_particles,
					      const int           a_numNewParticles,
					      const RealVect      a_pos,
					      const RealVect      a_lo,
					      const RealVect      a_hi,
					      const RealVect      a_bndryCentroid,
					      const RealVect      a_bndryNormal,
					      const Real          a_dx,
					      const Real          a_kappa) const {
  int weight, num, remainder;
  this->compute_particle_weights(weight, num, remainder, a_numNewParticles);
  
  for (int i = 0; i < num; i++){
    const RealVect p = this->random_position(a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
    a_particles.add(ito_particle(weight, p));
  }

  if(remainder > 0){ // Rest of weight in case we got superparticles
    const RealVect p = this->random_position(a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
    a_particles.add(ito_particle(weight + remainder, p));
  }
}

inline void ito_plasma_physics::remove_particles(List<ito_particle>& a_particles,
						 const int           a_numParticlesToRemove,
						 const RealVect      a_pos,
						 const RealVect      a_lo,
						 const RealVect      a_hi,
						 const RealVect      a_bndryCentroid,
						 const RealVect      a_normal,
						 const Real          a_dx,
						 const Real          a_kappa) const {
  MayDay::Abort("ito_plasma_physics::remove_particles - not implemented");
}

inline void ito_plasma_physics::add_photons(List<photon>&      a_photons,
					    const rte_species& a_species,
					    const int          a_num_photons,
					    const RealVect     a_pos,
					    const RealVect     a_lo,
					    const RealVect     a_hi,
					    const RealVect     a_bndryCentroid,
					    const RealVect     a_bndryNormal,
					    const Real         a_dx,
					    const Real         a_kappa) const {
  a_photons.clear();

  int weight, num, remainder;
  this->compute_particle_weights(weight, num, remainder, a_num_photons);
  
  for (int i = 0; i < num; i++){
    const RealVect p = this->random_position(a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
    const RealVect v = units::s_c0*random_direction();
      
    a_photons.add(photon(a_pos, v, a_species.get_kappa(p), weight));
  }

  // If we used superphotons the last photon gets some extra oomph. 
  if(remainder > 0){
    const RealVect p = this->random_position(a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
    const RealVect v = units::s_c0*random_direction();
    a_photons.add(photon(a_pos, v, a_species.get_kappa(p), weight + remainder));
  }
}

inline void ito_plasma_physics::add_photoionization(List<ito_particle>& a_particles1,
						    List<ito_particle>& a_particles2,
						    const List<photon>& a_photons) const {
  
  for (ListIterator<photon> lit(a_photons); lit.ok(); ++lit){
    const photon& phot  = lit();
    const RealVect pos  = phot.position();
    const Real mass     = phot.mass();

    a_particles1.add(ito_particle(mass, pos));
    a_particles2.add(ito_particle(mass, pos));
  }
}

inline int ito_plasma_physics::poisson_reaction(const Real a_propensity, const Real a_dt) const{
  
  int value = 0;
  const Real mean = a_propensity*a_dt;

  if(mean < m_poisson_switch){
    std::poisson_distribution<int> dist(mean);
    value = dist(m_rng);
  }
  else{
    std::normal_distribution<double> dist(mean, sqrt(mean));
    value = dist(m_rng);
  }

  return Max(0,value);
}

inline void ito_plasma_physics::draw_gaussian_particles(List<ito_particle>& a_particles1,
							List<ito_particle>& a_particles2,
							const int           a_numParticles,
							const RealVect      a_blobCenter,
							const Real          a_blobRadius,
							const Real          a_weight){

  // Set up the Gaussian RNG
  std::normal_distribution<Real> gauss(0.0, a_blobRadius);


  // Each MPI process draws the desired number of particles from a distribution
  const int quotient  = a_numParticles/numProc();
  const int remainder = a_numParticles % numProc();
  
  Vector<int> particlesPerRank(numProc(), quotient);
  
  for (int i = 0; i < remainder; i++){ 
    particlesPerRank[i] += 1;
  }

  a_particles1.clear();
  a_particles2.clear();

  // Now make the particles
  for (int i = 0; i < particlesPerRank[procID()]; i++){
    const RealVect pos = a_blobCenter + gauss(m_rng)*random_direction();
    
    a_particles1.add(ito_particle(a_weight, pos));
    a_particles2.add(ito_particle(a_weight, pos));
  }
}

inline Vector<int> ito_plasma_physics::get_particle_count(const Vector<List<ito_particle>* >& a_particles){
  Vector<int> ret(a_particles.size());

  for (int i = 0; i < a_particles.size(); i++){
    ret[i] = a_particles[i]->length();
  }

  return ret;
}

#endif
