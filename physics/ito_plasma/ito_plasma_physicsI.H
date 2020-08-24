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

#include <memory>

#include <PolyGeom.H>

using namespace physics::ito_plasma;

inline RealVect ito_plasma_physics::random_position(const RealVect a_pos,
						    const RealVect a_lo,
						    const RealVect a_hi,
						    const RealVect a_bndryCentroid,
						    const RealVect a_bndryNormal,
						    const Real     a_dx,
						    const Real     a_kappa) const {

  RealVect pos;
  if(a_kappa < 1.0){ // Rejection sampling. 
    pos = this->random_position(a_lo, a_hi, a_bndryCentroid, a_bndryNormal);
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
						    const RealVect a_bndryNormal) const {

  RealVect pos = this->random_position(a_lo, a_hi);
  bool valid   = PolyGeom::dot(pos-a_bndryCentroid, a_bndryNormal) >= 0.0;

  while(!valid){
    pos    = this->random_position(a_lo, a_hi);
    valid = PolyGeom::dot(pos-a_bndryCentroid, a_bndryNormal) >= 0.0;
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

inline void ito_plasma_physics::reconcile_particles(Vector<List<ito_particle>* >& a_particles,
						    const Vector<int>&            a_newNumParticles,
						    const Vector<int>&            a_oldNumParticles,
						    const RealVect                a_pos,
						    const RealVect                a_lo,
						    const RealVect                a_hi,
						    const RealVect                a_bndryCentroid,
						    const RealVect                a_bndryNormal,
						    const Real                    a_dx,
						    const Real                    a_kappa) const{

  for (int i = 0; i < a_particles.size(); i++){
    this->reconcile_particles(*a_particles[i], a_newNumParticles[i], a_oldNumParticles[i], a_pos, a_lo, a_hi, a_bndryCentroid,
			      a_bndryNormal, a_dx, a_kappa);
  }
}

inline void ito_plasma_physics::reconcile_particles(List<ito_particle>& a_particles,
						    const int           a_newNumParticles,
						    const int           a_oldNumParticles,
						    const RealVect      a_pos,
						    const RealVect      a_lo,
						    const RealVect      a_hi,
						    const RealVect      a_bndryCentroid,
						    const RealVect      a_bndryNormal,
						    const Real          a_dx,
						    const Real          a_kappa) const {

  const int diff = a_newNumParticles - a_oldNumParticles;

  if(diff > 0){ // Need to add 'diff' mass.
    this->add_particles(a_particles, diff, a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
  }
  else if(diff < 0){ // Need to remove 'diff' mass.
    this->remove_particles(a_particles, -diff, a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
  }
}

inline void ito_plasma_physics::compute_particle_weights(int&      a_weight,
							 int&      a_num,
							 int&      a_remainder,
							 const int a_numNewParticles) const {
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
  //  MayDay::Warning("ito_plasma_physics::remove_particles - not implemented");
}

void ito_plasma_physics::reconcile_photons(Vector<List<photon>* >& a_newPhotons,
					   const Vector<int>&      a_numNewPhotons,
					   const RealVect          a_pos,
					   const RealVect          a_lo,
					   const RealVect          a_hi,
					   const RealVect          a_bndryCentroid,
					   const RealVect          a_bndryNormal,
					   const Real              a_dx,
					   const Real              a_kappa) const {

  for (int i = 0; i < a_newPhotons.size(); i++){
    if(a_numNewPhotons[i] > 0){
      this->reconcile_photons(*a_newPhotons[i], *m_rte_species[i], a_numNewPhotons[i], a_pos, a_lo, a_hi,
			      a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
    }
  }
}

inline void ito_plasma_physics::reconcile_photons(List<photon>&      a_photons,
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

inline void ito_plasma_physics::reconcile_photoionization(Vector<List<ito_particle>* >& a_particles,
							   const Vector<List<photon>* >& a_photons) const {
  for (const auto& r : m_photo_reactions){

    // Source and targets
    const int& src                  = r.second.source();
    const std::vector<int>& targets = r.second.targets();

    for (ListIterator<photon> lit(*a_photons[src]); lit.ok(); ++lit){
      const photon& phot  = lit();
      const RealVect pos  = phot.position();
      const Real mass     = phot.mass();

      for (const auto& t : targets){
	a_particles[t]->add(ito_particle(mass, pos));
      }
    }
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

inline Vector<int> ito_plasma_physics::get_particle_count(const Vector<List<ito_particle>* >& a_particles) const {
  Vector<int> ret(a_particles.size(), 0);

  for (int i = 0; i < a_particles.size(); i++){
    ret[i] = this->get_particle_count(*a_particles[i]);
  }

  return ret;
}

inline int ito_plasma_physics::get_particle_count(const List<ito_particle>& a_particles) const {
  int ret = 0;

  for (ListIterator<ito_particle> lit(a_particles); lit.ok(); ++lit){
    ret += round(lit().mass());
  }

  return ret;
}

inline void ito_plasma_physics::tau_leap(Vector<int>& a_particles, Vector<int>& a_newPhotons, const Real a_dt) const {

  const Vector<int> particles  = a_particles; // Because propensities should not change. 
  
  for (const auto& r : m_reactions){
    const int num_reactions = this->poisson_reaction(r.second.propensity(particles), a_dt);

    r.second.jump_state(a_particles, a_newPhotons, num_reactions);
  }
}

inline void ito_plasma_physics::ssa_step(Vector<int>& a_particles, Vector<int>& a_newPhotons, const Real a_dt) const {

  Real t   = 0.0;
  Real tau = 0.0;
  Real A   = 0.0;
  Real u1  = 0.0;
  Real u2  = 0.0;
  int  r   = 0;

  const int N = m_reactions.size();


  // A bit nasty, but we need a data structure with random access for evaluating stuff like
  // a_(r-1) <= u2*A < a_r. So, take the internal map and put it into a vector. 
  std::vector<std::shared_ptr<ito_reaction> > reactions;
  for (const auto& r : m_reactions){
    reactions.emplace_back(std::make_shared<ito_reaction>(r.second));
  }

  std::vector<Real> a(N, 0.0);

  while(t <= a_dt && N > 0){

    // Compute rates
    A = 0.0;
    for (int i = 0; i < N; i++){
      a[i] = reactions[i]->propensity(a_particles);
      A   += a[i];
    }

    u1 = m_udist01(m_rng);
    u2 = m_udist01(m_rng);

    // Time to next reaction
    tau = log(1./u1)/A;
    t  += a_dt;

    // Ok, time to react
    if(t <= a_dt){

      r = 0;
      for (int i = 1; i < N; i++){
	if(a[i-1] <= u2*A && u2*A < a[i]){
	  r = i;
	  break;
	}
      }

      reactions[r]->jump_state(a_particles, a_newPhotons, 1);
    }
  }
}

inline void ito_plasma_physics::advance_reaction_network(Vector<List<ito_particle>* >& a_particles,
							 Vector<List<photon>* >&       a_photons,
							 Vector<List<photon>* >&       a_newPhotons,
							 const RealVect                a_E,           
							 const RealVect                a_pos,
							 const RealVect                a_centroidPos,
							 const RealVect                a_bndryCentroid,
							 const RealVect                a_bndryNormal,
							 const RealVect                a_lo,
							 const RealVect                a_hi,
							 const Real                    a_dx,
							 const Real                    a_kappa, 
							 const Real                    a_dt){
  this->update_reaction_rates(a_E, a_dx, a_kappa);

  // Get particle counts. 
  Vector<int> newPhotonCount   = Vector<int>(m_num_rte_species, 0);
  Vector<int> oldParticleCount = this->get_particle_count(a_particles);
  Vector<int> newParticleCount = oldParticleCount;

  // Do a tau-leaping step
  this->tau_leap(newParticleCount, newPhotonCount, a_dt); 

  // Reconcile particles
  this->reconcile_particles(a_particles,
			    newParticleCount,
			    oldParticleCount,
			    a_pos,
			    a_lo,
			    a_hi,
			    a_bndryCentroid,
			    a_bndryNormal,
			    a_dx,
			    a_kappa);

  // Generate new photons
  this->reconcile_photons(a_newPhotons, newPhotonCount, a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);

  // Run photo-reactions
  this->reconcile_photoionization(a_particles, a_photons);
}

inline void ito_plasma_physics::hybrid_advance(Vector<int>& a_particles, Vector<int>& a_newPhotons, const Real a_dt){
  this->tau_leap(a_particles, a_newPhotons, a_dt);
}

#endif
