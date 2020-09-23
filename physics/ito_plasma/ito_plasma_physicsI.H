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

#include <algorithm>
#include <unordered_set> 

#include <PolyGeom.H>

#define ITO_PLASMA_PHYSICS_DEBUG 0

using namespace physics::ito_plasma;

inline
RealVect ito_plasma_physics::random_position(const RealVect a_pos,
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

inline
RealVect ito_plasma_physics::random_position(const RealVect a_lo,
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

inline
RealVect ito_plasma_physics::random_position(const RealVect a_lo, const RealVect a_hi) const {

  RealVect pos = RealVect::Unit;

  for (int dir = 0; dir < SpaceDim; dir++){
    pos[dir] = a_lo[dir] + 0.5*(1.0 + m_udist11(m_rng))*(a_hi[dir] - a_lo[dir]);
  }

  return pos;
}

inline
RealVect ito_plasma_physics::random_direction() const {
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

inline
bool ito_plasma_physics::valid_state(const Vector<unsigned long long>& a_particles) const {
  bool ret = true;
  for (int i = 0; i < a_particles.size(); i++){
    if(a_particles[i] < 0){
      ret = false;
      break;
    }
  }

  return ret;
}

inline
void ito_plasma_physics::check_charge_conservation(const Vector<unsigned long long>& a_oldNumParticles,
						   const Vector<unsigned long long>& a_newNumParticles) const {
  constexpr long long zero = 0;
  
  long long charge_before = 0;
  long long charge_after = 0;

  for (int i = 0; i < a_oldNumParticles.size(); i++){
    const int q = m_ito_species[i]->get_charge();

    charge_before += a_oldNumParticles[i]*q;
    charge_after  += a_newNumParticles[i]*q;
  }

  const long long diff = charge_before - charge_after;

  if(diff != zero){
    MayDay::Abort("ito_plasma_physics::check_charge_conservation() - charge is not conserved!");
  }
}

inline
long long ito_plasma_physics::compute_charge(const Vector<List<ito_particle>* >& a_particles) const {

  long long ret = 0;
  for (int i = 0; i < a_particles.size(); i++){
    const long long q = (long long) m_ito_species[i]->get_charge();
    
    for (ListIterator<ito_particle> lit(*a_particles[i]); lit.ok(); ++lit){
      const long long w = llround(lit().mass());

      ret += w*q;
    }
  }

  return ret;
}

inline
void ito_plasma_physics::reconcile_particles(Vector<List<ito_particle>* >&     a_particles,
					     const Vector<unsigned long long>& a_newNumParticles,
					     const Vector<unsigned long long>& a_oldNumParticles,
					     const RealVect                    a_pos,
					     const RealVect                    a_lo,
					     const RealVect                    a_hi,
					     const RealVect                    a_bndryCentroid,
					     const RealVect                    a_bndryNormal,
					     const Real                        a_dx,
					     const Real                        a_kappa) const{

  for (int i = 0; i < a_particles.size(); i++){
    this->reconcile_particles(*a_particles[i], a_newNumParticles[i], a_oldNumParticles[i], a_pos, a_lo, a_hi, a_bndryCentroid,
			      a_bndryNormal, a_dx, a_kappa);
  }
}

inline
void ito_plasma_physics::reconcile_particles(List<ito_particle>&      a_particles,
					     const unsigned long long a_newNumParticles,
					     const unsigned long long a_oldNumParticles,
					     const RealVect           a_pos,
					     const RealVect           a_lo,
					     const RealVect           a_hi,
					     const RealVect           a_bndryCentroid,
					     const RealVect           a_bndryNormal,
					     const Real               a_dx,
					     const Real               a_kappa) const {

  const long long newNum = (long long) a_newNumParticles;
  const long long oldNum = (long long) a_oldNumParticles;
  
  const long long diff = newNum - oldNum;

  if(diff > 0){ // Need to add 'diff' mass.
    const unsigned long long Diff = (unsigned long long) (diff);
    this->add_particles(a_particles, Diff, a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
  }
  else if(diff < 0){ // Need to remove 'diff' mass.
    const unsigned long long Diff = (unsigned long long) (-diff);
    this->remove_particles(a_particles, Diff);
  }
}

inline
void ito_plasma_physics::compute_particle_weights(unsigned long long&      a_weight,
						  unsigned long long&      a_num,
						  unsigned long long&      a_remainder,
						  const unsigned long long a_numNewParticles) const {

  if(a_numNewParticles <= m_ppc){  
    a_weight    = 1;
    a_remainder = 0;
    a_num       = a_numNewParticles; 
  }
  else{ // Add superparticles
    a_weight    = a_numNewParticles/m_ppc;
    a_remainder = a_numNewParticles%m_ppc;
    a_num       = (a_remainder == 0) ? m_ppc : m_ppc - 1;

  }
}

inline
void ito_plasma_physics::add_particles(List<ito_particle>&      a_particles,
				       const unsigned long long a_numNewParticles,
				       const RealVect           a_pos,
				       const RealVect           a_lo,
				       const RealVect           a_hi,
				       const RealVect           a_bndryCentroid,
				       const RealVect           a_bndryNormal,
				       const Real               a_dx,
				       const Real               a_kappa) const {
  unsigned long long weight, num, remainder;

  constexpr unsigned long long zero = 0;
  
  this->compute_particle_weights(weight, num, remainder, a_numNewParticles);

  for (int i = 0; i < num; i++){
    const RealVect p = this->random_position(a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
    const Real w     = (Real) weight;
    a_particles.add(ito_particle(w, p));
  }

  if(remainder > zero){ // Rest of weight in case we got superparticles
    const RealVect p = this->random_position(a_pos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
    const Real w     = (Real) (weight + remainder);
    a_particles.add(ito_particle(w, p));
  }
}

inline
void ito_plasma_physics::remove_particles(List<ito_particle>&      a_particles,
					  const unsigned long long a_numParticlesToRemove) const {

#if ITO_PLASMA_PHYSICS_DEBUG // Debug
  this->sanity_check(a_particles);

  const unsigned long long weightBefore = this->get_total_weight(a_particles);
#endif

  constexpr unsigned long long zero = 0;
  
  ListIterator<ito_particle> lit(a_particles);
  Vector<unsigned long long> particleWeights(a_particles.length());
  
  unsigned long long minWeight  = std::numeric_limits<unsigned long long>::max();
  unsigned long long numRemoved = zero;

  // Compute weights and trim all particle weights down to the minimum weight. We shouldn't have weight = 0 when we
  // get here...
  this->compute_weights(lit, particleWeights, minWeight);                             // Computes weights
  this->trim_weights(numRemoved, particleWeights, a_numParticlesToRemove, minWeight); // Remove weight down to minWeight
  this->remove_weights(particleWeights, a_numParticlesToRemove - numRemoved);         // Uniformly removes rest of mass
  this->reconcile_weights(a_particles, particleWeights);                              // Reconcile particle weights

#if ITO_PLASMA_PHYSICS_DEBUG // Debug
  this->sanity_check(a_particles);
  const unsigned long long weightAfter = this->get_total_weight(a_particles);

  if(weightBefore - weightAfter != a_numParticlesToRemove){
    std::cout << "weight before = " << weightBefore << std::endl;
    std::cout << "weight after  = " << weightAfter << std::endl;
    std::cout << "Should remove = " << a_numParticlesToRemove << std::endl;
    std::cout << "numRemoved    = " << numRemoved << std::endl;
    
    MayDay::Abort("ito_plasma_physics::remove_particles - incorrect mass removed");
  }
#endif
}

inline
void ito_plasma_physics::set_mean_particle_energy(Vector<List<ito_particle>* >& a_particles, const Vector<Real> a_mean_energy) const{
  for (int i = 0; i < a_particles.size(); i++){
    this->set_mean_particle_energy(*a_particles[i], a_mean_energy[i]);
  }
}

inline
void ito_plasma_physics::set_mean_particle_energy(List<ito_particle>& a_particles, const Real a_mean_energy) const{
  for (ListIterator<ito_particle> lit(a_particles); lit.ok(); ++lit){
    lit().energy() = a_mean_energy;
  }
}

inline
void ito_plasma_physics::sanity_check(const List<ito_particle>& a_particles) const{

  ListIterator<ito_particle> lit(a_particles);
  
  for (lit.begin(); lit.ok(); ++lit){
    const long long w = llround(lit().mass());

    if(w < 1) MayDay::Abort("ito_plasma_physics::sanity_check - can't have w < 1!!!!");
  }
}

inline
unsigned long long ito_plasma_physics::get_total_weight(const List<ito_particle>& a_particles) const {

  unsigned long long ret = 0;
  for (ListIterator<ito_particle> lit(a_particles); lit.ok(); ++lit){
    const unsigned long long w = (unsigned long long) llround(lit().mass());

    ret += w;
  }

  return ret;
}

inline
void ito_plasma_physics::compute_weights(ListIterator<ito_particle>& a_lit,
					 Vector<unsigned long long>& a_particleWeights,
					 unsigned long long&         a_minWeight) const {
  int i = 0;
  for (a_lit.begin(); a_lit.ok(); ++a_lit, ++i){
    a_particleWeights[i] = (unsigned long long) llround(a_lit().mass());
    a_minWeight = std::min(a_particleWeights[i], a_minWeight);
  }
}

inline
void ito_plasma_physics::trim_weights(unsigned long long&         a_numRemoved,
				      Vector<unsigned long long>& a_particleWeights,
				      const unsigned long long    a_numParticlesToRemove,
				      const unsigned long long    a_minWeight) const{
  constexpr unsigned long long zero = 0;
  
  for (int i = 0; i < a_particleWeights.size(); i++){
    const unsigned long long diff1 = a_particleWeights[i]   - a_minWeight;  // Should always be > 0
    const unsigned long long diff2 = a_numParticlesToRemove - a_numRemoved; // Should always be > 0
    
    const unsigned long long r = std::max(zero, std::min(diff1, diff2));

    a_particleWeights[i] -= r;
    a_numRemoved         += r;
  }
}

inline
void ito_plasma_physics::remove_weights(Vector<unsigned long long>& a_particleWeights,
					const unsigned long long a_numParticlesToRemove) const {

  constexpr unsigned long long zero = 0;

  if(a_numParticlesToRemove > zero){
    const unsigned long long numCompParticles = (unsigned long long) a_particleWeights.size();
    const unsigned long long weight           = a_numParticlesToRemove/numCompParticles;
    const unsigned long long remainder        = a_numParticlesToRemove%numCompParticles; 

    // Remove weight from each particle. 
    if(weight > zero){
      for (int i = 0; i < a_particleWeights.size(); i++){
	a_particleWeights[i] -= weight;
      }
    }

    // Remove the remainder by cycling through the particles. 
    if(remainder > zero){
      unsigned long long W = zero;
      for (int i = 0; i < a_particleWeights.size(); i++){
	const unsigned long long w = std::min(a_particleWeights[i], remainder-W);
	
	a_particleWeights[i] -= w;
	W += w;
	
	if(W == remainder) { // OK, we have removed the remainder. Break out. 
	  break;
	}
      }
    }
  }
}

inline
void ito_plasma_physics::reconcile_weights(List<ito_particle>&               a_particles,
					   const Vector<unsigned long long>& a_particleWeights) const {

  constexpr unsigned long long zero = 0;  

  List<ito_particle> copyParticles(a_particles);
  a_particles.clear();

  int i = 0;
  for (ListIterator<ito_particle> lit(copyParticles); lit.ok(); ++lit, ++i){
    ito_particle& p = lit();

    if(a_particleWeights[i] > zero){
      p.mass() = (Real) a_particleWeights[i];
      a_particles.add(p);
    }
  }
}

void ito_plasma_physics::reconcile_photons(Vector<List<photon>* >&           a_newPhotons,
					   const Vector<unsigned long long>& a_numNewPhotons,
					   const RealVect                    a_pos,
					   const RealVect                    a_lo,
					   const RealVect                    a_hi,
					   const RealVect                    a_bndryCentroid,
					   const RealVect                    a_bndryNormal,
					   const Real                        a_dx,
					   const Real                        a_kappa) const {

  for (int i = 0; i < a_newPhotons.size(); i++){
    if(a_numNewPhotons[i] > 0){
      this->reconcile_photons(*a_newPhotons[i], *m_rte_species[i], a_numNewPhotons[i], a_pos, a_lo, a_hi,
			      a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
    }
  }
}

inline
void ito_plasma_physics::reconcile_photons(List<photon>&            a_photons,
					   const rte_species&       a_species,
					   const unsigned long long a_num_photons,
					   const RealVect           a_pos,
					   const RealVect           a_lo,
					   const RealVect           a_hi,
					   const RealVect           a_bndryCentroid,
					   const RealVect           a_bndryNormal,
					   const Real               a_dx,
					   const Real               a_kappa) const {
  a_photons.clear();

  unsigned long long weight, num, remainder;
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

inline
void ito_plasma_physics::reconcile_photoionization(Vector<List<ito_particle>* >& a_particles,
						   const Vector<List<photon>* >& a_photons) const {
  for (const auto& r : m_photo_reactions){

    // Source and targets
    const int& src                 = r.second.source();
    const std::deque<int>& targets = r.second.targets();

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

inline
unsigned long long ito_plasma_physics::poisson_reaction(const Real a_propensity, const Real a_dt) const{

  constexpr unsigned long long zero  = 0;
  unsigned long long value = 0;
  const Real mean = a_propensity*a_dt;

  if(mean < m_poisson_switch){
    std::poisson_distribution<unsigned long long> dist(mean);
    value = dist(m_rng);
  }
  else{
    std::normal_distribution<double> dist(mean, sqrt(mean));
    value = (unsigned long long) llround(dist(m_rng));
  }

  return Max(zero,value);
}

inline
int ito_plasma_physics::compute_particles_per_rank(const int a_numParticles) const {
  const int quotient  = a_numParticles/numProc();
  const int remainder = a_numParticles % numProc();
  
  Vector<int> particlesPerRank(numProc(), quotient);
  
  for (int i = 0; i < remainder; i++){ 
    particlesPerRank[i] += 1;
  }

  return particlesPerRank[procID()];
}

inline
void ito_plasma_physics::draw_gaussian_particles(List<ito_particle>& a_particles1,
						 List<ito_particle>& a_particles2,
						 const int           a_numParticles,
						 const RealVect      a_blobCenter,
						 const Real          a_blobRadius,
						 const Real          a_weight){

  // Set up the Gaussian RNG
  std::normal_distribution<Real> gauss(0.0, a_blobRadius);

  const int myNumParticles = this->compute_particles_per_rank(a_numParticles);

  // Now make the particles
  for (int i = 0; i < myNumParticles; i++){
    const RealVect pos = a_blobCenter + gauss(m_rng)*random_direction();
    
    a_particles1.add(ito_particle(a_weight, pos));
    a_particles2.add(ito_particle(a_weight, pos));
  }
}

inline
void ito_plasma_physics::draw_sphere_particles(List<ito_particle>& a_particles1,
					       List<ito_particle>& a_particles2,
					       const int           a_numParticles,
					       const RealVect      a_blobCenter,
					       const Real          a_blobRadius,
					       const Real          a_weight){

  const int myNumParticles = this->compute_particles_per_rank(a_numParticles);

  // Now make the particles
  for (int i = 0; i < myNumParticles; i++){
    bool accept  = false;
    RealVect pos = RealVect::Zero;
    
    while(!accept){
      for (int d = 0; d < SpaceDim; d++){
	pos[d] = m_udist11(m_rng);
      }
      
      if(pos.vectorLength() <= 1.0){
	accept = true;
	pos = pos*a_blobRadius + a_blobCenter;
      }
    }
    
    // Add the particle
    a_particles1.add(ito_particle(a_weight, pos));
    a_particles2.add(ito_particle(a_weight, pos));
  }
}

inline
void ito_plasma_physics::draw_box_particles(List<ito_particle>& a_particles1,
					    List<ito_particle>& a_particles2,
					    const int           a_numParticles,
					    const RealVect      a_boxLo,
					    const RealVect      a_boxHi,
					    const Real          a_weight) {

  const int myNumParticles = this->compute_particles_per_rank(a_numParticles);

  const RealVect diff = a_boxHi - a_boxLo;
  for (int i = 0; i < myNumParticles; i++){
    const RealVect pos = a_boxLo + m_udist01(m_rng)*diff;

    a_particles1.add(ito_particle(a_weight, pos));
    a_particles2.add(ito_particle(a_weight, pos));
  }
}

inline
Vector<unsigned long long> ito_plasma_physics::get_particle_count(const Vector<List<ito_particle>* >& a_particles) const {
  Vector<unsigned long long> ret(a_particles.size(), 0);

  for (int i = 0; i < a_particles.size(); i++){
    ret[i] = this->get_particle_count(*a_particles[i]);
  }

  return ret;
}

inline
unsigned long long ito_plasma_physics::get_particle_count(const List<ito_particle>& a_particles) const {
  unsigned long long ret = 0;

  for (ListIterator<ito_particle> lit(a_particles); lit.ok(); ++lit){
    ret += (unsigned long long) llround(lit().mass());
  }

  return ret;
}

inline
Vector<Real> ito_plasma_physics::get_particle_energy(const Vector<List<ito_particle>* >& a_particles) const {
  Vector<Real> ret(a_particles.size(), 0.0);

  for (int i = 0; i < a_particles.size(); i++){
    ret[i] = this->get_particle_energy(*a_particles[i]);
  }

  return ret;
}

inline
Real ito_plasma_physics::get_particle_energy(const List<ito_particle>& a_particles) const {
  Real ret = 0;

  for (ListIterator<ito_particle> lit(a_particles); lit.ok(); ++lit){
    ret += lit().energy();
  }

  return ret;
}

inline
Vector<Real> ito_plasma_physics::get_particle_mean_energy(const Vector<Real>& a_energies, const Vector<unsigned long long>& a_particle_count) const {
  Vector<Real> ret(a_energies.size(), 0.0);

  for (int i = 0; i < a_energies.size(); i++){
    ret[i] = this->get_particle_mean_energy(a_energies[i], a_particle_count[i]);
  }

  return ret;
}

inline
Real ito_plasma_physics::get_particle_mean_energy(const Real& a_energy, const unsigned long long& a_particle_count) const{
  const Real ret = (a_particle_count > 0) ? a_energy/a_particle_count : 0.0;
  return ret;
}

inline
void ito_plasma_physics::advance_reaction_network_lfa(Vector<List<ito_particle>* >& a_particles,
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
  
  // Get particle counts. 
  Vector<unsigned long long> newPhotonCount   = Vector<unsigned long long>(m_num_rte_species, 0);
  Vector<unsigned long long> oldParticleCount = this->get_particle_count(a_particles);
  Vector<unsigned long long> newParticleCount = oldParticleCount;
  Vector<Real> mean_energies(m_num_ito_species, 0.0);

#if ITO_PLASMA_PHYSICS_DEBUG
  const long long Qbefore = this->compute_charge(a_particles);
#endif

  // Number advance
  this->advance_particles(newParticleCount, newPhotonCount, a_dt);

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

  // Run photo-reactions and increment the number of particles as needed
  this->reconcile_photoionization(a_particles, a_photons);

  // Set new mean energies
  this->set_mean_particle_energy(a_particles, mean_energies);


#if ITO_PLASMA_PHYSICS_DEBUG // Debug code
  this->check_charge_conservation(oldParticleCount, newParticleCount);

  const long long Qafter = this->compute_charge(a_particles);

  constexpr long long zero = 0;
  if(Qbefore - Qafter != zero){
    std::cout << "Qbefore = " << Qbefore << std::endl;
    std::cout << "Qafter  = " << Qafter << std::endl;
    MayDay::Abort("ito_plasma_physics::advance_reaction_network_lfa(...) - charge is not conserved!!!");
  }
#endif
}

void ito_plasma_physics::advance_reaction_network_lea(Vector<List<ito_particle>* >& a_particles,
						      Vector<List<photon>* >&       a_photons,
						      Vector<List<photon>* >&       a_newPhotons,
						      const Vector<Real>&           a_sources, 
						      const RealVect                a_E,           
						      const RealVect                a_cellPos,
						      const RealVect                a_centroidPos,
						      const RealVect                a_bndryCentroid,
						      const RealVect                a_bndryNormal,
						      const RealVect                a_lo,
						      const RealVect                a_hi,
						      const Real                    a_dx,
						      const Real                    a_kappa, 
						      const Real                    a_dt){

  // Get particle counts. 
  Vector<unsigned long long> newPhotonCount   = Vector<unsigned long long>(m_num_rte_species, 0);
  Vector<unsigned long long> oldParticleCount = this->get_particle_count(a_particles);
  Vector<Real>               energies         = this->get_particle_energy(a_particles);
  Vector<Real>               mean_energies    = this->get_particle_mean_energy(energies, oldParticleCount);
  Vector<unsigned long long> newParticleCount = oldParticleCount;
  

  // Doa  number advance in the LEA. 
  this->advance_particles(newParticleCount, newPhotonCount, mean_energies, a_sources, a_dt, a_E, a_dx, a_kappa);
  
  // Reconcile particles
  this->reconcile_particles(a_particles,
  			    newParticleCount,
  			    oldParticleCount,
  			    a_cellPos,
  			    a_lo,
  			    a_hi,
			    a_bndryCentroid,
  			    a_bndryNormal,
  			    a_dx,
  			    a_kappa);

  // Generate new photons
  this->reconcile_photons(a_newPhotons, newPhotonCount, a_cellPos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);

  // Run photo-reactions and increment the number of particles as needed
  this->reconcile_photoionization(a_particles, a_photons);

  // Set new mean energies for all the particles. 
  this->set_mean_particle_energy(a_particles, mean_energies);
}

inline
ReactionTypes ito_plasma_physics::get_critical_reactions(const Vector<unsigned long long>& a_particles) const{

  ReactionPtrs criticalReactions;    // std::deque<shared_ptr<ito_reaction> > under the hood
  ReactionPtrs nonCriticalReactions; // std::deque<shared_ptr<ito_reaction> > under the hood

  const long long Ncrit = (long long) m_Ncrit; // Need long long version of this because because a_particles is long long
  
  for (const auto& r : m_reactions){

    const std::map<int, int>& stateChange = r.second.get_state_change();

    long long Lj    = std::numeric_limits<long long>::max();
    bool consuming  = false;

    // Check how many times reaction r can fire before consuming all its species
    for (const auto& s: stateChange){

      const int reactant = s.first;
      const int nuIJ     = s.second;

      if(nuIJ < 0){
	consuming = true;
	
	const long long s = (long long) std::abs(nuIJ);
	const long long b = (long long) (a_particles[reactant] + s-1)/s;
	
	Lj = std::min(Lj, b);
      }
    }

    // Check if reaction is a critical reaction
    bool crit = false;
    if(Lj < Ncrit && consuming){
      crit = true;
    }

    // Move reactions to where they should be
    if(crit){
      criticalReactions.emplace_back(std::make_shared<ito_reaction>(r.second));
    }
    else{
      nonCriticalReactions.emplace_back(std::make_shared<ito_reaction>(r.second));
    }
  }

  // Return the critical and non-critical reactions
  ReactionTypes ret = ReactionTypes();
  ret.first  = criticalReactions;
  ret.second = nonCriticalReactions;

  return ret;
}

inline
Real ito_plasma_physics::get_critical_step(const Vector<unsigned long long>& a_particles,
					   const ReactionPtrs&               a_critical_reactions) const {
						 
  // TLDR: a_critical_reaction is a vector under the hood.
  Real dt = std::numeric_limits<Real>::max();
  
  if(a_critical_reactions.size() > 0){
    Real A = 0.0;
    for (const auto& r : a_critical_reactions){
      A += r->propensity(a_particles);
    }
    
    const Real u1 = m_udist01(m_rng);
    dt = log(1./u1)/A;
  }

  return dt;
}

inline
Real ito_plasma_physics::get_non_critical_step(const Vector<unsigned long long>& a_particles,
					       const ReactionPtrs&               a_non_critical_reactions) const {
						 
  // TLDR: a_non_critical_reaction is a vector under the hood.
  
  Real dt = std::numeric_limits<Real>::max();

  if(a_non_critical_reactions.size() > 0) {


    std::deque<Real> propensities;
    for (const auto& r : a_non_critical_reactions){
      const Real p = r->propensity(a_particles);

      propensities.emplace_back(p);
    }

    // 1. Compute propensities of all reactions. It is possible for propensities to be zero. Order matters,
    //    so use a deque
    std::deque<std::pair<Real, std::shared_ptr<ito_reaction> > > reactions;
    for(const auto& r : a_non_critical_reactions){
      const Real p = r->propensity(a_particles);

      reactions.emplace_back(std::pair<Real, std::shared_ptr<ito_reaction> >(p, r));
    }
      

    // 2. Get the list of all reactants
    std::unordered_set<int> allReactants;
    for (const auto& reaction : a_non_critical_reactions){

      const auto& cur_reactants = reaction->get_reactants();
      for (const auto& l : cur_reactants){
	allReactants.insert(l); 
      }
    }


    // 3. Iterate through all reactants and compute deviations
    for (const auto& reactant : allReactants){

      Real mu    = 0.0;
      Real sigma = 0.0;

      constexpr int gi = 1; // Set gi to 1 for now, should think about this later. 

      for (const auto& R : reactions){
	const Real& ap        = R.first;     // Reaction propensity.
	const ito_reaction& r = *(R.second); // Reaction. 
	const int muIJ        = r.get_state_change(reactant); // State change of reactant due to current reaction

	mu    += muIJ * ap;
	sigma += muIJ * muIJ * ap;

      }

      Real dt1 = std::numeric_limits<Real>::max();
      Real dt2 = std::numeric_limits<Real>::max();

      if(mu != 0.0) {
	dt1 = std::max(m_eps*a_particles[reactant]/gi, 1.)/std::abs(mu);
      }
      if(sigma != 0.0) {
	const Real f = std::max(m_eps*a_particles[reactant]/gi,1.0);
	dt2 = f*f/std::abs(sigma);
      }

      dt = std::min(dt, std::min(dt1, dt2));
    }

  }


  return dt;
}

inline
Real ito_plasma_physics::propensity(const Vector<unsigned long long>& a_particles) const{

  Real A = 0.0;
  for (const auto& r : m_reactions){
    A += r.second.propensity(a_particles);
  }

  return A;
}


inline
Real ito_plasma_physics::propensity(const Vector<unsigned long long>& a_particles, const ReactionPtrs& a_critical_reactions) const{

  Real A = 0.0;
  std::deque<Real> propensities = this->propensities(a_particles, a_critical_reactions);
  
  for (const auto& p : propensities){
    A += p;
  }

  return A;
}

inline
std::deque<Real> ito_plasma_physics::propensities(const Vector<unsigned long long>& a_particles) const {
  std::deque<Real> propensities;
  
  for (const auto& r : m_reactions){
    const Real p = r.second.propensity(a_particles);
    propensities.emplace_back(p);
  }

  return propensities;
}

inline
std::deque<Real> ito_plasma_physics::propensities(const Vector<unsigned long long>& a_particles,
						  const ReactionPtrs&               a_critical_reactions) const{
  std::deque<Real> propensities;

  for (const auto& r : a_critical_reactions){
    const Real p = r->propensity(a_particles);
    propensities.emplace_back(p);
  }

  return propensities;
}
inline
void ito_plasma_physics::ssa_step(Vector<unsigned long long>& a_particles,
				  Vector<unsigned long long>& a_newPhotons,
				  const std::deque<Real>&     a_propensities,
				  const ReactionPtrs&         a_reactions,
				  const Real                  a_A) const {
  const Real u2 = m_udist01(m_rng);
  const int M   = a_reactions.size();
  
  constexpr unsigned long long one = 1;
  
  if(a_reactions.size() > 0){
    int r = 0;

    for (int i = 1; i < M; i++){
      if(a_propensities[i-1] <= u2*a_A && u2*a_A < a_propensities[i]){
	r = i;
	break;
      }
    }

    a_reactions[r]->jump_state(a_particles, a_newPhotons, one);
  }
}

inline
void ito_plasma_physics::ssa_step(Vector<unsigned long long>& a_particles,
				  Vector<unsigned long long>& a_newPhotons,
				  const std::deque<Real>&     a_propensities,
				  const Real                  a_totalPropensity) const {
  
  // Map the map onto a deque. 
  std::deque<std::shared_ptr<ito_reaction> > reactions;
  for (const auto& r : m_reactions){
    reactions.emplace_back(std::make_shared<ito_reaction>(r.second));
  }

  this->ssa_step(a_particles, a_newPhotons, a_propensities, reactions, a_totalPropensity);
}

inline
void ito_plasma_physics::tau_step(Vector<unsigned long long>& a_particles,
				  Vector<unsigned long long>& a_newPhotons,
				  const std::deque<Real>&     a_propensities,
				  const ReactionPtrs&         a_reactions,
				  const Real                  a_dt) const {
  
  const Vector<unsigned long long> particles = a_particles;
  for (const auto& r : a_reactions){
    const unsigned long long num_reactions = this->poisson_reaction(r->propensity(particles), a_dt);
    
    r->jump_state(a_particles, a_newPhotons, num_reactions);
  }
}

inline
void ito_plasma_physics::tau_step(Vector<unsigned long long>& a_particles,
				  Vector<unsigned long long>& a_newPhotons,
				  const std::deque<Real>&     a_propensities,
				  const Real                  a_dt) const {

  // Map the map onto a deque
  std::deque<std::shared_ptr<ito_reaction> > reactions;
  for (const auto& r : m_reactions){
    reactions.emplace_back(std::make_shared<ito_reaction>(r.second));
  }

  this->tau_step(a_particles, a_newPhotons, a_propensities, reactions, a_dt);
}

inline
void ito_plasma_physics::advance_ssa(Vector<unsigned long long>& a_particles,
				     Vector<unsigned long long>& a_newPhotons,
				     const Real                  a_dt) const {
  constexpr unsigned long long one = 1;
  if(m_reactions.size() > 0){
    Real curDt = 0.0;

    std::deque<std::shared_ptr<ito_reaction> > reactions;
    for (const auto& r : m_reactions){
      reactions.emplace_back(std::make_shared<ito_reaction>(r.second));
    }

    const int M  = reactions.size();

    Real u1;
    Real u2;
    Real A;
    Real dt;

    int r;

    while(curDt <= a_dt){

      A  = this->propensity(a_particles);
      u1 = m_udist01(m_rng);
      dt = log(1./u1)/A;

      curDt += dt;

      if(curDt <= a_dt){

	const std::deque<Real> propensities = this->propensities(a_particles);
    
	u2 = m_udist01(m_rng);

	r = 0;
	for (int i = 1; i < M; i++){
	  if(propensities[i-1] <= u2*A && u2*A < propensities[i]){
	    r = i;
	    break;
	  }
	}

	reactions[r]->jump_state(a_particles, a_newPhotons, one);
      }
    }
  }
}

inline
void ito_plasma_physics::advance_tau(Vector<unsigned long long>& a_particles,
				     Vector<unsigned long long>& a_newPhotons,
				     const Real                  a_dt) const {
  std::deque<Real> propensities = this->propensities(a_particles);
  this->tau_step(a_particles, a_newPhotons, propensities, a_dt);
}

inline
void ito_plasma_physics::advance_hybrid(Vector<unsigned long long>& a_particles,
					Vector<unsigned long long>& a_newPhotons,
					const Real                  a_dt) const {

  Real curTime = 0.0;

  while(curTime < a_dt){

    // Get critical and non-critical reactions. 
    ReactionTypes rTypes = this->get_critical_reactions(a_particles);

    ReactionPtrs criticalReactions    = rTypes.first;
    ReactionPtrs nonCriticalReactions = rTypes.second;

    // Compute time steps for critical and non-critical reactions
    Real dtCrit    = this->get_critical_step(a_particles,     criticalReactions);
    Real dtNonCrit = this->get_non_critical_step(a_particles, nonCriticalReactions);

    bool valid = false;
    while(!valid){ 
      
      const Real curDt   = std::min(a_dt-curTime, std::min(dtCrit, dtNonCrit));
      const bool nonCrit = dtNonCrit < dtCrit || criticalReactions.size() == 0 || dtCrit < (a_dt-curTime);
	    
      // How will we advance?
      if(nonCrit){// No critical reactions fire. Advance only noncritical reactions.

	// Compute the total propensity function in case it is faster to run SSA steps
	Real A = this->propensity(a_particles);
	
	if(A*curDt < m_SSAlim){ // Tau-leaping is inefficient, do some SSA on the whole reaction set. 
	  for (int i = 0; i < m_NSSA; i++){

	    const std::deque<Real> propensities = this->propensities(a_particles);

	    // Draw time to next reaction. 
	    const Real u1     = m_udist01(m_rng);
	    const Real dtNext = log(1./u1)/A;

	    valid    = true;
	    curTime += dtNext;
	    if(curTime <= a_dt){ // Ok, we are inside [0, a_dt] and we can do the reaction. 
	      this->ssa_step(a_particles, a_newPhotons, propensities, A);

	      A = this->propensity(a_particles); // Recompute propensities
	    }
	    else{
	      break; // No need to perform the rest. 
	    }
	  }
	}
	else{ // Perform tau-leaping with no-critical reactions only
	  // Need a copy in case step is rejected
	  Vector<unsigned long long> particles = a_particles;
	  Vector<unsigned long long> photons   = a_newPhotons;
      
	  const std::deque<Real> propensitiesNC = this->propensities(particles, nonCriticalReactions);
	
	  this->tau_step(particles, photons, propensitiesNC, nonCriticalReactions, curDt);

	  // Make sure step was valid
	  valid = this->valid_state(particles);
	  if(valid){
	    curTime     += curDt;
	    a_particles  = particles;
	    a_newPhotons = photons;
	  }
	  else{
	    dtNonCrit = dtNonCrit*0.5;
	    continue;
	  }
	}
      } 
      else{ // One critical reaction fired. Advance that reaction and use tau-leaping for the rest.
	// Need a copy in case step is rejected
	Vector<unsigned long long> particles = a_particles;
	Vector<unsigned long long> photons   = a_newPhotons;

	const std::deque<Real> propensitiesC  = this->propensities(particles, criticalReactions);
	const std::deque<Real> propensitiesNC = this->propensities(particles, nonCriticalReactions);
	
	const Real Ac = this->propensity(particles, criticalReactions);

	// SSA for critical reactions
	this->ssa_step(particles, photons, propensitiesC, criticalReactions, Ac);

	// Tau-leaping for non-critical reactions
	this->tau_step(particles, photons, propensitiesNC, nonCriticalReactions, curDt);

	valid = this->valid_state(particles);

	if(valid){
	  curTime     += curDt;
	  a_particles  = particles;
	  a_newPhotons = photons;
	}
	else{
	  dtNonCrit = dtNonCrit * 0.5;
	  continue;
	}
      }
    }
  }
}

inline
void ito_plasma_physics::advance_particles(Vector<unsigned long long>& a_particles,
					   Vector<unsigned long long>& a_newPhotons,
					   const Real a_dt) const {

  // Which algorithm?
  if(m_algorithm == algorithm::ssa){
    this->advance_ssa(a_particles, a_newPhotons, a_dt); 
  }
  else if(m_algorithm == algorithm::tau){
    this->advance_tau(a_particles, a_newPhotons, a_dt); 
  }
  else if(m_algorithm == algorithm::hybrid){
    this->advance_hybrid(a_particles, a_newPhotons, a_dt); 
  }
}

inline
void ito_plasma_physics::update_reaction_rates(const RealVect a_E, const Vector<Real> a_mean_energies, const Real a_dx, const Real a_kappa) const {
  
  // Update the rate based on the mean energies or the electric field, depending on how we use the LFA/LEA. 
  if(m_coupling == coupling::LFA){
    this->update_reaction_rates_lfa(a_E, a_dx, a_kappa);
  }
  else if(m_coupling == coupling::LEA){
    this->update_reaction_rates_lea(a_E, a_mean_energies, a_dx, a_kappa);
  }
  
}

inline
void ito_plasma_physics::tau_step(Vector<unsigned long long>& a_particles,
				  Vector<unsigned long long>& a_newPhotons,
				  Vector<Real>&               a_mean_energies,
				  const Vector<Real>&         a_sources,
				  const ReactionPtrs&         a_reactions,
				  const Real                  a_dx,
				  const Real                  a_kappa,
				  const Real                  a_dt) const {

  // mean_energies -> total energies
  Vector<Real> energies(a_mean_energies);
  for (int i = 0; i < a_particles.size(); i++){
    energies[i] *= a_particles[i];
  }

  // Do a tau step
  for (const auto& r : a_reactions){
    const unsigned long long num_reactions = this->poisson_reaction(r->propensity(a_particles), a_dt);

    // Jump state and energy
    r->jump_state(a_particles, a_newPhotons, num_reactions);    // 
    r->jump_energy(energies,   a_mean_energies, num_reactions); // Not necessary for LFA, but not expensive either....
  }

  // Recompute mean energies afterwards.
  const Real dV = a_kappa*pow(a_dx, SpaceDim);
  for (int i = 0; i < a_particles.size(); i++){
    if(a_particles[i] > 0){
      a_mean_energies[i] = Max(0.0, (energies[i] + a_sources[i]*dV*a_dt)/(a_particles[i]));
    }
    else{
      a_mean_energies[i] = 0.0;
    }
  }
}

inline
void ito_plasma_physics::tau_step(Vector<unsigned long long>& a_particles,
				  Vector<unsigned long long>& a_newPhotons,
				  Vector<Real>&               a_mean_energies,
				  const Vector<Real>&         a_sources,
				  const Real                  a_dx,
				  const Real                  a_kappa,
				  const Real                  a_dt) const {

  std::deque<std::shared_ptr<ito_reaction> > reactions;
  for (const auto& r : m_reactions){
    reactions.emplace_back(std::make_shared<ito_reaction>(r.second));
  }

  this->tau_step(a_particles, a_newPhotons, a_mean_energies, a_sources, a_dx, a_kappa, a_dt);
}

inline
void ito_plasma_physics::advance_tau(Vector<unsigned long long>& a_particles,
				     Vector<unsigned long long>& a_newPhotons,
				     Vector<Real>&               a_mean_energies,
				     const Vector<Real>&         a_sources,
				     const Real                  a_dt,
				     const RealVect              a_E,
				     const Real                  a_dx,
				     const Real                  a_kappa) const {

  // Update reaction rates in case we use LEA. 
  this->update_reaction_rates(a_E, a_mean_energies, a_dx, a_kappa);

  // Do a tau step
  this->tau_step(a_particles, a_newPhotons, a_mean_energies, a_sources, a_dt, a_dx, a_kappa);
}


inline
void ito_plasma_physics::ssa_step(Vector<unsigned long long>& a_particles,
				  Vector<unsigned long long>& a_newPhotons,
				  Vector<Real>&               a_mean_energies,
				  const Vector<Real>&         a_sources,
				  const Real                  a_dx,
				  const Real                  a_kappa,
				  const Real                  a_dt) const {

  
  std::deque<std::shared_ptr<ito_reaction> > reactions;
  for (const auto& r : m_reactions){
    reactions.emplace_back(std::make_shared<ito_reaction>(r.second));
  }

  // Compute all propensities
  Real A = 0.0;
  std::deque<Real> propensities;
  for (const auto& r : reactions){
    const Real a = r->propensity(a_particles);
    propensities.emplace_back(a);
    A += a;
  }

  this->ssa_step(a_particles, a_newPhotons, a_mean_energies, a_sources, propensities, reactions, A, a_dx, a_kappa, a_dt);
}

inline
void ito_plasma_physics::ssa_step(Vector<unsigned long long>& a_particles,
				  Vector<unsigned long long>& a_newPhotons,
				  Vector<Real>&               a_mean_energies,
				  const Vector<Real>&         a_sources,
				  const std::deque<Real>&     a_propensities,
				  const Real                  a_A,
				  const Real                  a_dx,
				  const Real                  a_kappa,
				  const Real                  a_dt) const {

  std::deque<std::shared_ptr<ito_reaction> > reactions;
  for (const auto& r : m_reactions){
    reactions.emplace_back(std::make_shared<ito_reaction>(r.second));
  }

  this->ssa_step(a_particles, a_newPhotons, a_mean_energies, a_sources, a_propensities, reactions, a_A, a_dx, a_kappa, a_dt);
}

inline
void ito_plasma_physics::ssa_step(Vector<unsigned long long>& a_particles,
				  Vector<unsigned long long>& a_newPhotons,
				  Vector<Real>&               a_mean_energies,
				  const Vector<Real>&         a_sources,
				  const std::deque<Real>&     a_propensities,
				  const ReactionPtrs&         a_reactions,
				  const Real                  a_A,
				  const Real                  a_dx,
				  const Real                  a_kappa,
				  const Real                  a_dt) const {

  const int M = a_reactions.size();
  
  constexpr unsigned long long one = 1;
  
  if(M > 0){
    
    // Mean_energies -> total energies
    Vector<Real> energies(a_mean_energies);
    for (int i = 0; i < a_particles.size(); i++){
      energies[i] *= a_particles[i];
    }


    // RNG. 
    const Real u2 = m_udist01(m_rng);

    // Determine the reaction type. 
    int r = 0; 
    for (int i = 1; i < M; i++){
      if(a_propensities[i-1] <= u2*a_A && u2*a_A < a_propensities[i]){
	r = i;
	break;
      }
    }

    // Jump state and energy
    a_reactions[r]->jump_state(a_particles, a_newPhotons,    one);
    a_reactions[r]->jump_energy(energies,   a_mean_energies, one); 


    // Recompute mean energies afterwards.
    const Real dV = a_kappa*pow(a_dx, SpaceDim);
    for (int i = 0; i < a_particles.size(); i++){
      if(a_particles[i] > 0){
	a_mean_energies[i] = Max(0.0, (energies[i] + a_sources[i]*dV*a_dt)/(a_particles[i]));
      }
      else{
	a_mean_energies[i] = 0.0;
      }
    }
  }
}



inline
void ito_plasma_physics::advance_ssa(Vector<unsigned long long>& a_particles,
				     Vector<unsigned long long>& a_newPhotons,
				     Vector<Real>&               a_mean_energies,
				     const Vector<Real>&         a_sources,
				     const Real                  a_dt,
				     const RealVect              a_E,
				     const Real                  a_dx,
				     const Real                  a_kappa) const {

  constexpr unsigned long long one = 1;

  const int M  = m_reactions.size();
  
  if(M > 0){
    Real curDt = 0.0; // Simulated time within SSA

    // All the reactions
    std::deque<std::shared_ptr<ito_reaction> > reactions;
    for (const auto& r : m_reactions){
      reactions.emplace_back(std::make_shared<ito_reaction>(r.second));
    }


    Real u1;     // Uniform random number on [0,1]
    Real A;      // Total propensity
    Real nextDt; // Time to next collision

    while(curDt <= a_dt){

      // Update reaction rates.
      this->update_reaction_rates(a_E, a_mean_energies, a_dx, a_kappa);

      // Get all propensities and compute the total propensity. 
      const std::deque<Real> propensities = this->propensities(a_particles);

      // Total propensity
      A = 0.0;
      for (const auto& a : propensities){
	A += a;
      }

      // Random numebr and next collision time. 
      u1     = m_udist01(m_rng);
      nextDt = log(1./u1)/A;

      // Time after next collision. 
      curDt += nextDt;

      if(curDt <= a_dt){ // Collision within a_dt, do an SSA step. 
	this->ssa_step(a_particles, a_newPhotons, a_mean_energies, a_sources, propensities, reactions, A, a_dx, a_kappa, nextDt);
      }
      else{ // No reaction in this time interval, but we should still increment/decrement the energy with the source.
	const Real dt = nextDt - (curDt - a_dt);
	this->energy_step(a_mean_energies, a_particles, a_sources, a_dx, a_kappa, dt);
      }
    }
  }
}

inline
void ito_plasma_physics::energy_step(Vector<Real>&                     a_mean_energies,
				     const Vector<unsigned long long>& a_particles,
				     const Vector<Real>&               a_sources,
				     const Real                        a_dx,
				     const Real                        a_kappa,
				     const Real                        a_dt) const {

  const Real dV = a_kappa*pow(a_dx, SpaceDim);
  for (int i = 0; i < a_particles.size(); i++){
    if(a_particles[i] > 0){
      a_mean_energies[i] = Max(0.0, a_mean_energies[i] + a_sources[i]*dV*a_dt/a_particles[i]);
    }
    else{
      a_mean_energies[i] = 0.0;
    }
  }
}


inline
void ito_plasma_physics::advance_particles(Vector<unsigned long long>& a_particles,
					   Vector<unsigned long long>& a_newPhotons,
					   Vector<Real>&               a_mean_energies,
					   const Vector<Real>&         a_sources,
					   const Real                  a_dt,
					   const RealVect              a_E,
					   const Real                  a_dx,
					   const Real                  a_kappa) const{
  // Which algorithm?
  if(m_algorithm == algorithm::ssa){
    this->advance_ssa(a_particles, a_newPhotons, a_mean_energies, a_sources, a_dt, a_E, a_dx, a_kappa); 
  }
  else if(m_algorithm == algorithm::tau){
    this->advance_tau(a_particles, a_newPhotons, a_mean_energies, a_sources, a_dt, a_E, a_dx, a_kappa); 
  }
  else if(m_algorithm == algorithm::hybrid){
    MayDay::Abort("ito_plasma_physics::advance_particles_lea - hybrid not implemented");
  }
}

inline
void ito_plasma_physics::advance_hybrid(Vector<unsigned long long>& a_particles,
					Vector<unsigned long long>& a_newPhotons,
					Vector<Real>&               a_mean_energies,
					const Vector<Real>&         a_sources,
					const Real                  a_dt,
					const RealVect              a_E,
					const Real                  a_dx,
					const Real                  a_kappa) const {
  Real curTime = 0.0;

  while(curTime < a_dt){

    // Update reaction rates
    this->update_reaction_rates(a_E, a_mean_energies, a_dx, a_kappa);

    // Get critical and non-critical reactions. 
    ReactionTypes rTypes = this->get_critical_reactions(a_particles);

    ReactionPtrs criticalReactions    = rTypes.first;
    ReactionPtrs nonCriticalReactions = rTypes.second;

    // Compute time steps for critical and non-critical reactions
    Real dtCrit    = this->get_critical_step(a_particles,     criticalReactions);
    Real dtNonCrit = this->get_non_critical_step(a_particles, nonCriticalReactions);

    bool valid = false;
    while(!valid){ 
      
      const Real curDt   = std::min(a_dt-curTime, std::min(dtCrit, dtNonCrit));
      const bool nonCrit = dtNonCrit < dtCrit || criticalReactions.size() == 0 || dtCrit < (a_dt-curTime);
	    
      // How will we advance?
      if(nonCrit){// No critical reactions fire. Advance only noncritical reactions.

	// Compute the total propensity function in case it is faster to run SSA steps
	Real A = this->propensity(a_particles);
	
	if(A*curDt < m_SSAlim){ // Tau-leaping is inefficient, do some SSA on the WHOLE reaction set. 
	  for (int i = 0; i < m_NSSA; i++){

	    // Update reaction rates, propensities, and the total propensity. 
	    const std::deque<Real> propensities = this->propensities(a_particles);
	    A = this->propensity(a_particles); // Recompute propensities
	      
	    // Draw time to next reaction. 
	    const Real u1     = m_udist01(m_rng);
	    const Real dtNext = log(1./u1)/A;

	    valid    = true;
	    curTime += dtNext;
	    if(curTime <= a_dt){ // Ok, we are inside [0, a_dt] and we can do the reaction. 
	      this->ssa_step(a_particles, a_newPhotons, a_mean_energies, a_sources, propensities, A, a_dx, a_kappa, a_dt);
	    }
	    else{ // Next reaction is outside the interval, but we should still increment with the energy.
	      const Real dt = dtNext - (curDt - a_dt);
	      this->energy_step(a_mean_energies, a_particles, a_sources, a_dx, a_kappa, dt);
	      
	      break; // No need to perform the rest. 
	    }

	    // Update reaction rates for next SSA step. 
	    this->update_reaction_rates(a_E, a_mean_energies, a_dx, a_kappa);
	  }
	}
	else{ // Perform tau-leaping with no-critical reactions only
	  // Need a copy in case step is rejected
	  Vector<unsigned long long> particles = a_particles;
	  Vector<unsigned long long> photons   = a_newPhotons;
	  Vector<Real> mean_energies           = a_mean_energies;
	
	  this->tau_step(particles, photons, mean_energies, a_sources, nonCriticalReactions, a_dx, a_kappa, curDt);

	  // Make sure step was valid
	  valid = this->valid_state(particles);
	  if(valid){
	    curTime        += curDt;
	    a_particles     = particles;
	    a_newPhotons    = photons;
	    a_mean_energies = mean_energies;
	  }
	  else{
	    dtNonCrit = dtNonCrit*0.5;
	    continue;
	  }
	}
      } 
      else{ // One critical reaction fired. Advance that reaction and use tau-leaping for the rest.
	// Need a copy in case step is rejected
	Vector<unsigned long long> particles = a_particles;
	Vector<unsigned long long> photons   = a_newPhotons;
	Vector<Real> mean_energies           = a_mean_energies;

	const std::deque<Real> propensitiesC  = this->propensities(particles, criticalReactions);
	const std::deque<Real> propensitiesNC = this->propensities(particles, nonCriticalReactions);
	
	const Real Ac = this->propensity(particles, criticalReactions);

	// SSA for critical reactions
	//	this->ssa_step(particles, photons, propensitiesC, criticalReactions, Ac);
	this->ssa_step(particles, photons, mean_energies, a_sources, propensitiesC, criticalReactions, Ac, a_dx, a_kappa, curDt);

	// Tau-leaping for non-critical reactions
	//	this->tau_step(particles, photons, propensitiesNC, nonCriticalReactions, curDt);
	this->tau_step(particles, photons, mean_energies, a_sources, nonCriticalReactions, a_dx, a_kappa, curDt);

	valid = this->valid_state(particles);

	if(valid){
	  curTime        += curDt;
	  a_particles     = particles;
	  a_newPhotons    = photons;
	  a_mean_energies = mean_energies;
	}
	else{
	  dtNonCrit = dtNonCrit * 0.5;
	  continue;
	}
      }
    }
  }
}

#endif
