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
bool ito_plasma_physics::valid_advance(const Vector<List<ito_particle>* >& a_particles,
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

inline
void ito_plasma_physics::reconcile_particles(Vector<List<ito_particle>* >& a_particles,
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

inline
void ito_plasma_physics::reconcile_particles(List<ito_particle>& a_particles,
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

inline
void ito_plasma_physics::compute_particle_weights(int&      a_weight,
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

inline
void ito_plasma_physics::add_particles(List<ito_particle>& a_particles,
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

inline
void ito_plasma_physics::remove_particles(List<ito_particle>& a_particles,
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

inline
void ito_plasma_physics::reconcile_photons(List<photon>&      a_photons,
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

inline
void ito_plasma_physics::reconcile_photoionization(Vector<List<ito_particle>* >& a_particles,
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

inline
int ito_plasma_physics::poisson_reaction(const Real a_propensity, const Real a_dt) const{
  
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

inline
void ito_plasma_physics::draw_gaussian_particles(List<ito_particle>& a_particles1,
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

inline
Vector<int> ito_plasma_physics::get_particle_count(const Vector<List<ito_particle>* >& a_particles) const {
  Vector<int> ret(a_particles.size(), 0);

  for (int i = 0; i < a_particles.size(); i++){
    ret[i] = this->get_particle_count(*a_particles[i]);
  }

  return ret;
}

inline
int ito_plasma_physics::get_particle_count(const List<ito_particle>& a_particles) const {
  int ret = 0;

  for (ListIterator<ito_particle> lit(a_particles); lit.ok(); ++lit){
    ret += round(lit().mass());
  }

  return ret;
}

inline
void ito_plasma_physics::tau_leap(Vector<int>& a_particles, Vector<int>& a_newPhotons, const Real a_dt) const {

  const Vector<int> particles  = a_particles; // Because propensities should not change. 
  
  for (const auto& r : m_reactions){
    const int num_reactions = this->poisson_reaction(r.second.propensity(particles), a_dt);

    r.second.jump_state(a_particles, a_newPhotons, num_reactions);
  }
}

inline
void ito_plasma_physics::advance_reaction_network(Vector<List<ito_particle>* >& a_particles,
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

inline
void ito_plasma_physics::hybrid_advance(Vector<int>& a_particles, Vector<int>& a_newPhotons, const Real a_dt){
  this->tau_leap(a_particles, a_newPhotons, a_dt);
}

inline
ReactionTypes ito_plasma_physics::get_critical_reactions(const Vector<int>& a_particles) const{

  ReactionPtrs criticalReactions;    // Maybe these should be maps, too?
  ReactionPtrs nonCriticalReactions; // Maybe these should be maps, too?
  
  for (const auto& r : m_reactions){

    const std::map<int, int>& stateChange = r.second.get_state_change();

    // Cheap trick for setting to largest representable integer
    int Lj           = std::numeric_limits<int>::max();
    bool consuming  = false;

    // Check how many times reaction r can fire before consuming all its species
    for (const auto& s: stateChange){

      const int reactant = s.first;
      const int nuIJ     = s.second;

      if(nuIJ < 0){
	consuming = true;
	const int b = (a_particles[reactant] + std::abs(nuIJ) - 1)/std::abs(nuIJ); // Upwards division of Xi/|nuIJ|
	Lj = std::min(Lj, b);
      }
    }

    // Check if reaction is a critical reaction
    bool crit = false;
    if(Lj < m_Ncrit && consuming){
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
Real ito_plasma_physics::get_critical_step(const Vector<int>&  a_particles,
					   const ReactionPtrs& a_critical_reactions) const {
						 
  // TLDR: a_critical_reaction is a vector under the hood.
  
  Real A   = 0.0;
  for (const auto& r : a_critical_reactions){
    A += r->propensity(a_particles);
  }

  const Real u1 = m_udist01(m_rng);
  const Real dt = log(1./u1)/A;

  return dt;
}

inline
Real ito_plasma_physics::get_non_critical_step(const Vector<int>&  a_particles,
					       const ReactionPtrs& a_non_critical_reactions) const {
						 
  // TLDR: a_non_critical_reaction is a vector under the hood.
  
  Real dt = std::numeric_limits<Real>::max();

  //  const int M = a_non_critical_reactions.size();

  //  if(M > 0){
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
	const Real& ap        = R.first;     // Propensity
	const ito_reaction& r = *(R.second); // Reaction object
	const int muIJ        = r.get_state_change(reactant);

	mu    += muIJ * ap;
	sigma += muIJ * muIJ * ap;

      }

      Real dt1 = std::numeric_limits<Real>::max();
      Real dt2 = std::numeric_limits<Real>::max();

      if(mu    != 0.0) dt1 = std::max(m_eps*a_particles[reactant]/gi, 1.)/mu;
      if(sigma != 0.0) dt2 = std::max(m_eps*a_particles[reactant]/gi, 1.)/sigma;

      dt = std::min(dt, std::min(dt1, dt2));
    }

  }


  return dt;
}

inline
Real ito_plasma_physics::total_propensity(const Vector<int>& a_particles) const{

  Real A = 0.0;
  for (const auto& r : m_reactions){
    A += r.second.propensity(a_particles);
  }

  return A;
}


inline
Real ito_plasma_physics::critical_propensity(const Vector<int>& a_particles, const ReactionPtrs& a_critical_reactions) const{

  Real A = 0.0;
  std::deque<Real> propensities = this->critical_propensities(a_particles, a_critical_reactions);
  
  for (const auto& p : propensities){
    A += p;
  }

  return A;
}

inline
std::deque<Real> ito_plasma_physics::critical_propensities(const Vector<int>&  a_particles,
							   const ReactionPtrs& a_critical_reactions) const{

  std::deque<Real> propensities;

  for (const auto& r : a_critical_reactions){
    const Real p = r->propensity(a_particles);
    propensities.emplace_back(p);
  }

  return propensities;
}



inline
void ito_plasma_physics::ssa_step(Vector<int>&             a_particles,
				  Vector<int>&             a_newPhotons,
				  const ReactionPtrs&      a_critical_reactions,
				  const std::deque<Real>&  a_propensities,
				  const Real               a_A) const {


  const Real u2 = m_udist01(m_rng);
  const int M = a_critical_reactions.size();
  
  if(a_critical_reactions.size() > 0){
    int r = 0;

    for (int i = 1; i < M; i++){
      if(a_propensities[i-1] <= u2*a_A && u2*a_A < a_propensities[i]){
	r = i;
	break;
      }
    }

    a_critical_reactions[r]->jump_state(a_particles, a_newPhotons, 1);
  }
}

inline
void ito_plasma_physics::ssa_step(Vector<int>&             a_particles,
				  Vector<int>&             a_newPhotons,
				  const std::deque<Real>&  a_propensities,
				  const Real               a_totalPropensity) const {
  
  // A bit nasty, but we need a data structure with random access for evaluating stuff like
  // a_(r-1) <= u2*A < a_r. So, take the internal map and put it into a vector. 
  std::deque<std::shared_ptr<ito_reaction> > reactions;
  for (const auto& r : m_reactions){
    reactions.emplace_back(std::make_shared<ito_reaction>(r.second));
  }

  this->ssa_step(a_particles, a_newPhotons, reactions, a_propensities, a_totalPropensity);

}

inline
void ito_plasma_physics::tau_step(Vector<int>&             a_particles,
				  Vector<int>&             a_newPhotons,
				  const ReactionPtrs&      a_non_critical_reactions,
				  const std::deque<Real>&  a_propensities,
				  const Real               a_dt) const {

  Vector<int> particles = a_particles;
  for (const auto& r : a_non_critical_reactions){
    const int num_reactions = this->poisson_reaction(r->propensity(particles), a_dt);

    r->jump_state(a_particles, a_newPhotons, num_reactions);
  }
}

inline
void ito_plasma_physics::tau_step(Vector<int>&             a_particles,
				  Vector<int>&             a_newPhotons,
				  const std::deque<Real>&  a_propensities,
				  const Real               a_dt) const {
  // A bit nasty, but we need a data structure with random access for evaluating stuff like
  // a_(r-1) <= u2*A < a_r. So, take the internal map and put it into a vector. 
  std::deque<std::shared_ptr<ito_reaction> > reactions;
  for (const auto& r : m_reactions){
    reactions.emplace_back(std::make_shared<ito_reaction>(r.second));
  }

  this->tau_step(a_particles, a_newPhotons, reactions, a_propensities, a_dt);
}

#endif
