/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaPhysicsImplem.H 
  @brief  Implementation of CD_ItoPlasmaPhysics.H
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaPhysicsImplem_H
#define CD_ItoPlasmaPhysicsImplem_H

// Std includes
#include <algorithm>
#include <unordered_set>

// Chombo includes
#include <PolyGeom.H>

// Our includes
#include <CD_ItoPlasmaPhysics.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

#define ITO_PLASMA_PHYSICS_DEBUG 0

using namespace Physics::ItoPlasma;

inline RealVect
ItoPlasmaPhysics::randomPosition(const RealVect a_cellPos,
                                 const RealVect a_centroidPos,
                                 const RealVect a_lo,
                                 const RealVect a_hi,
                                 const RealVect a_bndryCentroid,
                                 const RealVect a_bndryNormal,
                                 const Real     a_dx,
                                 const Real     a_kappa) const
{

  RealVect pos;

  constexpr Real thresh = 1.E-3;
  if (a_kappa < thresh) {
    pos = a_centroidPos;
  }
  else if (a_kappa < 1.0) { // Rejection sampling.
    pos = this->randomPosition(a_lo, a_hi, a_bndryCentroid, a_bndryNormal);
  }
  else { // Regular cell. Get a position.
    pos = this->randomPosition(a_lo, a_hi);
  }

  pos = a_cellPos + pos * a_dx;

  return pos;
}

inline RealVect
ItoPlasmaPhysics::randomPosition(const RealVect a_lo,
                                 const RealVect a_hi,
                                 const RealVect a_bndryCentroid,
                                 const RealVect a_bndryNormal) const
{

  RealVect pos   = this->randomPosition(a_lo, a_hi);
  bool     valid = PolyGeom::dot(pos - a_bndryCentroid, a_bndryNormal) >= 0.0;

  while (!valid) {
    pos   = this->randomPosition(a_lo, a_hi);
    valid = PolyGeom::dot(pos - a_bndryCentroid, a_bndryNormal) >= 0.0;
  }

  return pos;
}

inline RealVect
ItoPlasmaPhysics::randomPosition(const RealVect a_lo, const RealVect a_hi) const
{

  RealVect pos = RealVect::Unit;

  for (int dir = 0; dir < SpaceDim; dir++) {
    pos[dir] = a_lo[dir] + m_udist01(m_rng) * (a_hi[dir] - a_lo[dir]);
  }

  return pos;
}

inline RealVect
ItoPlasmaPhysics::randomDirection() const
{
  const Real EPS = 1.E-8;
#if CH_SPACEDIM == 2
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1 * x1 + x2 * x2;
  while (r >= 1.0 || r < EPS) {
    x1 = m_udist11(m_rng);
    x2 = m_udist11(m_rng);
    r  = x1 * x1 + x2 * x2;
  }

  return RealVect(x1, x2) / sqrt(r);
#elif CH_SPACEDIM == 3
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1 * x1 + x2 * x2;
  while (r >= 1.0 || r < EPS) {
    x1 = m_udist11(m_rng);
    x2 = m_udist11(m_rng);
    r  = x1 * x1 + x2 * x2;
  }

  const Real x = 2 * x1 * sqrt(1 - r);
  const Real y = 2 * x2 * sqrt(1 - r);
  const Real z = 1 - 2 * r;

  return RealVect(x, y, z);
#endif
}

inline bool
ItoPlasmaPhysics::isStateValid(const Vector<long long>& a_particles) const
{
  bool ret = true;

  for (int i = 0; i < a_particles.size(); i++) {
    if (a_particles[i] < 0) {
      ret = false;
      break;
    }
  }

  return ret;
}

inline void
ItoPlasmaPhysics::checkChargeConservation(const Vector<long long>& a_oldNumParticles,
                                          const Vector<long long>& a_newNumParticles) const
{
  constexpr long long zero = 0;

  long long charge_before = 0;
  long long charge_after  = 0;

  for (int i = 0; i < a_oldNumParticles.size(); i++) {
    const int q = m_ItoSpecies[i]->getChargeNumber();

    charge_before += a_oldNumParticles[i] * q;
    charge_after += a_newNumParticles[i] * q;
  }

  const long long diff = charge_before - charge_after;

  if (diff != zero) {
    MayDay::Abort("ItoPlasmaPhysics::checkChargeConservation() - charge is not conserved!");
  }
}

inline long long
ItoPlasmaPhysics::computeCharge(const Vector<List<ItoParticle>*>& a_particles) const
{

  long long ret = 0;
  for (int i = 0; i < a_particles.size(); i++) {
    const long long q = (long long)m_ItoSpecies[i]->getChargeNumber();

    for (ListIterator<ItoParticle> lit(*a_particles[i]); lit.ok(); ++lit) {
      const long long w = llround(lit().weight());

      ret += w * q;
    }
  }

  return ret;
}

inline void
ItoPlasmaPhysics::reconcileParticles(Vector<List<ItoParticle>*>& a_particles,
                                     const Vector<long long>&    a_newNumParticles,
                                     const Vector<long long>&    a_oldNumParticles,
                                     const RealVect              a_cellPos,
                                     const RealVect              a_centroidPos,
                                     const RealVect              a_lo,
                                     const RealVect              a_hi,
                                     const RealVect              a_bndryCentroid,
                                     const RealVect              a_bndryNormal,
                                     const Real                  a_dx,
                                     const Real                  a_kappa) const
{

#if ITO_PLASMA_PHYSICS_DEBUG
  const bool isStateValid = this->isStateValid(a_newNumParticles);
  if (!isStateValid) {
    MayDay::Abort("ito_plmasma_physics::reconcileParticles - not a valid state!!!");
  }
#endif

  for (int i = 0; i < a_particles.size(); i++) {
    this->reconcileParticles(*a_particles[i],
                             a_newNumParticles[i],
                             a_oldNumParticles[i],
                             a_cellPos,
                             a_centroidPos,
                             a_lo,
                             a_hi,
                             a_bndryCentroid,
                             a_bndryNormal,
                             a_dx,
                             a_kappa);
  }
}

inline void
ItoPlasmaPhysics::reconcileParticles(List<ItoParticle>& a_particles,
                                     const long long    a_newNumParticles,
                                     const long long    a_oldNumParticles,
                                     const RealVect     a_cellPos,
                                     const RealVect     a_centroidPos,
                                     const RealVect     a_lo,
                                     const RealVect     a_hi,
                                     const RealVect     a_bndryCentroid,
                                     const RealVect     a_bndryNormal,
                                     const Real         a_dx,
                                     const Real         a_kappa) const
{

  const long long     diff = a_newNumParticles - a_oldNumParticles;
  constexpr long long zero = 0;

  if (diff > zero) { // Need to add 'diff' mass.
    const long long Diff = (long long)(diff);
    this->addParticles(a_particles,
                       diff,
                       a_cellPos,
                       a_centroidPos,
                       a_lo,
                       a_hi,
                       a_bndryCentroid,
                       a_bndryNormal,
                       a_dx,
                       a_kappa);
  }
  else if (diff < zero) { // Need to remove 'diff' mass.
    this->remove_particles(a_particles, -diff);
  }
}

inline void
ItoPlasmaPhysics::computeParticleWeights(long long&      a_weight,
                                         long long&      a_num,
                                         long long&      a_remainder,
                                         const long long a_numNewParticles) const
{

  if (a_numNewParticles <= m_ppc) {
    a_weight    = 1;
    a_remainder = 0;
    a_num       = a_numNewParticles;
  }
  else { // Add superparticles
    a_weight    = a_numNewParticles / m_ppc;
    a_remainder = a_numNewParticles % m_ppc;
    a_num       = (a_remainder == 0) ? m_ppc : m_ppc - 1;
  }
}

inline void
ItoPlasmaPhysics::addParticles(List<ItoParticle>& a_particles,
                               const long long    a_numNewParticles,
                               const RealVect     a_cellPos,
                               const RealVect     a_centroidPos,
                               const RealVect     a_lo,
                               const RealVect     a_hi,
                               const RealVect     a_bndryCentroid,
                               const RealVect     a_bndryNormal,
                               const Real         a_dx,
                               const Real         a_kappa) const
{
  long long weight, num, remainder;

  constexpr long long zero = 0;

  this->computeParticleWeights(weight, num, remainder, a_numNewParticles);

  for (int i = 0; i < num; i++) {
    const RealVect p =
      this->randomPosition(a_cellPos, a_centroidPos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
    const Real w = (Real)weight;
    a_particles.add(ItoParticle(w, p));
  }

  if (remainder > zero) { // Rest of weight in case we got superparticles
    const RealVect p =
      this->randomPosition(a_cellPos, a_centroidPos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
    const Real w = (Real)(weight + remainder);
    a_particles.add(ItoParticle(w, p));
  }
}

inline void
ItoPlasmaPhysics::remove_particles(List<ItoParticle>& a_particles, const long long a_numParticlesToRemove) const
{

  constexpr long long zero = 0;

  if (a_numParticlesToRemove > zero) {
#if ITO_PLASMA_PHYSICS_DEBUG // Debug
    this->sanityCheck(a_particles);

    const long long weightBefore = this->getTotalWeight(a_particles);
#endif

    ListIterator<ItoParticle> lit(a_particles);
    Vector<long long>         particleWeights(a_particles.length());

    long long minWeight  = std::numeric_limits<long long>::max();
    long long numRemoved = zero;

    // Compute weights and trim all particle weights down to the minimum weight. We shouldn't have weight = 0 when we
    // get here...
    this->computeWeights(lit, particleWeights, minWeight); // Computes weights
    this->trimWeights(numRemoved,
                      particleWeights,
                      a_numParticlesToRemove,
                      minWeight);                                              // Remove weight down to minWeight
    this->removeWeights(particleWeights, a_numParticlesToRemove - numRemoved); // Uniformly removes rest of mass
    this->reconcileWeights(a_particles, particleWeights);                      // Reconcile particle weights

#if ITO_PLASMA_PHYSICS_DEBUG // Debug
    this->sanityCheck(a_particles);
    const long long weightAfter = this->getTotalWeight(a_particles);

    if (weightBefore - weightAfter != a_numParticlesToRemove) {
      std::cout << "weight before = " << weightBefore << std::endl;
      std::cout << "weight after  = " << weightAfter << std::endl;
      std::cout << "Should remove = " << a_numParticlesToRemove << std::endl;
      std::cout << "numRemoved    = " << numRemoved << std::endl;

      MayDay::Abort("ItoPlasmaPhysics::remove_particles - incorrect mass removed");
    }
#endif
  }
}

inline void
ItoPlasmaPhysics::setMeanParticleEnergy(Vector<List<ItoParticle>*>& a_particles, const Vector<Real> a_mean_energy) const
{
  for (int i = 0; i < a_particles.size(); i++) {
    this->setMeanParticleEnergy(*a_particles[i], a_mean_energy[i]);
  }
}

inline void
ItoPlasmaPhysics::setMeanParticleEnergy(List<ItoParticle>& a_particles, const Real a_mean_energy) const
{
#if ITO_PLASMA_PHYSICS_DEBUG
  if (a_mean_energy < 0.0) {
    pout() << "ItoPlasmaPhysics::setMeanParticleEnergy - energy = " << a_mean_energy << endl;
  }
  else if (a_mean_energy != a_mean_energy) {
    MayDay::Abort("ItoPlasmaPhysics::setMeanParticleEnergy - a_mean_energy is NaN!");
  }
#endif

  for (ListIterator<ItoParticle> lit(a_particles); lit.ok(); ++lit) {
    lit().energy() = a_mean_energy;
  }
}

inline void
ItoPlasmaPhysics::sanityCheck(const List<ItoParticle>& a_particles) const
{

  ListIterator<ItoParticle> lit(a_particles);

  for (lit.begin(); lit.ok(); ++lit) {
    const Real w = lit().weight();

    if (w < 1.0)
      MayDay::Abort("ItoPlasmaPhysics::sanityCheck - can't have w < 1!!!!");
  }
}

inline long long
ItoPlasmaPhysics::getTotalWeight(const List<ItoParticle>& a_particles) const
{

  long long ret = 0;
  for (ListIterator<ItoParticle> lit(a_particles); lit.ok(); ++lit) {
    const long long w = (long long)llround(lit().weight());

    ret += w;
  }

  return ret;
}

inline void
ItoPlasmaPhysics::computeWeights(ListIterator<ItoParticle>& a_lit,
                                 Vector<long long>&         a_particleWeights,
                                 long long&                 a_minWeight) const
{
  int i = 0;
  for (a_lit.begin(); a_lit.ok(); ++a_lit, ++i) {
    a_particleWeights[i] = (long long)llround(a_lit().weight());
    a_minWeight          = std::min(a_particleWeights[i], a_minWeight);
  }
}

inline void
ItoPlasmaPhysics::trimWeights(long long&         a_numRemoved,
                              Vector<long long>& a_particleWeights,
                              const long long    a_numParticlesToRemove,
                              const long long    a_minWeight) const
{
  constexpr long long zero = 0;

  for (int i = 0; i < a_particleWeights.size(); i++) {
    const long long diff1 = a_particleWeights[i] - a_minWeight;    // Should always be > 0
    const long long diff2 = a_numParticlesToRemove - a_numRemoved; // Should always be > 0

    const long long r = std::max(zero, std::min(diff1, diff2));

    a_particleWeights[i] -= r;
    a_numRemoved += r;
  }
}

inline void
ItoPlasmaPhysics::removeWeights(Vector<long long>& a_particleWeights, const long long a_numParticlesToRemove) const
{

  constexpr long long zero = 0;

  if (a_numParticlesToRemove > zero) {
    const long long numCompParticles = (long long)a_particleWeights.size();
    const long long weight           = a_numParticlesToRemove / numCompParticles;
    const long long remainder        = a_numParticlesToRemove % numCompParticles;

    // Remove weight from each particle.
    if (weight > zero) {
      for (int i = 0; i < a_particleWeights.size(); i++) {
        a_particleWeights[i] -= weight;
      }
    }

    // Remove the remainder by cycling through the particles.
    if (remainder > zero) {
      long long W = zero;
      for (int i = 0; i < a_particleWeights.size(); i++) {
        const long long w = std::min(a_particleWeights[i], remainder - W);

        a_particleWeights[i] -= w;
        W += w;

        if (W == remainder) { // OK, we have removed the remainder. Break out.
          break;
        }
      }
    }
  }
}

inline void
ItoPlasmaPhysics::reconcileWeights(List<ItoParticle>& a_particles, const Vector<long long>& a_particleWeights) const
{

  constexpr long long zero = 0;

  List<ItoParticle> copyParticles(a_particles);
  a_particles.clear();

  int i = 0;
  for (ListIterator<ItoParticle> lit(copyParticles); lit.ok(); ++lit, ++i) {
    ItoParticle& p = lit();

    if (a_particleWeights[i] > zero) {
      p.weight() = (Real)a_particleWeights[i];
      a_particles.add(p);
    }
  }
}

void
ItoPlasmaPhysics::reconcilePhotons(Vector<List<Photon>*>&   a_newPhotons,
                                   const Vector<long long>& a_numNewPhotons,
                                   const RealVect           a_cellPos,
                                   const RealVect           a_centroidPos,
                                   const RealVect           a_lo,
                                   const RealVect           a_hi,
                                   const RealVect           a_bndryCentroid,
                                   const RealVect           a_bndryNormal,
                                   const Real               a_dx,
                                   const Real               a_kappa) const
{

  for (int i = 0; i < a_newPhotons.size(); i++) {
    if (a_numNewPhotons[i] > 0) {
      this->reconcilePhotons(*a_newPhotons[i],
                             *m_rtSpecies[i],
                             a_numNewPhotons[i],
                             a_cellPos,
                             a_centroidPos,
                             a_lo,
                             a_hi,
                             a_bndryCentroid,
                             a_bndryNormal,
                             a_dx,
                             a_kappa);
    }
  }
}

inline void
ItoPlasmaPhysics::reconcilePhotons(List<Photon>&    a_Photons,
                                   const RtSpecies& a_species,
                                   const long long  a_num_Photons,
                                   const RealVect   a_cellPos,
                                   const RealVect   a_centroidPos,
                                   const RealVect   a_lo,
                                   const RealVect   a_hi,
                                   const RealVect   a_bndryCentroid,
                                   const RealVect   a_bndryNormal,
                                   const Real       a_dx,
                                   const Real       a_kappa) const
{
  a_Photons.clear();

  long long weight, num, remainder;
  this->computeParticleWeights(weight, num, remainder, a_num_Photons);

  for (int i = 0; i < num; i++) {
    const RealVect p =
      this->randomPosition(a_cellPos, a_centroidPos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
    const RealVect v = Units::c * randomDirection();

    a_Photons.add(Photon(p, v, a_species.getAbsorptionCoefficient(p), weight));
  }

  // If we used superPhotons the last Photon gets some extra oomph.
  if (remainder > 0) {
    const RealVect p =
      this->randomPosition(a_cellPos, a_centroidPos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
    const RealVect v = Units::c * randomDirection();
    a_Photons.add(Photon(p, v, a_species.getAbsorptionCoefficient(p), weight + remainder));
  }
}

inline void
ItoPlasmaPhysics::reconcilePhotoionization(Vector<List<ItoParticle>*>&  a_particles,
                                           Vector<Real>&                a_mean_energies,
                                           const Vector<long long>&     a_numParticles,
                                           const Vector<List<Photon>*>& a_Photons) const
{

  // TLDR: This routine runs through all the photoreactions. For each reaction we look up the target species (e.g. electrons, ions)
  //       and compute the total energy for those species. We then add the photoreaction products, including the added energy, and then
  //       we finally recompute the mean energy.

  //
  Vector<long long> newNumParticles = a_numParticles;
  Vector<Real>      energies(a_mean_energies);
  for (int i = 0; i < energies.size(); i++) {
    energies[i] *= newNumParticles[i];
  }

  for (const auto& r : m_photoReactions) {

    // Source and targets
    const int&             src     = r.second.source();
    const std::deque<int>& targets = r.second.targets();

    const std::deque<std::pair<int, Real>>& target_energies = r.second.getTargetEnergies();

    for (ListIterator<Photon> lit(*a_Photons[src]); lit.ok(); ++lit) {
      const Photon&  phot = lit();
      const RealVect pos  = phot.position();
      const Real     mass = phot.weight(); // Basically the number of reactions, in real format.

      // Increment energies and update total particles.
      for (const auto& t : targets) {
        a_particles[t]->add(ItoParticle(mass, pos));

        newNumParticles[t] += llround(mass); // Now have this many particles in the cell.
      }

      // Increment energies
      for (const auto& te : target_energies) {
        const int  idx = te.first;
        const Real dE  = te.second;
        energies[idx] += mass * dE;
      }
    }
  }

  // Recompute the mean energies. If dE < a_mean_energies then the addition of Photons lead to cooling, pretty cool :)
  for (int i = 0; i < a_mean_energies.size(); i++) {
    if (a_numParticles[i] > 0) {
      a_mean_energies[i] = Max(0.0, energies[i] / newNumParticles[i]);
    }
    else {
      a_mean_energies[i] = 0.0;
    }
  }
}

inline long long
ItoPlasmaPhysics::poissonReaction(const Real a_propensity, const Real a_dt) const
{

  constexpr long long zero  = 0;
  long long           value = 0;
  const Real          mean  = a_propensity * a_dt;

  if (mean < m_fieldSolver_switch) {
    std::poisson_distribution<long long> dist(mean);
    value = dist(m_rng);
  }
  else {
    std::normal_distribution<double> dist(mean, sqrt(mean));
    value = (long long)llround(dist(m_rng));
  }

  return Max(zero, value);
}

inline int
ItoPlasmaPhysics::computeParticlesPerRank(const int a_numParticles) const
{
  const int quotient  = a_numParticles / numProc();
  const int remainder = a_numParticles % numProc();

  Vector<int> particlesPerRank(numProc(), quotient);

  for (int i = 0; i < remainder; i++) {
    particlesPerRank[i] += 1;
  }

  return particlesPerRank[procID()];
}

inline void
ItoPlasmaPhysics::drawGaussianParticles(List<ItoParticle>& a_particles1,
                                        List<ItoParticle>& a_particles2,
                                        const int          a_numParticles,
                                        const RealVect     a_blobCenter,
                                        const Real         a_blobRadius,
                                        const Real         a_weight,
                                        const Real         a_energy1,
                                        const Real         a_energy2)
{

  // Set up the Gaussian RNG
  std::normal_distribution<Real> gauss(0.0, a_blobRadius);

  const int myNumParticles = this->computeParticlesPerRank(a_numParticles);

  // Now make the particles
  for (int i = 0; i < myNumParticles; i++) {
    const RealVect pos = a_blobCenter + gauss(m_rng) * randomDirection();

    a_particles1.add(ItoParticle(a_weight, pos, RealVect::Zero, 0.0, 0.0, a_energy1));
    a_particles2.add(ItoParticle(a_weight, pos, RealVect::Zero, 0.0, 0.0, a_energy2));
  }
}

inline void
ItoPlasmaPhysics::drawSphereParticles(List<ItoParticle>& a_particles1,
                                      List<ItoParticle>& a_particles2,
                                      const int          a_numParticles,
                                      const RealVect     a_blobCenter,
                                      const Real         a_blobRadius,
                                      const Real         a_weight,
                                      const Real         a_energy1,
                                      const Real         a_energy2)
{

  const int myNumParticles = this->computeParticlesPerRank(a_numParticles);

  // Now make the particles
  for (int i = 0; i < myNumParticles; i++) {
    bool     accept = false;
    RealVect pos    = RealVect::Zero;

    while (!accept) {
      for (int d = 0; d < SpaceDim; d++) {
        pos[d] = m_udist11(m_rng);
      }

      if (pos.vectorLength() <= 1.0) {
        accept = true;
        pos    = pos * a_blobRadius + a_blobCenter;
      }
    }

    // Add the particle
    a_particles1.add(ItoParticle(a_weight, pos, RealVect::Zero, 0.0, 0.0, a_energy1));
    a_particles2.add(ItoParticle(a_weight, pos, RealVect::Zero, 0.0, 0.0, a_energy2));
  }
}

inline void
ItoPlasmaPhysics::drawBoxParticles(List<ItoParticle>& a_particles1,
                                   List<ItoParticle>& a_particles2,
                                   const int          a_numParticles,
                                   const RealVect     a_boxLo,
                                   const RealVect     a_boxHi,
                                   const Real         a_weight,
                                   const Real         a_energy1,
                                   const Real         a_energy2)
{

  const int myNumParticles = this->computeParticlesPerRank(a_numParticles);

  const RealVect diff = a_boxHi - a_boxLo;
  for (int i = 0; i < myNumParticles; i++) {
    const RealVect pos = a_boxLo + m_udist01(m_rng) * diff;

    a_particles1.add(ItoParticle(a_weight, pos, RealVect::Zero, 0.0, 0.0, a_energy1));
    a_particles2.add(ItoParticle(a_weight, pos, RealVect::Zero, 0.0, 0.0, a_energy2));
  }
}

inline Vector<long long>
ItoPlasmaPhysics::getParticleCount(const Vector<List<ItoParticle>*>& a_particles) const
{
  Vector<long long> ret(a_particles.size(), 0);

  for (int i = 0; i < a_particles.size(); i++) {
    ret[i] = this->getParticleCount(*a_particles[i]);
  }

  return ret;
}

inline long long
ItoPlasmaPhysics::getParticleCount(const List<ItoParticle>& a_particles) const
{
  long long ret = 0;

  for (ListIterator<ItoParticle> lit(a_particles); lit.ok(); ++lit) {
    ret += (long long)llround(lit().weight());
  }

  return ret;
}

inline Vector<Real>
ItoPlasmaPhysics::getParticleEnergy(const Vector<List<ItoParticle>*>& a_particles) const
{
  Vector<Real> ret(a_particles.size(), 0.0);

  for (int i = 0; i < a_particles.size(); i++) {
    ret[i] = this->getParticleEnergy(*a_particles[i]);
  }

  return ret;
}

inline Real
ItoPlasmaPhysics::getParticleEnergy(const List<ItoParticle>& a_particles) const
{
  Real ret = 0;

  for (ListIterator<ItoParticle> lit(a_particles); lit.ok(); ++lit) {
    const Real& m = lit().weight();
    const Real& e = lit().energy();

    ret += m * e;
  }

  return ret;
}

inline Vector<Real>
ItoPlasmaPhysics::getParticleMeanEnergy(const Vector<Real>& a_energies, const Vector<long long>& a_particle_count) const
{
  Vector<Real> ret(a_energies.size(), 0.0);

  for (int i = 0; i < a_energies.size(); i++) {
    ret[i] = this->getParticleMeanEnergy(a_energies[i], a_particle_count[i]);
  }

  return ret;
}

inline Real
ItoPlasmaPhysics::getParticleMeanEnergy(const Real& a_energy, const long long& a_particle_count) const
{
  const Real ret = (a_particle_count > 0) ? a_energy / a_particle_count : 0.0;
  return ret;
}

inline void
ItoPlasmaPhysics::advanceReactionNetwork(Vector<List<ItoParticle>*>& a_particles,
                                         Vector<List<Photon>*>&      a_Photons,
                                         Vector<List<Photon>*>&      a_newPhotons,
                                         const Vector<Real>&         a_sources,
                                         const RealVect              a_E,
                                         const RealVect              a_cellPos,
                                         const RealVect              a_centroidPos,
                                         const RealVect              a_bndryCentroid,
                                         const RealVect              a_bndryNormal,
                                         const RealVect              a_lo,
                                         const RealVect              a_hi,
                                         const Real                  a_dx,
                                         const Real                  a_kappa,
                                         const Real                  a_dt)
{
#if ITO_PLASMA_PHYSICS_DEBUG
  const long long Qbefore = this->computeCharge(a_particles);
#endif

  // Get particle counts.
  Vector<long long> newPhotonCount   = Vector<long long>(m_numRtSpecies, 0);
  Vector<long long> oldParticleCount = this->getParticleCount(a_particles);
  Vector<long long> newParticleCount = oldParticleCount;

  // Energies, mean_energies, and sources. Make sure LFA gives dummy arguments.
  Vector<Real> sources, energies, mean_energies;
  if (m_coupling == coupling::LFA) {
    energies.resize(m_num_ItoSpecies, 0.0);
    mean_energies.resize(m_num_ItoSpecies, 0.0);
    sources.resize(m_num_ItoSpecies, 0.0);
  }
  else {
    energies      = this->getParticleEnergy(a_particles);
    mean_energies = this->getParticleMeanEnergy(energies, oldParticleCount);
    sources       = a_sources;

    // // The one that comes in through ItoPlasmaStepper is in Joule/(dV*dt), make it eV/(cell*dt)
    const Real dV = a_kappa * pow(a_dx, SpaceDim);
    for (int i = 0; i < sources.size(); i++) {
      sources[i] *= dV / Units::Qe;
    }
  }

#if ITO_PLASMA_PHYSICS_DEBUG
  for (int i = 0; i < a_sources.size(); i++) {
    if (a_sources[i] != a_sources[i]) {
      MayDay::Abort("ItoPlasmaPhysics::advanceReactionNetwork - got nan source!");
    }
    if (mean_energies[i] != mean_energies[i]) {
      MayDay::Abort("ItoPlasmaPhysics::advanceReactionNetwork - got nan energy!");
    }
  }
#endif

#if 0 // Debug, set sources to 0
  for (int i = 0; i < sources.size(); i++){
    sources[i] = 0.0;
  }
#endif

  // Do particle advancement
  this->advanceParticles(newParticleCount, newPhotonCount, mean_energies, sources, a_dt, a_E, a_dx, a_kappa);

  // Reconcile particles
  this->reconcileParticles(a_particles,
                           newParticleCount,
                           oldParticleCount,
                           a_cellPos,
                           a_centroidPos,
                           a_lo,
                           a_hi,
                           a_bndryCentroid,
                           a_bndryNormal,
                           a_dx,
                           a_kappa);

  // Generate new Photons
  this->reconcilePhotons(a_newPhotons,
                         newPhotonCount,
                         a_cellPos,
                         a_centroidPos,
                         a_lo,
                         a_hi,
                         a_bndryCentroid,
                         a_bndryNormal,
                         a_dx,
                         a_kappa);

  // Run photo-reactions and increment the number of particles as needed
  this->reconcilePhotoionization(a_particles, mean_energies, newParticleCount, a_Photons);

  // Set new mean energies for all the particles.
  this->setMeanParticleEnergy(a_particles, mean_energies);

#if ITO_PLASMA_PHYSICS_DEBUG // Debug code
  this->checkChargeConservation(oldParticleCount, newParticleCount);

  const long long Qafter = this->computeCharge(a_particles);

  constexpr long long zero = 0;
  if (Qbefore - Qafter != zero) {
    std::cout << "Qbefore = " << Qbefore << std::endl;
    std::cout << "Qafter  = " << Qafter << std::endl;
    MayDay::Abort("ItoPlasmaPhysics::advanceReactionNetwork_lfa(...) - charge is not conserved!!!");
  }
#endif
}

inline ReactionTypes
ItoPlasmaPhysics::getCriticalReactions(const Vector<long long>& a_particles) const
{

  ReactionPtrs criticalReactions;    // std::deque<shared_ptr<ItoPlasmaReaction> > under the hood
  ReactionPtrs nonCriticalReactions; // std::deque<shared_ptr<ItoPlasmaReaction> > under the hood

  const long long Ncrit = (long long)m_Ncrit; // Need long long version of this because because a_particles is long long

  for (const auto& r : m_reactions) {

    const std::map<int, int>& stateChange = r.second.getStateChange();

    long long Lj        = std::numeric_limits<long long>::max();
    bool      consuming = false;

    // Check how many times reaction r can fire before consuming all its species
    for (const auto& s : stateChange) {

      const int reactant = s.first;
      const int nuIJ     = s.second;

      if (nuIJ < 0) {
        consuming = true;

        const long long s = (long long)std::abs(nuIJ);
        const long long b = (long long)(a_particles[reactant] + s - 1) / s;

        Lj = std::min(Lj, b);
      }
    }

    // Check if reaction is a critical reaction
    bool crit = false;
    if (Lj < Ncrit && consuming) {
      crit = true;
    }

    // Move reactions to where they should be
    if (crit) {
      criticalReactions.emplace_back(std::make_shared<ItoPlasmaReaction>(r.second));
    }
    else {
      nonCriticalReactions.emplace_back(std::make_shared<ItoPlasmaReaction>(r.second));
    }
  }

  // Return the critical and non-critical reactions
  ReactionTypes ret = ReactionTypes();
  ret.first         = criticalReactions;
  ret.second        = nonCriticalReactions;

  return ret;
}

inline Real
ItoPlasmaPhysics::getCriticalTimeStep(const Vector<long long>& a_particles,
                                      const ReactionPtrs&      a_critical_reactions) const
{

  // TLDR: a_critical_reaction is a vector under the hood.
  Real dt = std::numeric_limits<Real>::max();

  if (a_critical_reactions.size() > 0) {
    Real A = 0.0;
    for (const auto& r : a_critical_reactions) {
      A += r->propensity(a_particles);
    }

    const Real u1 = m_udist01(m_rng);
    dt            = log(1. / u1) / A;
  }

  return dt;
}

inline Real
ItoPlasmaPhysics::getNonCriticalTimeStep(const Vector<long long>& a_particles,
                                         const ReactionPtrs&      a_non_critical_reactions) const
{

  // TLDR: a_non_critical_reaction is a vector under the hood.

  Real dt = std::numeric_limits<Real>::max();

  if (a_non_critical_reactions.size() > 0) {

    std::deque<Real> propensities;
    for (const auto& r : a_non_critical_reactions) {
      const Real p = r->propensity(a_particles);

      propensities.emplace_back(p);
    }

    // 1. Compute propensities of all reactions. It is possible for propensities to be zero. Order matters,
    //    so use a deque
    std::deque<std::pair<Real, std::shared_ptr<ItoPlasmaReaction>>> reactions;
    for (const auto& r : a_non_critical_reactions) {
      const Real p = r->propensity(a_particles);

      reactions.emplace_back(std::pair<Real, std::shared_ptr<ItoPlasmaReaction>>(p, r));
    }

    // 2. Get the list of ALL reactants
    std::unordered_set<int> allReactants;
    //    for (const auto& reaction : m_reactions){
    for (const auto& reaction : a_non_critical_reactions) {

      const auto& cur_reactants = reaction->getReactants();
      for (const auto& l : cur_reactants) {
        allReactants.insert(l);
      }
    }

    // 3. Iterate through all reactants and compute deviations
    for (const auto& reactant : allReactants) {
      if (a_particles[reactant] > 0) {
        Real mu     = 0.0;
        Real sigma2 = 0.0;

        constexpr int gi = 1; // Set gi to 1 for now, should think about this later.

        for (const auto& R : reactions) {
          const Real&              ap = R.first;                    // Reaction propensity.
          const ItoPlasmaReaction& r  = *(R.second);                // Reaction.
          const int muIJ              = r.getStateChange(reactant); // State change of reactant due to current reaction

          mu += muIJ * ap;
          sigma2 += muIJ * muIJ * ap;
        }

        Real dt1 = std::numeric_limits<Real>::max();
        Real dt2 = std::numeric_limits<Real>::max();

        const Real f = std::max(m_eps * a_particles[reactant] / gi, 1.0);
        if (mu != 0.0) {
          dt1 = f / std::abs(mu);
        }
        if (sigma2 != 0.0) {
          dt2 = (f * f) / std::abs(sigma2);
        }

        dt = std::min(dt, std::min(dt1, dt2));
      }
    }
  }

  return dt;
}

inline Real
ItoPlasmaPhysics::propensity(const Vector<long long>& a_particles) const
{

  Real A = 0.0;
  for (const auto& r : m_reactions) {
    A += r.second.propensity(a_particles);
  }

  return A;
}

inline Real
ItoPlasmaPhysics::propensity(const Vector<long long>& a_particles, const ReactionPtrs& a_critical_reactions) const
{

  Real             A            = 0.0;
  std::deque<Real> propensities = this->propensities(a_particles, a_critical_reactions);

  for (const auto& p : propensities) {
    A += p;
  }

  return A;
}

inline std::deque<Real>
ItoPlasmaPhysics::propensities(const Vector<long long>& a_particles) const
{
  std::deque<Real> propensities;

  for (const auto& r : m_reactions) {
    const Real p = r.second.propensity(a_particles);
    propensities.emplace_back(p);
  }

  return propensities;
}

inline std::deque<Real>
ItoPlasmaPhysics::propensities(const Vector<long long>& a_particles, const ReactionPtrs& a_critical_reactions) const
{
  std::deque<Real> propensities;

  for (const auto& r : a_critical_reactions) {
    const Real p = r->propensity(a_particles);
    propensities.emplace_back(p);
  }

  return propensities;
}

inline void
ItoPlasmaPhysics::updateReactionRates(const RealVect     a_E,
                                      const Vector<Real> a_mean_energies,
                                      const Real         a_dx,
                                      const Real         a_kappa) const
{

  // Update the rate based on the mean energies or the electric field, depending on how we use the LFA/LEA.
  if (m_coupling == coupling::LFA) {
    this->updateReactionRatesLFA(a_E, a_dx, a_kappa);
  }
  else if (m_coupling == coupling::LEA) {
    this->updateReactionRatesLEA(a_E, a_mean_energies, a_dx, a_kappa);
  }
}

inline void
ItoPlasmaPhysics::stepTau(Vector<long long>&  a_particles,
                          Vector<long long>&  a_newPhotons,
                          Vector<Real>&       a_mean_energies,
                          const Vector<Real>& a_sources,
                          const ReactionPtrs& a_reactions,
                          const Real          a_dx,
                          const Real          a_kappa,
                          const Real          a_dt) const
{

  // These are the propensities.
  const std::deque<Real> propensities = this->propensities(a_particles, a_reactions);

  // Call the other version.
  this->stepTau(a_particles, a_newPhotons, a_mean_energies, a_sources, a_reactions, propensities, a_dx, a_kappa, a_dt);
}

inline void
ItoPlasmaPhysics::stepTau(Vector<long long>&      a_particles,
                          Vector<long long>&      a_newPhotons,
                          Vector<Real>&           a_mean_energies,
                          const Vector<Real>&     a_sources,
                          const ReactionPtrs&     a_reactions,
                          const std::deque<Real>& a_propensities,
                          const Real              a_dx,
                          const Real              a_kappa,
                          const Real              a_dt) const
{
  // mean_energies -> total energies
  Vector<Real> energies(m_num_ItoSpecies);
  for (int i = 0; i < a_particles.size(); i++) {
    energies[i] = Max(0.0, a_mean_energies[i] * a_particles[i]);
  }

  // Go through list of reactions.
  for (int i = 0; i < a_reactions.size(); i++) {
    const long long num_reactions = this->poissonReaction(a_propensities[i], a_dt);

    // Jump state and energy
    a_reactions[i]->jumpState(a_particles, a_newPhotons, num_reactions); //
    a_reactions[i]->jumpEnergy(energies, a_mean_energies, a_sources, num_reactions, a_dt);
  }

  // Recompute mean energies afterwards.
  const Real dV = a_kappa * pow(a_dx, SpaceDim);

  for (int i = 0; i < a_particles.size(); i++) {
    if (a_particles[i] > 0LL) {
      a_mean_energies[i] = Max(0.0, energies[i] / (1.0 * a_particles[i]));
    }
    else {
      a_mean_energies[i] = 0.0;
    }
  }
}

inline void
ItoPlasmaPhysics::stepTau(Vector<long long>&  a_particles,
                          Vector<long long>&  a_newPhotons,
                          Vector<Real>&       a_mean_energies,
                          const Vector<Real>& a_sources,
                          const Real          a_dx,
                          const Real          a_kappa,
                          const Real          a_dt) const
{

  std::deque<std::shared_ptr<ItoPlasmaReaction>> reactions;
  for (const auto& r : m_reactions) {
    reactions.emplace_back(std::make_shared<ItoPlasmaReaction>(r.second));
  }

  this->stepTau(a_particles, a_newPhotons, a_mean_energies, a_sources, reactions, a_dx, a_kappa, a_dt);
}

inline void
ItoPlasmaPhysics::advanceTau(Vector<long long>&  a_particles,
                             Vector<long long>&  a_newPhotons,
                             Vector<Real>&       a_mean_energies,
                             const Vector<Real>& a_sources,
                             const Real          a_dt,
                             const RealVect      a_E,
                             const Real          a_dx,
                             const Real          a_kappa) const
{
  // Update reaction rates in case we use LEA.
  this->updateReactionRates(a_E, a_mean_energies, a_dx, a_kappa);

  // Do a tau step. Use rejection sampling so we only get valid states. But don't update reaction rates.
  Real curTime = 0.0;
  Real nextDt  = a_dt;

  while (curTime < a_dt) {

    // Try big step first.
    nextDt = a_dt - curTime;

    bool valid = false;
    while (!valid) {

      Vector<long long> particles     = a_particles;
      Vector<long long> Photons       = a_newPhotons;
      Vector<Real>      mean_energies = a_mean_energies;

      // Do the tau-leaping
      this->stepTau(particles, Photons, mean_energies, a_sources, a_dx, a_kappa, nextDt);

      // Accept step if valid, else, reduce dt.
      valid = this->isStateValid(particles);
      if (valid) {
        a_particles     = particles;
        a_newPhotons    = Photons;
        a_mean_energies = mean_energies;

        curTime += nextDt;

        // Update reaction rates in case we use LEA.
        this->updateReactionRates(a_E, a_mean_energies, a_dx, a_kappa);
      }
      else {
        nextDt = nextDt * 0.5;
      }
    }
  }
}

inline void
ItoPlasmaPhysics::stepSSA(Vector<long long>&  a_particles,
                          Vector<long long>&  a_newPhotons,
                          Vector<Real>&       a_mean_energies,
                          const Vector<Real>& a_sources,
                          const Real          a_dx,
                          const Real          a_kappa,
                          const Real          a_dt) const
{

  std::deque<std::shared_ptr<ItoPlasmaReaction>> reactions;
  for (const auto& r : m_reactions) {
    reactions.emplace_back(std::make_shared<ItoPlasmaReaction>(r.second));
  }

  // Compute all propensities
  Real             A = 0.0;
  std::deque<Real> propensities;
  for (const auto& r : reactions) {
    const Real a = r->propensity(a_particles);
    propensities.emplace_back(a);
    A += a;
  }

  this->stepSSA(a_particles, a_newPhotons, a_mean_energies, a_sources, propensities, reactions, A, a_dx, a_kappa, a_dt);
}

inline void
ItoPlasmaPhysics::stepSSA(Vector<long long>&      a_particles,
                          Vector<long long>&      a_newPhotons,
                          Vector<Real>&           a_mean_energies,
                          const Vector<Real>&     a_sources,
                          const std::deque<Real>& a_propensities,
                          const Real              a_A,
                          const Real              a_dx,
                          const Real              a_kappa,
                          const Real              a_dt) const
{

  std::deque<std::shared_ptr<ItoPlasmaReaction>> reactions;
  for (const auto& r : m_reactions) {
    reactions.emplace_back(std::make_shared<ItoPlasmaReaction>(r.second));
  }

  this->stepSSA(a_particles,
                a_newPhotons,
                a_mean_energies,
                a_sources,
                a_propensities,
                reactions,
                a_A,
                a_dx,
                a_kappa,
                a_dt);
}

inline void
ItoPlasmaPhysics::stepSSA(Vector<long long>&      a_particles,
                          Vector<long long>&      a_newPhotons,
                          Vector<Real>&           a_mean_energies,
                          const Vector<Real>&     a_sources,
                          const std::deque<Real>& a_propensities,
                          const ReactionPtrs&     a_reactions,
                          const Real              a_A,
                          const Real              a_dx,
                          const Real              a_kappa,
                          const Real              a_dt) const
{

  const int M = a_reactions.size();

  constexpr long long one = 1LL;

  if (M > 0) {

    // Mean_energies -> total energies
    Vector<Real> energies(a_mean_energies);
    for (int i = 0; i < a_particles.size(); i++) {
      energies[i] *= a_particles[i];
    }

    // RNG.
    const Real u2 = m_udist01(m_rng);

    // Determine the reaction type.
    int r = 0;
    for (int i = 1; i < M; i++) {
      if (a_propensities[i - 1] <= u2 * a_A && u2 * a_A < a_propensities[i]) {
        r = i;
        break;
      }
    }

    // Jump state and energy
    a_reactions[r]->jumpState(a_particles, a_newPhotons, one);
    a_reactions[r]->jumpEnergy(energies, a_mean_energies, a_sources, one, a_dt);

    // Recompute mean energies afterwards.
    for (int i = 0; i < a_particles.size(); i++) {
      if (a_particles[i] > 0LL) {
        a_mean_energies[i] = Max(0.0, energies[i] / (1.0 * a_particles[i]));
      }
      else {
        a_mean_energies[i] = 0.0;
      }
    }
  }
}

inline void
ItoPlasmaPhysics::advanceSSA(Vector<long long>&  a_particles,
                             Vector<long long>&  a_newPhotons,
                             Vector<Real>&       a_mean_energies,
                             const Vector<Real>& a_sources,
                             const Real          a_dt,
                             const RealVect      a_E,
                             const Real          a_dx,
                             const Real          a_kappa) const
{

  constexpr long long one = 1;

  const int M = m_reactions.size();

  if (M > 0) {
    Real curDt = 0.0; // Simulated time within SSA

    // All the reactions
    std::deque<std::shared_ptr<ItoPlasmaReaction>> reactions;
    for (const auto& r : m_reactions) {
      reactions.emplace_back(std::make_shared<ItoPlasmaReaction>(r.second));
    }

    Real u1;     // Uniform random number on [0,1]
    Real A;      // Total propensity
    Real nextDt; // Time to next collision

    while (curDt <= a_dt) {

      // Update reaction rates.
      this->updateReactionRates(a_E, a_mean_energies, a_dx, a_kappa);

      // Get all propensities and compute the total propensity.
      const std::deque<Real> propensities = this->propensities(a_particles);

      // Total propensity
      A = 0.0;
      for (const auto& a : propensities) {
        A += a;
      }

      // Random numebr and next collision time.
      u1     = m_udist01(m_rng);
      nextDt = log(1. / u1) / A;

      if (curDt + nextDt <= a_dt) { // Collision within a_dt, do an SSA step.
        curDt += nextDt;
        this->stepSSA(a_particles,
                      a_newPhotons,
                      a_mean_energies,
                      a_sources,
                      propensities,
                      reactions,
                      A,
                      a_dx,
                      a_kappa,
                      nextDt);
      }
      else { // No reaction in this time interval, but we should still increment/decrement the energy with the source.
        curDt += nextDt;
        this->stepEnergy(a_mean_energies, a_particles, a_sources, a_dx, a_kappa, a_dt - curDt);
      }
    }
  }
}

inline void
ItoPlasmaPhysics::stepEnergy(Vector<Real>&            a_mean_energies,
                             const Vector<long long>& a_particles,
                             const Vector<Real>&      a_sources,
                             const Real               a_dx,
                             const Real               a_kappa,
                             const Real               a_dt) const
{
  constexpr long long zero = 0;

  for (int i = 0; i < a_particles.size(); i++) {
    if (a_particles[i] > zero) {
      a_mean_energies[i] = Max(0.0, a_mean_energies[i] + a_sources[i] * a_dt / a_particles[i]);
    }
    else {
      a_mean_energies[i] = 0.0;
    }

#if ITO_PLASMA_PHYSICS_DEBUG
    if (a_mean_energies[i] != a_mean_energies[i]) {
      MayDay::Abort("ItoPlasmaPhysics::stepEnergy - got NaN energy...");
    }
#endif
  }
}

inline void
ItoPlasmaPhysics::advanceParticles(Vector<long long>&  a_particles,
                                   Vector<long long>&  a_newPhotons,
                                   Vector<Real>&       a_mean_energies,
                                   const Vector<Real>& a_sources,
                                   const Real          a_dt,
                                   const RealVect      a_E,
                                   const Real          a_dx,
                                   const Real          a_kappa) const
{
  // // Which algorithm?
  if (m_algorithm == algorithm::ssa) {
    this->advanceSSA(a_particles, a_newPhotons, a_mean_energies, a_sources, a_dt, a_E, a_dx, a_kappa);
  }
  else if (m_algorithm == algorithm::tau) {
    this->advanceTau(a_particles, a_newPhotons, a_mean_energies, a_sources, a_dt, a_E, a_dx, a_kappa);
  }
  else if (m_algorithm == algorithm::hybrid) {
    this->advanceHybrid(a_particles, a_newPhotons, a_mean_energies, a_sources, a_dt, a_E, a_dx, a_kappa);
  }
}

inline void
ItoPlasmaPhysics::advanceHybrid(Vector<long long>&  a_particles,
                                Vector<long long>&  a_newPhotons,
                                Vector<Real>&       a_mean_energies,
                                const Vector<Real>& a_sources,
                                const Real          a_dt,
                                const RealVect      a_E,
                                const Real          a_dx,
                                const Real          a_kappa) const
{
  Real curTime = 0.0;

  while (curTime < a_dt) {

    // Update reaction rates
    this->updateReactionRates(a_E, a_mean_energies, a_dx, a_kappa);

    // Get critical and non-critical reactions.
    ReactionTypes rTypes = this->getCriticalReactions(a_particles);

    ReactionPtrs criticalReactions    = rTypes.first;
    ReactionPtrs nonCriticalReactions = rTypes.second;

    // Compute time steps for critical and non-critical reactions
    Real dtCrit    = this->getCriticalTimeStep(a_particles, criticalReactions);
    Real dtNonCrit = this->getNonCriticalTimeStep(a_particles, nonCriticalReactions);

    bool valid = false;
    while (!valid) {

      const Real curDt   = std::min(a_dt - curTime, std::min(dtCrit, dtNonCrit));
      const bool nonCrit = dtNonCrit < dtCrit || criticalReactions.size() == 0 || dtCrit > (a_dt - curTime);

      // How will we advance?
      if (nonCrit) { // No critical reactions fire. Advance only noncritical reactions.

        // Compute the total propensity function in case it is faster to run SSA steps
        Real A = this->propensity(a_particles);

        if (A * curDt < m_SSAlim) { // Tau-leaping is inefficient, do some SSA on the WHOLE reaction set.
          Real dtSSA = 0.;
          int  nSSA  = 0;

          while (dtSSA < curDt && nSSA < m_NSSA) {

            // Update reaction rates, propensities, and the total propensity.
            const std::deque<Real> propensities = this->propensities(a_particles);
            A                                   = this->propensity(a_particles);

            // Draw time to next reaction.
            const Real u1     = m_udist01(m_rng);
            const Real dtColl = log(1. / u1) / A;
            const Real dtNext = Min(dtColl, curDt - dtSSA); // Don't exceed curDt with the next time step.

            if (dtSSA + dtColl < curDt) { // Collision happened inside curDt
              this->stepSSA(a_particles,
                            a_newPhotons,
                            a_mean_energies,
                            a_sources,
                            propensities,
                            A,
                            a_dx,
                            a_kappa,
                            dtColl);
            }
            else { // Next reaction is outside curDt, but we should still increment with the energy.
              this->stepEnergy(a_mean_energies, a_particles, a_sources, a_dx, a_kappa, dtNext);
            }

            // Update reaction rates and propensity for the next SSA step.
            this->updateReactionRates(a_E, a_mean_energies, a_dx, a_kappa);

            dtSSA += dtNext;
            nSSA += 1;
          }

          valid = true;
          curTime += dtSSA;
        }
        else { // Perform tau-leaping with non-critical reactions only

          // Need a copy in case step is rejected
          Vector<long long> particles     = a_particles;
          Vector<long long> Photons       = a_newPhotons;
          Vector<Real>      mean_energies = a_mean_energies;

          // Do the tau-leaping
          this->stepTau(particles, Photons, mean_energies, a_sources, nonCriticalReactions, a_dx, a_kappa, curDt);

          // Make sure step was valid
          valid = this->isStateValid(particles);
          if (valid) {
            curTime += curDt;
            a_particles     = particles;
            a_newPhotons    = Photons;
            a_mean_energies = mean_energies;
          }
          else {
            dtNonCrit = dtNonCrit * 0.5;
            continue;
          }
        }
      }
      else { // One critical reaction fired. Advance that reaction and use tau-leaping for the rest.
        // Need a copy in case step is rejected
        Vector<long long> particles     = a_particles;
        Vector<long long> Photons       = a_newPhotons;
        Vector<Real>      mean_energies = a_mean_energies;

        const std::deque<Real> propensitiesC  = this->propensities(particles, criticalReactions);
        const std::deque<Real> propensitiesNC = this->propensities(particles, nonCriticalReactions);

        const Real Ac = this->propensity(particles, criticalReactions);

        // SSA for critical reactions
        this->stepSSA(particles,
                      Photons,
                      mean_energies,
                      a_sources,
                      propensitiesC,
                      criticalReactions,
                      Ac,
                      a_dx,
                      a_kappa,
                      curDt);

        // Tau-leaping for non-critical reactions
        this->stepTau(particles,
                      Photons,
                      mean_energies,
                      a_sources,
                      nonCriticalReactions,
                      propensitiesNC,
                      a_dx,
                      a_kappa,
                      curDt);

        valid = this->isStateValid(particles);

        if (valid) {
          curTime += curDt;
          a_particles     = particles;
          a_newPhotons    = Photons;
          a_mean_energies = mean_energies;
        }
        else {
          dtNonCrit = dtNonCrit * 0.5;
          continue;
        }
      }
    }

    // Update reaction rates for next step
    this->updateReactionRates(a_E, a_mean_energies, a_dx, a_kappa);
  }
}

#include <CD_NamespaceFooter.H>

#endif
