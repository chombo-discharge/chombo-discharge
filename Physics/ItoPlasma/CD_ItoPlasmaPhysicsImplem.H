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
#include <CD_ParticleOps.H>
#include <CD_ParticleManagement.H>
#include <CD_ItoPlasmaPhysics.H>
#include <CD_Random.H>
#include <CD_DataParser.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

inline ItoPlasmaPhysics::ItoPlasmaPhysics()
{
  CH_TIME("ItoPlasmaPhysics::ItoPlasmaPhysics");

  m_coupling  = coupling::LFA;  

  m_className = "ItoPlasmaPhysics";

  m_reactions.clear();
  m_photoReactions.clear();

  // Some default settings in case user forgets to call the parsing algorithms.
  m_debug     = true;  
  m_ppc       = 32;
  m_Ncrit     = 5;
  m_eps       = 0.1;
  m_NSSA      = 10;
  m_SSAlim    = 5.0;
  m_algorithm = Algorithm::Hybrid;
}

inline ItoPlasmaPhysics::~ItoPlasmaPhysics() { CH_TIME("ItoPlasmaPhysics::~ItoPlasmaPhysics"); }

inline void
ItoPlasmaPhysics::defineKMC() noexcept
{
  CH_TIME("ItoPlasmaPhysics::defineKMC");

  m_kmcState.define(m_plasmaSpecies.size(), m_rtSpecies.size());

  m_kmcSolver.define(m_kmcReactions);

}

inline void
ItoPlasmaPhysics::parsePPC() noexcept
{
  CH_TIME("ItoPlasmaPhysics::parsePPC");

  ParmParse pp(m_className.c_str());

  pp.get("react_ppc", m_ppc);
}

inline void
ItoPlasmaPhysics::parseDebug() noexcept
{
  CH_TIME("ItoPlasmaPhysics::parseDebug");

  ParmParse pp(m_className.c_str());

  pp.get("debug", m_debug);
}

inline void
ItoPlasmaPhysics::parseAlgorithm() noexcept
{
  CH_TIME("ItoPlasmaPhysics::parseAlgorithm");

  ParmParse pp(m_className.c_str());

  std::string str;

  pp.get("algorithm", str);
  pp.get("Ncrit", m_Ncrit);
  pp.get("NSSA", m_NSSA);
  pp.get("prop_eps", m_eps);
  pp.get("SSAlim", m_SSAlim);

  if (str == "hybrid") {
    m_algorithm = Algorithm::Hybrid;
  }
  else if (str == "tau") {
    m_algorithm = Algorithm::Tau;
  }
  else if (str == "ssa") {
    m_algorithm = Algorithm::SSA;
  }
  else {
    MayDay::Abort("ItoPlasmaPhysics::parseAlgorithm - unknown algorithm requested");
  }
}

inline const Vector<RefCountedPtr<ItoSpecies>>&
ItoPlasmaPhysics::getItoSpecies() const
{
  return m_plasmaSpecies;
}

inline const Vector<RefCountedPtr<RtSpecies>>&
ItoPlasmaPhysics::getRtSpecies() const
{
  return m_rtSpecies;
}

inline int
ItoPlasmaPhysics::getNumPlasmaSpecies() const
{
  return m_plasmaSpecies.size();
}

inline int
ItoPlasmaPhysics::getNumPhotonSpecies() const
{
  return m_rtSpecies.size();
}

inline ItoPlasmaPhysics::coupling
ItoPlasmaPhysics::getCoupling() const
{
  return m_coupling;
}

inline Real
ItoPlasmaPhysics::initialSigma(const Real a_time, const RealVect a_pos) const
{
  return 0.0;
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
                                     const Real                  a_kappa) const noexcept
{
  CH_TIME("ItoPlasmaPhysics::reconcileParticles(Vector<List<ItoParticle>*>, ...)");

  for (int i = 0; i < a_particles.size(); i++) {
    const long long diff = a_newNumParticles[i] - a_oldNumParticles[i];

    if (diff > 0LL) {
      // Adding particles, which is fairly simple. Just choose weights for the particles and go.
      std::vector<long long> particleWeights = ParticleManagement::partitionParticleWeights(diff, (long long)m_ppc);

      for (const auto& w : particleWeights) {
        const RealVect x = Random::randomPosition(a_cellPos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);

        a_particles[i]->add(ItoParticle(w, x));
      }
    }
    else if (diff < 0LL) {
      // Removing particles is a bit more difficult because we need to manipulate weights.
      this->removeParticles(*a_particles[i], -diff);
    }
  }
}

inline void
ItoPlasmaPhysics::removeParticles(List<ItoParticle>& a_particles, const long long a_numParticlesToRemove) const
{
  CH_TIME("ItoPlasmaPhysics::removeParticles");

  constexpr long long zero = 0LL;

  CH_assert(a_numParticlesToRemove >= zero);

  // Quick lambda for getting total particle weight. Used for debugging.
  auto getTotalWeight = [&]() -> long long {
    long long W = zero;

    for (ListIterator<ItoParticle> lit(a_particles); lit.ok(); ++lit) {
      W += llround(lit().weight());

      if (lit().weight() < 1.0) {
        MayDay::Error("ItoPlasmaPhysics::removeParticles -- bad particle mass!");
      }
    }

    return W;
  };

  if (a_numParticlesToRemove > zero) {

    // For debugging only.
    long long totalWeightBefore = 0;
    long long totalWeightAfter  = 0;

    // Debug hook, compute the total particle weight before we start removing weights.
    if (m_debug) {
      totalWeightBefore = getTotalWeight();

      if (totalWeightBefore < a_numParticlesToRemove) {
        MayDay::Error("ItoPlasmaPhysics::removeParticles: logic bust (trying to remove too many particles)");
      }
    }

    // Remove physical particles.
    ParticleManagement::removePhysicalParticles(a_particles, a_numParticlesToRemove);

    // Remove particles with too low weight.
    ParticleManagement::deleteParticles(a_particles, std::numeric_limits<Real>::min());

    // Debug hook, make sure that particle weights are > 0 AND we've removed the desired
    // particle weight.
    if (m_debug) {
      totalWeightAfter = getTotalWeight();

      if (totalWeightBefore - totalWeightAfter != a_numParticlesToRemove) {
        pout() << "ItoPlasmaPhysics::removeParticles: Total weight before = " << totalWeightBefore << endl;
        pout() << "ItoPlasmaPhysics::removeParticles: Total weight after  = " << totalWeightAfter << endl;
        pout() << "ItoPlasmaPhysics::removeParticles: Should have removed = " << a_numParticlesToRemove << endl;

        MayDay::Abort("ItoPlasmaPhysics::removeParticles - incorrect mass removed");
      }
    }
  }
}

inline void
ItoPlasmaPhysics::reconcilePhotons(Vector<List<Photon>*>&   a_newPhotons,
                                   const Vector<long long>& a_numNewPhotons,
                                   const RealVect           a_cellPos,
                                   const RealVect           a_centroidPos,
                                   const RealVect           a_lo,
                                   const RealVect           a_hi,
                                   const RealVect           a_bndryCentroid,
                                   const RealVect           a_bndryNormal,
                                   const Real               a_dx,
                                   const Real               a_kappa) const noexcept
{
  CH_TIME("ItoPlasmaPhysics::reconcilePhotons");

  for (int i = 0; i < a_newPhotons.size(); i++) {
    if (a_numNewPhotons[i] > 0LL) {

      const std::vector<long long> photonWeights = ParticleManagement::partitionParticleWeights(a_numNewPhotons[i],
                                                                                                (long long)m_ppc);

      for (const auto& w : photonWeights) {
        const RealVect x = Random::randomPosition(a_cellPos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
        const RealVect v = Units::c * Random::getDirection();

        a_newPhotons[i]->add(Photon(x, v, m_rtSpecies[i]->getAbsorptionCoefficient(x), 1.0 * w));
      }
    }
  }
}

inline void
ItoPlasmaPhysics::reconcilePhotoionization(Vector<List<ItoParticle>*>&  a_particles,
                                           const Vector<List<Photon>*>& a_absorbedPhotons) const noexcept
{
  CH_TIME("ItoPlasmaPhysics::reconcilePhotoionization");

  for (const auto& r : m_photoReactions) {

    // Source and targets
    const int&             src     = r.second.source();
    const std::deque<int>& targets = r.second.targets();

    for (ListIterator<Photon> lit(*a_absorbedPhotons[src]); lit.ok(); ++lit) {
      const RealVect x = lit().position();
      const Real     w = lit().weight();

      for (const auto& t : targets) {
        a_particles[t]->add(ItoParticle(w, x));
      }
    }
  }
}

inline Vector<long long>
ItoPlasmaPhysics::getNumPhysicalParticles(const Vector<List<ItoParticle>*>& a_particles) const
{
  Vector<long long> ret(a_particles.size(), 0);

  for (int i = 0; i < a_particles.size(); i++) {
    for (ListIterator<ItoParticle> lit(*a_particles[i]); lit.ok(); ++lit) {
      ret[i] = llround(lit().weight());
    }
  }

  return ret;
}

inline void
ItoPlasmaPhysics::advanceReactionNetwork(Vector<List<ItoParticle>*>& a_particles,
                                         Vector<List<Photon>*>&      a_absorbedPhotons,
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
  CH_TIME("ItoPlasmaPhysics::advanceReactionNetwork");

  MayDay::Abort("ItoPlasmaPhysics::advanceReactionNetwork -- logic bust. Deprecated routine");

  // Get particle counts.
  Vector<long long> newPhotonCount   = Vector<long long>(m_numRtSpecies, 0);
  Vector<long long> oldParticleCount = this->getNumPhysicalParticles(a_particles);
  Vector<long long> newParticleCount = oldParticleCount;

  // Energies, mean_energies, and sources. Make sure LFA gives dummy arguments.
  Vector<Real> sources, energies, mean_energies;

  // Populate the KMC solver state.
  std::vector<long long>& reactives    = m_kmcState.getReactiveState();
  std::vector<long long>& nonReactives = m_kmcState.getNonReactiveState();

  for (size_t i = 0; i < a_particles.size(); i++) {
    reactives[i] = 0LL;
    for (ListIterator<ItoParticle> lit(*a_particles[i]); lit.ok(); ++lit) {
      reactives[i] += llround(lit().weight());
    }
  }

  // Run the Kinetic Monte Carlo solver.
  this->advanceParticles(newParticleCount, newPhotonCount, mean_energies, sources, a_dt, a_E, a_dx, a_kappa);

  // Reconcile particles -- this turns the new particle count into computational particles. This requires modifications
  // to the particle weights, removal of particles, and addition of particles.
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

  // Reconcile photons -- this creates new photons that the solver will use.
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

  // Run through photo-reactions and add photo-ionization products.
  this->reconcilePhotoionization(a_particles, a_absorbedPhotons);
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

    const Real u1 = Random::getUniformReal01();
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
  Vector<Real> energies(m_numPlasmaSpecies);
  for (int i = 0; i < a_particles.size(); i++) {
    energies[i] = Max(0.0, a_mean_energies[i] * a_particles[i]);
  }

  // Go through list of reactions.
  for (int i = 0; i < a_reactions.size(); i++) {
    const long long num_reactions = Random::getPoisson<long long>(a_propensities[i] * a_dt);

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
  this->updateReactionRates(a_E, a_dx, a_kappa);

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
        this->updateReactionRates(a_E, a_dx, a_kappa);
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
    const Real u2 = Random::getUniformReal01();

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
      this->updateReactionRates(a_E, a_dx, a_kappa);

      // Get all propensities and compute the total propensity.
      const std::deque<Real> propensities = this->propensities(a_particles);

      // Total propensity
      A = 0.0;
      for (const auto& a : propensities) {
        A += a;
      }

      // Random numebr and next collision time.
      u1     = Random::getUniformReal01();
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
  CH_TIME("ItoPlasmaPhysics::advanceParticles");

  // Define solver parameters in case they've been changed since last time. 
  m_kmcSolver.setSolverParameters(m_Ncrit, m_NSSA, m_eps, m_SSAlim);  

  // Populate the KMC solver state.
  std::vector<long long>& kmcParticles = m_kmcState.getReactiveState();
  std::vector<long long>& kmcPhotons   = m_kmcState.getNonReactiveState();

  for (size_t i = 0; i < a_particles.size(); i++) {
    kmcParticles[i] = a_particles[i];
  }

  for (auto& p : kmcPhotons) {
    p = 0LL;
  }

  // Update the reaction rates and run the KMC solver.
  this->updateReactionRates(a_E, a_dx, a_kappa);
  
  // // Which algorithm?
  if (m_algorithm == Algorithm::SSA) {
#if 0
    this->advanceSSA(a_particles, a_newPhotons, a_mean_energies, a_sources, a_dt, a_E, a_dx, a_kappa);
#else
    m_kmcSolver.advanceSSA(m_kmcState, a_dt);
#endif
  }
  else if (m_algorithm == Algorithm::Tau) {
#if 0
    this->advanceTau(a_particles, a_newPhotons, a_mean_energies, a_sources, a_dt, a_E, a_dx, a_kappa);
#else
    m_kmcSolver.advanceTau(m_kmcState, a_dt);
#endif
  }
  else if (m_algorithm == Algorithm::Hybrid) {
#if 0
    this->advanceHybrid(a_particles, a_newPhotons, a_mean_energies, a_sources, a_dt, a_E, a_dx, a_kappa);
#else
    m_kmcSolver.advanceHybrid(m_kmcState, a_dt);
#endif
  }

  // Put KMC back into ItoPlasma
  for (size_t i = 0; i < a_particles.size(); i++) {
    a_particles[i] = kmcParticles[i];
  }
  for (size_t i = 0; i < a_newPhotons.size(); i++) {
    a_newPhotons[i] = kmcPhotons[i];
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
    this->updateReactionRates(a_E, a_dx, a_kappa);

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
            const Real u1     = Random::getUniformReal01();
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
            this->updateReactionRates(a_E, a_dx, a_kappa);

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
    this->updateReactionRates(a_E, a_dx, a_kappa);
  }
}

#include <CD_NamespaceFooter.H>

#endif
