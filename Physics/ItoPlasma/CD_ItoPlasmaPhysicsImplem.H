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

// Our includes
#include <CD_ItoPlasmaPhysics.H>
#include <CD_ParticleManagement.H>
#include <CD_Random.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

inline ItoPlasmaPhysics::ItoPlasmaPhysics()
{
  CH_TIME("ItoPlasmaPhysics::ItoPlasmaPhysics");

  m_className = "ItoPlasmaPhysics";

  m_kmcReactions.clear();
  m_photoReactions.clear();

  // Some default settings in case user forgets to call the parsing algorithms.
  m_debug     = true;
  m_ppc       = 32LL;
  m_Ncrit     = 5;
  m_eps       = 0.1;
  m_NSSA      = 10;
  m_SSAlim    = 5.0;
  m_algorithm = Algorithm::Hybrid;
}

inline ItoPlasmaPhysics::~ItoPlasmaPhysics() { CH_TIME("ItoPlasmaPhysics::~ItoPlasmaPhysics"); }

inline Real
ItoPlasmaPhysics::computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_densities) const noexcept
{
  CH_TIME("ItoPlasmaPhysics::computeDt");

  return std::numeric_limits<Real>::max();
}

inline void
ItoPlasmaPhysics::defineKMC() noexcept
{
  CH_TIME("ItoPlasmaPhysics::defineKMC");

  m_kmcState.define(m_plasmaSpecies.size(), m_rtSpecies.size());

  m_kmcSolver.define(m_kmcReactions);
}

inline void
ItoPlasmaPhysics::parseRuntimeOptions() noexcept
{
  CH_TIME("ItoPlasmaPhysics::parseRuntimeOptions");

  this->parsePPC();
  this->parseDebug();
  this->parseAlgorithm();
}

inline void
ItoPlasmaPhysics::parsePPC() noexcept
{
  CH_TIME("ItoPlasmaPhysics::parsePPC");

  ParmParse pp(m_className.c_str());

  int ppc;

  pp.get("react_ppc", ppc);

  m_ppc = (long long)ppc;
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
    MayDay::Error("ItoPlasmaPhysics::parseAlgorithm - unknown algorithm requested");
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

inline Real
ItoPlasmaPhysics::initialSigma(const Real a_time, const RealVect a_pos) const
{
  return 0.0;
}

inline void
ItoPlasmaPhysics::advanceKMC(Vector<long long>& a_numParticles,
                             Vector<long long>& a_numNewPhotons,
                             const Real         a_dt,
                             const RealVect     a_E,
                             const Real         a_dx,
                             const Real         a_kappa) const
{
  CH_TIME("ItoPlasmaPhysics::advanceKMC");

  // Define solver parameters in case they've been changed since last time.
  m_kmcSolver.setSolverParameters(m_Ncrit, m_NSSA, m_eps, m_SSAlim);

  // Populate the KMC solver state.
  std::vector<long long>& kmcParticles = m_kmcState.getReactiveState();
  std::vector<long long>& kmcPhotons   = m_kmcState.getNonReactiveState();

  for (size_t i = 0; i < a_numParticles.size(); i++) {
    kmcParticles[i] = a_numParticles[i];
  }

  for (auto& p : kmcPhotons) {
    p = 0LL;
  }

  // Update the reaction rates and run the KMC solver.
  this->updateReactionRates(a_E, a_dx, a_kappa);

  // // Which algorithm?
  if (m_algorithm == Algorithm::SSA) {
    m_kmcSolver.advanceSSA(m_kmcState, a_dt);
  }
  else if (m_algorithm == Algorithm::Tau) {
    m_kmcSolver.advanceTau(m_kmcState, a_dt);
  }
  else if (m_algorithm == Algorithm::Hybrid) {
    m_kmcSolver.advanceHybrid(m_kmcState, a_dt);
  }

  // Put KMC back into ItoPlasma
  for (size_t i = 0; i < a_numParticles.size(); i++) {
    a_numParticles[i] = kmcParticles[i];
  }
  for (size_t i = 0; i < a_numNewPhotons.size(); i++) {
    a_numNewPhotons[i] = kmcPhotons[i];
  }
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

    if (m_debug && a_newNumParticles[i] < 0LL) {
      MayDay::Abort("ItoPlasmaPhysics::reconcileParticles - crap");
    }

    if (diff > 0LL) {
      // Adding particles, which is fairly simple. Just choose weights for the particles and go.
      std::vector<long long> particleWeights = ParticleManagement::partitionParticleWeights(diff, m_ppc);

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
                                                                                                m_ppc);

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
    const size_t&            src     = r->getSourcePhoton();
    const std::list<size_t>& targets = r->getTargetSpecies();

    for (ListIterator<Photon> lit(*a_absorbedPhotons[src]); lit.ok(); ++lit) {
      const RealVect x = lit().position();
      const Real     w = lit().weight();

      for (const auto& t : targets) {
        a_particles[t]->add(ItoParticle(w, x));
      }
    }
  }
}

#include <CD_NamespaceFooter.H>

#endif
