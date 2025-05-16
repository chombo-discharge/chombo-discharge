/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCPhysicsImplem.H 
  @brief  Implementation of CD_ItoKMCPhysics.H
  @author Robert Marskar
*/

#ifndef CD_ItoKMCPhysicsImplem_H
#define CD_ItoKMCPhysicsImplem_H

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_ItoKMCPhysics.H>
#include <CD_ParticleManagement.H>
#include <CD_Random.H>
#include <CD_Units.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

inline ItoKMCPhysics::ItoKMCPhysics() noexcept
{
  CH_TIME("ItoKMCPhysics::ItoKMCPhysics");

  m_className = "ItoKMCPhysics";

  m_kmcReactions.clear();
  m_photoReactions.clear();

  // Some default settings (mostly in case user forgets to call the parsing algorithms).
  m_isDefined         = false;
  m_debug             = true;
  m_hasKMCSolver      = false;
  m_downstreamSpecies = -1;
  m_maxNewParticles   = 32;
  m_maxNewPhotons     = 32;
  m_Ncrit             = 5;
  m_eps               = 2.0;
  m_NSSA              = 5;
  m_maxIter           = 10;
  m_SSAlim            = 5.0;
  m_exitTol           = 1.E-6;
  m_algorithm         = Algorithm::ExplicitEuler;
  m_particlePlacement = ParticlePlacement::Random;
}

inline ItoKMCPhysics::~ItoKMCPhysics() noexcept
{
  CH_TIME("ItoKMCPhysics::~ItoKMCPhysics");
}

inline void
ItoKMCPhysics::define() noexcept
{
  CH_TIME("ItoKMCPhysics::define");

  this->defineSpeciesMap();
  this->definePhotoPathways();

  // Safety hook -- make sure no one defines reactions using an out-of-range species index.
#ifndef NDEBUG
  for (const auto& R : m_kmcReactions) {
    const auto& lhsReactants = R.getReactants();
    const auto& rhsReactants = R.getReactiveProducts();
    const auto& rhsPhotons   = R.getNonReactiveProducts();

    for (const auto& r : lhsReactants) {
      CH_assert(r < m_itoSpecies.size() + m_cdrSpecies.size());
    }
    for (const auto& r : rhsReactants) {
      CH_assert(r < m_itoSpecies.size() + m_cdrSpecies.size());
    }
    for (const auto& r : rhsPhotons) {
      CH_assert(r < m_rtSpecies.size());
    }
  }
#endif

  m_isDefined = true;
}

inline void
ItoKMCPhysics::defineSpeciesMap() noexcept
{
  CH_TIME("ItoKMCPhysics::defineSpeciesMap");

  const int numItoSpecies = this->getNumItoSpecies();
  const int numCdrSpecies = this->getNumCdrSpecies();

  int species = 0;
  for (int i = 0; i < numItoSpecies; i++, species++) {
    m_speciesMap.emplace(species, std::make_pair(SpeciesType::Ito, i));
  }

  for (int i = 0; i < numCdrSpecies; i++, species++) {
    m_speciesMap.emplace(species, std::make_pair(SpeciesType::CDR, i));
  }
}

inline void
ItoKMCPhysics::defineKMC() const noexcept
{
  CH_TIME("ItoKMCPhysics::defineKMC");

  CH_assert(!m_hasKMCSolver);

  // Deep copy of reaction rates
  m_kmcReactionsThreadLocal.resize(0);
  for (const auto& r : m_kmcReactions) {
    m_kmcReactionsThreadLocal.emplace_back(std::make_shared<const KMCReaction>(r));
  }

  m_kmcSolver.define(m_kmcReactionsThreadLocal);
  m_kmcSolver.setSolverParameters(m_Ncrit, m_NSSA, m_maxIter, m_eps, m_SSAlim, m_exitTol);
  m_kmcState.define(m_itoSpecies.size() + m_cdrSpecies.size(), m_rtSpecies.size());

  m_hasKMCSolver = true;
}

inline void
ItoKMCPhysics::killKMC() const noexcept
{
  CH_TIME("ItoKMCPhysics::defineKMC");

  CH_assert(m_hasKMCSolver);

  m_kmcReactionsThreadLocal.resize(0);
  m_kmcSolver.define(m_kmcReactionsThreadLocal);
  m_kmcState.define(0, 0);

  m_hasKMCSolver = false;
}

inline void
ItoKMCPhysics::definePhotoPathways() noexcept
{
  CH_TIME("ItoKMCPhysics::definePhotoPathways");

  // Build a temporary list of pathways. I.e. restructure the list of reactions
  //
  // Y1 -> A
  // Y1 -> B
  // Y2 -> C
  // Y2 -> D
  //
  // into separate lists for Y1, Y2, ....
  //
  std::map<int, std::vector<std::pair<int, Real>>> pathways;

  for (int i = 0; i < m_photoReactions.size(); i++) {
    const ItoKMCPhotoReaction& r = m_photoReactions[i];

    const size_t& src        = r.getSourcePhoton();
    const Real    efficiency = r.getEfficiency();

    pathways[src].emplace_back(std::make_pair(i, efficiency));
  }

  // Go through the temporary pathways list and compute the relative efficiencies of one of the
  // photons triggering a reaction. The relative efficiencies are given by
  //
  // p(i) = R(i)/sum_j R(j).
  //
  for (const auto& p : pathways) {
    const int                               photoSpecies             = p.first;
    const std::vector<std::pair<int, Real>> reactionsAndEfficiencies = p.second;

    std::map<int, int> localToGlobalMap;
    std::list<double>  efficiencies;
    double             sumEfficiencies = 0.0;

    for (int i = 0; i < reactionsAndEfficiencies.size(); i++) {
      sumEfficiencies += (double)reactionsAndEfficiencies[i].second;
    }

    for (int i = 0; i < reactionsAndEfficiencies.size(); i++) {
      localToGlobalMap.emplace(i, reactionsAndEfficiencies[i].first);
      efficiencies.emplace_back((double)reactionsAndEfficiencies[i].second / sumEfficiencies);
    }

    std::discrete_distribution<int> distribution(efficiencies.begin(), efficiencies.end());

    m_photoPathways.insert(std::make_pair((int)photoSpecies, std::make_pair(distribution, localToGlobalMap)));
  }
}

inline const std::map<int, std::pair<SpeciesType, int>>&
ItoKMCPhysics::getSpeciesMap() const noexcept
{
  CH_TIME("ItoKMCPhysics::getSpeciesMap");

  return m_speciesMap;
}

inline void
ItoKMCPhysics::parseRuntimeOptions() noexcept
{
  CH_TIME("ItoKMCPhysics::parseRuntimeOptions");

  this->parsePPC();
  this->parseDebug();
  this->parseAlgorithm();
}

inline void
ItoKMCPhysics::parsePPC() noexcept
{
  CH_TIME("ItoKMCPhysics::parsePPC");

  ParmParse pp(m_className.c_str());

  pp.get("max_new_particles", m_maxNewParticles);
  pp.get("max_new_photons", m_maxNewPhotons);
}

inline void
ItoKMCPhysics::parseDebug() noexcept
{
  CH_TIME("ItoKMCPhysics::parseDebug");

  ParmParse pp(m_className.c_str());

  pp.get("debug", m_debug);
}

inline void
ItoKMCPhysics::parseAlgorithm() noexcept
{
  CH_TIME("ItoKMCPhysics::parseAlgorithm");

  ParmParse pp(m_className.c_str());

  std::string str;

  pp.get("algorithm", str);
  pp.get("crit_num", m_Ncrit);
  pp.get("SSA_num", m_NSSA);
  pp.get("prop_eps", m_eps);
  pp.get("SSA_lim", m_SSAlim);
  pp.get("max_iter", m_maxIter);
  pp.get("exit_tolerance", m_exitTol);

  if (str == "ssa") {
    m_algorithm = Algorithm::SSA;
  }
  else if (str == "explicit_euler") {
    m_algorithm = Algorithm::ExplicitEuler;
  }
  else if (str == "midpoint") {
    m_algorithm = Algorithm::Midpoint;
  }
  else if (str == "prc") {
    m_algorithm = Algorithm::PRC;
  }
  else if (str == "implicit_euler") {
    m_algorithm = Algorithm::ImplicitEuler;
  }
  else if (str == "hybrid_explicit_euler") {
    m_algorithm = Algorithm::HybridExplicitEuler;
  }
  else if (str == "hybrid_midpoint") {
    m_algorithm = Algorithm::HybridMidpoint;
  }
  else if (str == "hybrid_prc") {
    m_algorithm = Algorithm::HybridPRC;
  }
  else if (str == "hybrid_implicit_euler") {
    m_algorithm = Algorithm::HybridImplicitEuler;
  }
  else {
    MayDay::Error("ItoKMCPhysics::parseAlgorithm - unknown algorithm requested");
  }
}

inline const Vector<RefCountedPtr<ItoSpecies>>&
ItoKMCPhysics::getItoSpecies() const
{
  return m_itoSpecies;
}

inline const Vector<RefCountedPtr<CdrSpecies>>&
ItoKMCPhysics::getCdrSpecies() const
{
  return m_cdrSpecies;
}

inline const Vector<RefCountedPtr<RtSpecies>>&
ItoKMCPhysics::getRtSpecies() const
{
  return m_rtSpecies;
}

inline int
ItoKMCPhysics::getNumItoSpecies() const
{
  return m_itoSpecies.size();
}

inline int
ItoKMCPhysics::getNumCdrSpecies() const
{
  return m_cdrSpecies.size();
}

inline int
ItoKMCPhysics::getNumPlasmaSpecies() const
{
  return m_itoSpecies.size() + m_cdrSpecies.size();
}

inline int
ItoKMCPhysics::getNumPhotonSpecies() const
{
  return m_rtSpecies.size();
}

inline Real
ItoKMCPhysics::initialSigma(const Real a_time, const RealVect a_pos) const
{
  return 0.0;
}

inline void
ItoKMCPhysics::advanceKMC(Vector<FPR>&            a_numParticles,
                          Vector<FPR>&            a_numNewPhotons,
                          Real&                   a_criticalDt,
                          Real&                   a_nonCriticalDt,
                          const Vector<Real>&     a_phi,
                          const Vector<RealVect>& a_gradPhi,
                          const Real              a_dt,
                          const RealVect          a_E,
                          const RealVect          a_pos,
                          const Real              a_dx,
                          const Real              a_kappa) const
{
  CH_TIME("ItoKMCPhysics::advanceKMC");

  // Note: This is called PER GRID CELL, i.e. within OpenMP parallel regions. For this reason the KMC solver
  //       must be defined through defineKMC() (which must be later killed).
  CH_assert(m_isDefined);
  CH_assert(m_hasKMCSolver);

  std::vector<FPR>& kmcParticles = m_kmcState.getReactiveState();
  std::vector<FPR>& kmcPhotons   = m_kmcState.getNonReactiveState();

  for (size_t i = 0; i < a_numParticles.size(); i++) {
    kmcParticles[i] = a_numParticles[i];

    if(a_numParticles[i] < 0.0) {
      MayDay::Abort("negative particles on the way in");
    }
  }

  for (auto& p : kmcPhotons) {
    p = 0LL;
  }

  // Lambda function used for computing charge before and after reactions. Used only in debug mode
  // for ensuring that nothing goes wrong with charge conservation in the chemistry integration.
  auto computeCharge = [&]() -> Real {
    Real Q = 0.0;
    for (int i = 0; i < kmcParticles.size(); i++) {
      const SpeciesType& speciesType  = m_speciesMap.at(i).first;
      const int&         speciesIndex = m_speciesMap.at(i).second;

      int Z = 0;

      switch (speciesType) {
      case SpeciesType::Ito: {
        Z = m_itoSpecies[speciesIndex]->getChargeNumber();

        break;
      }
      case SpeciesType::CDR: {
        Z = m_cdrSpecies[speciesIndex]->getChargeNumber();

        break;
      }
      default: {
        MayDay::Abort("ItoKMCPhysics::advanceKMC -- logic bust in computeCharge()");

        break;
      }
      }

      Q += kmcParticles[i] * Z;
    }

    return Q;
  };

  // In debug mode, compute the total charge.
  const Real chargeBefore = m_debug ? computeCharge() : 0LL;

  // Update the reaction rates to be used by the KMC solver.
  this->updateReactionRates(m_kmcReactionsThreadLocal, a_E, a_pos, a_phi, a_gradPhi, a_dx, a_kappa);

  // Run the KMC solver.
  switch (m_algorithm) {
  case Algorithm::SSA: {
    m_kmcSolver.advanceSSA(m_kmcState, a_dt);

    break;
  }
  case Algorithm::ExplicitEuler: {
    m_kmcSolver.advanceTau(m_kmcState, a_dt, KMCLeapPropagator::ExplicitEuler);

    break;
  }
  case Algorithm::Midpoint: {
    m_kmcSolver.advanceTau(m_kmcState, a_dt, KMCLeapPropagator::Midpoint);

    break;
  }
  case Algorithm::PRC: {
    m_kmcSolver.advanceTau(m_kmcState, a_dt, KMCLeapPropagator::PRC);

    break;
  }
  case Algorithm::ImplicitEuler: {
    m_kmcSolver.advanceTau(m_kmcState, a_dt, KMCLeapPropagator::ImplicitEuler);

    break;
  }
  case Algorithm::HybridExplicitEuler: {
    m_kmcSolver.advanceHybrid(m_kmcState, a_dt, KMCLeapPropagator::ExplicitEuler);

    break;
  }
  case Algorithm::HybridMidpoint: {
    m_kmcSolver.advanceHybrid(m_kmcState, a_dt, KMCLeapPropagator::Midpoint);

    break;
  }
  case Algorithm::HybridPRC: {
    m_kmcSolver.advanceHybrid(m_kmcState, a_dt, KMCLeapPropagator::PRC);

    break;
  }
  case Algorithm::HybridImplicitEuler: {
    m_kmcSolver.advanceHybrid(m_kmcState, a_dt, KMCLeapPropagator::ImplicitEuler);

    break;
  }
  default: {
    MayDay::Error("ItoKMCPhysics::advanceKMC - logic bust");
  }
  }

  // Put KMC back into ItoKMC
  for (size_t i = 0; i < a_numParticles.size(); i++) {
    a_numParticles[i] = (FPR)kmcParticles[i];

    if(a_numParticles[i] < 0.0) {
      MayDay::Abort("negative particles on the way out");
    }    
  }
  for (size_t i = 0; i < a_numNewPhotons.size(); i++) {
    a_numNewPhotons[i] = (FPR)kmcPhotons[i];

    if(a_numNewPhotons[i] < 0.0) {
      MayDay::Abort("Cannot have < 0 photons");
    }
  }

  const Real chargeAfter = m_debug ? computeCharge() : 0LL;

  if (chargeAfter != chargeBefore) {

    std::cout << "charge before = " << chargeBefore << std::endl;
    std::cout << "charge after = " << chargeAfter << std::endl;
    std::cout << "diff = " << chargeAfter - chargeBefore << std::endl;
    std::cout << "rel diff = " << (chargeAfter/chargeBefore - 1.0) << std::endl;
    std::cout << kmcParticles[0] << std::endl;
    //    MayDay::Warning("ItoKMCPhysics::advanceKMC -- charge not conserved!");
  }

  // Compute the critical and non-critical time steps
  const auto partitionedReactions = m_kmcSolver.partitionReactions(m_kmcState);

  const auto& criticalReactions    = partitionedReactions.first;
  const auto& nonCriticalReactions = partitionedReactions.second;

  const auto criticalPropensities    = m_kmcSolver.propensities(m_kmcState, criticalReactions);
  const auto nonCriticalPropensities = m_kmcSolver.propensities(m_kmcState, nonCriticalReactions);

  a_criticalDt    = m_kmcSolver.getCriticalTimeStep(criticalPropensities);
  a_nonCriticalDt = m_kmcSolver.getNonCriticalTimeStep(m_kmcState, nonCriticalReactions, nonCriticalPropensities);
}

inline void
ItoKMCPhysics::reconcileParticles(Vector<List<ItoParticle>*>& a_particles,
                                  const Vector<FPR>&          a_newNumParticles,
                                  const Vector<FPR>&          a_oldNumParticles,
                                  const RealVect              a_electricField,
                                  const RealVect              a_cellPos,
                                  const RealVect              a_centroidPos,
                                  const RealVect              a_lo,
                                  const RealVect              a_hi,
                                  const RealVect              a_bndryCentroid,
                                  const RealVect              a_bndryNormal,
                                  const Real                  a_dx,
                                  const Real                  a_kappa) const noexcept
{
  CH_TIMERS("ItoKMCPhysics::reconcileParticles");
  CH_TIMER("ItoKMCPhysics::reconcileParticles::compute_downstream", t1);
  CH_TIMER("ItoKMCPhysics::reconcileParticles::particle_placement", t2);

  CH_assert(m_isDefined);
  CH_assert(a_particles.size() == a_newNumParticles.size());
  CH_assert(a_oldNumParticles.size() == a_newNumParticles.size());

  if (m_debug) {
    for (int i = 0; i < a_particles.size(); i++) {
      const FPR& numNew = a_newNumParticles[i];
      const FPR& numOld = a_oldNumParticles[i];

      if (numNew < (FPR)0.0) {
        MayDay::Warning("ItoKMCPhysics::reconcileParticles - new number of particles is < 0 (overflow issue?)");
      }
      else if ((long long)numNew < 0LL) {
        MayDay::Warning("ItoKMCPhysics::reconcileParticles - integer overflow!");
      }

      if (numOld < (FPR)0.0) {
        MayDay::Warning("ItoKMCPhysics::reconcileParticles - old number of particles is < 0");
      }
      else if ((long long)numOld < 0LL) {
        MayDay::Warning("ItoKMCPhysics::reconcileParticles - integer overflow for old particles!");
      }
    }
  }

  // Compute the upstream position of the particles (which is usually the electrons).
  bool     hasDownstream    = false;
  RealVect upstreamPosition = RealVect::Zero;
  RealVect upstreamLo       = -0.5 * RealVect::Unit;
  RealVect upstreamHi       = +0.5 * RealVect::Unit;
  RealVect v                = RealVect::Zero;
  int      Z                = 0;

  if (m_particlePlacement == ParticlePlacement::Downstream) {
    CH_assert(m_downstreamSpecies >= 0);

    Z = m_itoSpecies[m_downstreamSpecies]->getChargeNumber();

    if (Z != 0) {
      v = Z * a_electricField;
      v = v / v.vectorLength();
    }

    hasDownstream = this->computeUpstreamPosition(upstreamPosition,
                                                  upstreamLo,
                                                  upstreamHi,
                                                  Z,
                                                  *a_particles[m_downstreamSpecies],
                                                  a_electricField,
                                                  a_cellPos,
                                                  a_dx);
  }

  CH_START(t2);
  for (int i = 0; i < a_particles.size(); i++) {
    const long long diff = (long long)(a_newNumParticles[i] - a_oldNumParticles[i]);

    if (diff > 0LL) {
      // Adding particles, which is fairly simple. Just choose weights for the particles and go.
      const std::vector<long long> particleWeights = ParticleManagement::partitionParticleWeights(
        diff,
        (long long)m_maxNewParticles);

      for (const auto& w : particleWeights) {
        RealVect x = RealVect::Zero;

        // Figure out where to place the particles.
        switch (m_particlePlacement) {
        case ParticlePlacement::Centroid: {
          x = a_cellPos + a_centroidPos * a_dx;

          break;
        }
        case ParticlePlacement::Random: {
          x = Random::randomPosition(a_cellPos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);

          break;
        }
        case ParticlePlacement::Downstream: {
          if (hasDownstream) {
            x = Random::randomPosition(upstreamLo, upstreamHi, upstreamPosition, v);

            if ((x - a_bndryCentroid).dotProduct(a_bndryNormal) < 0.0) {
              x = upstreamPosition;
            }

            x = a_cellPos + x * a_dx;
          }
          else {

            x = Random::randomPosition(a_cellPos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
          }

          break;
        }
        default: {
          MayDay::Error("ItoKMCPhysics::reconcileParticles - logic bust");

          break;
        }
        }

        a_particles[i]->add(ItoParticle(1.0 * w, x));
      }
    }
    else if (diff < 0LL) {
      // Removing particles is a bit more difficult because we need to manipulate weights.
      this->removeParticles(*a_particles[i], -diff);
    }
  }
  CH_STOP(t2);
}

inline bool
ItoKMCPhysics::computeUpstreamPosition(RealVect&                a_pos,
                                       RealVect&                a_lo,
                                       RealVect&                a_hi,
                                       const int&               a_Z,
                                       const List<ItoParticle>& a_particles,
                                       const RealVect&          a_electricField,
                                       const RealVect&          a_cellPos,
                                       const Real&              a_dx) const noexcept
{
  CH_TIME("ItoKMCPhysics::computeUpstreamPosition");

  CH_assert(a_dx > 0.0);

  a_pos = RealVect::Zero;
  a_lo  = -0.5 * RealVect::Unit;
  a_hi  = +0.5 * RealVect::Unit;

  // Quick lambda for turning physical position into unit cell position.
  auto unitCellPosition = [&](const ItoParticle& p) -> RealVect {
    return (p.position() - a_cellPos) / a_dx;
  };

  const int      Z = (a_Z > 0) ? 1 : (a_Z < 0) ? -1 : 0;
  const RealVect E = a_electricField / a_electricField.vectorLength();
  const RealVect v = Z * E;

  // Upstream can only exist if there is a velocity direction and there are particles.
  const bool hasDownstream = (a_particles.length() > 0) && (v.vectorLength() > 0.0);

  if (hasDownstream) {

    Real D = std::numeric_limits<Real>::max();

    for (ListIterator<ItoParticle> lit(a_particles); lit.ok(); ++lit) {
      const RealVect x = unitCellPosition(lit());
      const Real     d = x.dotProduct(v);

      if (d < D) {
        D     = d;
        a_pos = x;
      }
    }

    DataOps::computeMinValidBox(a_lo, a_hi, v, a_pos);
  }

  if (m_debug) {
    for (int dir = 0; dir < SpaceDim; dir++) {
      if (a_pos[dir] > 0.5 || a_pos[dir] < -0.5) {
        MayDay::Abort("ItoKMCPhysics::computeUpstreamPosition - logic bust");
      }
    }
  }

  return hasDownstream;
}

inline void
ItoKMCPhysics::removeParticles(List<ItoParticle>& a_particles, const long long a_numParticlesToRemove) const
{
  CH_TIME("ItoKMCPhysics::removeParticles");

  constexpr long long zero = 0LL;

  CH_assert(m_isDefined);
  CH_assert(a_numParticlesToRemove >= zero);

  // Quick lambda for getting total particle weight. Used for debugging.
  auto getTotalWeight = [&]() -> long long {
    long long W = zero;

    for (ListIterator<ItoParticle> lit(a_particles); lit.ok(); ++lit) {
      W += llround(lit().weight());

      if (lit().weight() < 1.0) {
        MayDay::Error("ItoKMCPhysics::removeParticles -- bad particle mass!");
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
        MayDay::Error("ItoKMCPhysics::removeParticles: logic bust (trying to remove too many particles)");
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

      const long long errDiff = std::abs(totalWeightBefore - totalWeightAfter) - a_numParticlesToRemove;
      if (std::abs(errDiff) != zero) {

        pout() << "ItoKMCPhysics::removeParticles: Total weight before = " << totalWeightBefore << endl;
        pout() << "ItoKMCPhysics::removeParticles: Total weight after  = " << totalWeightAfter << endl;
        pout() << "ItoKMCPhysics::removeParticles: Should have removed = " << a_numParticlesToRemove << endl;
        pout() << "ItoKMCPhysics::removeParticles: Error               = " << errDiff << endl;

        MayDay::Abort("ItoKMCPhysics::removeParticles - incorrect mass removed");
      }
    }
  }
}

inline void
ItoKMCPhysics::reconcilePhotons(Vector<List<Photon>*>& a_newPhotons,
                                const Vector<FPR>&     a_numNewPhotons,
                                const RealVect         a_cellPos,
                                const RealVect         a_centroidPos,
                                const RealVect         a_lo,
                                const RealVect         a_hi,
                                const RealVect         a_bndryCentroid,
                                const RealVect         a_bndryNormal,
                                const Real             a_dx,
                                const Real             a_kappa) const noexcept
{
  CH_TIME("ItoKMCPhysics::reconcilePhotons");

  CH_assert(m_isDefined);

  for (int i = 0; i < a_newPhotons.size(); i++) {
    if (a_numNewPhotons[i] > 0LL) {

      const std::vector<long long> photonWeights = ParticleManagement::partitionParticleWeights(
        (long long)a_numNewPhotons[i],
        (long long)m_maxNewPhotons);

      for (const auto& w : photonWeights) {
        const RealVect x = Random::randomPosition(a_cellPos, a_lo, a_hi, a_bndryCentroid, a_bndryNormal, a_dx, a_kappa);
        const RealVect v = Units::c * Random::getDirection();

        a_newPhotons[i]->add(Photon(x, v, m_rtSpecies[i]->getAbsorptionCoefficient(x), 1.0 * w));
      }
    }
  }
}

inline void
ItoKMCPhysics::reconcilePhotoionization(Vector<List<ItoParticle>*>&   a_itoParticles,
                                        Vector<List<PointParticle>*>& a_cdrParticles,
                                        const Vector<List<Photon>*>&  a_absorbedPhotons) const noexcept
{
  CH_TIME("ItoKMCPhysics::reconcilePhotoionization");

  CH_assert(m_isDefined);
  CH_assert(a_itoParticles.size() == m_itoSpecies.size());
  CH_assert(a_cdrParticles.size() == m_cdrSpecies.size());

  for (int i = 0; i < a_absorbedPhotons.size(); i++) {
    if (m_photoPathways.find(i) != m_photoPathways.end()) {
      std::discrete_distribution<int> d                = m_photoPathways.at(i).first;
      const std::map<int, int>&       localToGlobalMap = m_photoPathways.at(i).second;

      for (ListIterator<Photon> lit(*a_absorbedPhotons[i]); lit.ok(); ++lit) {
        const RealVect x = lit().position();
        const Real     w = lit().weight();

        // Determine the photo-reaction type.
        const int localReaction  = Random::get(d);
        const int globalReaction = localToGlobalMap.at(localReaction);

        const ItoKMCPhotoReaction& photoReaction = m_photoReactions[globalReaction];
        const std::list<size_t>&   plasmaTargets = photoReaction.getTargetSpecies();

        for (const auto& t : plasmaTargets) {
          const SpeciesType& type       = m_speciesMap.at(t).first;
          const int&         localIndex = m_speciesMap.at(t).second;

          if (type == SpeciesType::Ito) {
            a_itoParticles[localIndex]->add(ItoParticle(w, x));
          }
          else if (type == SpeciesType::CDR) {
            a_cdrParticles[localIndex]->add(PointParticle(x, w));
          }
          else {
            MayDay::Error("CD_ItoKMCPhysics.H - logic bust in reconcilePhotoionization");
          }
        }
      }
    }
  }
}

#include <CD_NamespaceFooter.H>

#endif
