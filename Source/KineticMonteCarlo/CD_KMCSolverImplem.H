/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCSolverImplem.H
  @brief  Implementation of CD_KMCSolver.H
  @author Robert Marskar
*/

#ifndef CD_KMCSolverImplem_H
#define CD_KMCSolverImplem_H

// Std includes
#include <limits>
#include <unordered_set>

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_Random.H>
#include <CD_KMCSolver.H>
#include <CD_LaPackUtils.H>
#include <CD_NamespaceHeader.H>

template <typename R, typename State, typename T>
inline KMCSolver<R, State, T>::KMCSolver() noexcept
{
  CH_TIME("KMCSolver::KMCSolver");

  this->setSolverParameters(0, 0, 100, std::numeric_limits<Real>::max(), 0.0, 1.E-6);
}

template <typename R, typename State, typename T>
inline KMCSolver<R, State, T>::KMCSolver(const ReactionList& a_reactions) noexcept
{
  CH_TIME("KMCSolver::KMCSolver");

  this->define(a_reactions);
}

template <typename R, typename State, typename T>
inline KMCSolver<R, State, T>::~KMCSolver() noexcept
{
  CH_TIME("KMCSolver::~KMCSolver");
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::define(const ReactionList& a_reactions) noexcept
{
  CH_TIME("KMCSolver::define");

  m_reactions = a_reactions;

  // Default settings. These are equivalent to ALWAYS using tau-leaping.
  this->setSolverParameters(0, 0, 100, std::numeric_limits<Real>::max(), 0.0, 1.E-6);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::setSolverParameters(const T    a_Ncrit,
                                            const T    a_numSSA,
                                            const T    a_maxIter,
                                            const Real a_eps,
                                            const Real a_SSAlim,
                                            const Real a_exitTol) noexcept
{
  CH_TIME("KMCSolver::setSolverParameters");

  m_Ncrit   = a_Ncrit;
  m_numSSA  = a_numSSA;
  m_maxIter = a_maxIter;
  m_eps     = a_eps;
  m_SSAlim  = a_SSAlim;
  m_exitTol = a_exitTol;
}

template <typename R, typename State, typename T>
inline std::vector<std::vector<T>>
KMCSolver<R, State, T>::getNu(const State& a_state, const ReactionList& a_reactions) const noexcept
{
  CH_TIME("KMCSolver::getNu(State)");

  std::vector<std::vector<T>> ret(a_reactions.size());

  for (int j = 0; j < a_reactions.size(); j++) {
    std::vector<T>& nu = ret[j];

    // Linearize a trivial state
    State state = a_state;

    std::vector<T> preState = state.linearOut();
    for (auto& x : preState) {
      x = static_cast<T>(0);
    }
    state.linearIn(preState);

    // Advance with exactly one reaction.
    a_reactions[j]->advanceState(state, static_cast<T>(1));

    std::vector<T> postState = state.linearOut();

    // Compute the state change vector.
    nu.resize(preState.size());
    for (int i = 0; i < preState.size(); i++) {
      nu[i] = postState[i] - preState[i];
    }
  }

  return ret;
}

template <typename R, typename State, typename T>
inline std::vector<Real>
KMCSolver<R, State, T>::propensities(const State& a_state) const noexcept
{
  CH_TIME("KMCSolver::propensities(State)");

  return this->propensities(a_state, m_reactions);
}

template <typename R, typename State, typename T>
inline std::vector<Real>
KMCSolver<R, State, T>::propensities(const State& a_state, const ReactionList& a_reactions) const noexcept
{
  CH_TIME("KMCSolver::propensities(State, ReactionList)");

  std::vector<Real> A(a_reactions.size());

  const size_t numReactions = a_reactions.size();

  for (size_t i = 0; i < numReactions; i++) {
    A[i] = a_reactions[i]->propensity(a_state);
  }

  return A;
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::totalPropensity(const State& a_state) const noexcept
{
  CH_TIME("KMCSolver::totalPropensity(State)");

  return this->totalPropensity(a_state, m_reactions);
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::totalPropensity(const State& a_state, const ReactionList& a_reactions) const noexcept
{
  CH_TIME("KMCSolver::totalPropensity(State, ReactionList)");

  Real A = 0.0;

  const size_t numReactions = a_reactions.size();

  for (size_t i = 0; i < numReactions; i++) {
    A += a_reactions[i]->propensity(a_state);
  }

  return A;
}

template <typename R, typename State, typename T>
inline std::pair<typename KMCSolver<R, State, T>::ReactionList, typename KMCSolver<R, State, T>::ReactionList>
KMCSolver<R, State, T>::partitionReactions(const State& a_state) const noexcept
{
  CH_TIME("KMCSolver::partitionReactions(State)");

  return this->partitionReactions(a_state, m_reactions);
}

template <typename R, typename State, typename T>
inline std::pair<typename KMCSolver<R, State, T>::ReactionList, typename KMCSolver<R, State, T>::ReactionList>
KMCSolver<R, State, T>::partitionReactions(const State& a_state, const ReactionList& a_reactions) const noexcept
{
  CH_TIME("KMCSolver::partitionReactions(State, ReactionList)");

  ReactionList criticalReactions;
  ReactionList nonCriticalReactions;

  const size_t numReactions = a_reactions.size();

  for (size_t i = 0; i < numReactions; i++) {
    const T Lj = a_reactions[i]->computeCriticalNumberOfReactions(a_state);

    if (Lj < m_Ncrit) {
      criticalReactions.emplace_back(a_reactions[i]);
    }
    else {
      nonCriticalReactions.emplace_back(a_reactions[i]);
    }
  }

  return std::make_pair(criticalReactions, nonCriticalReactions);
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::getCriticalTimeStep(const State& a_state) const noexcept

{
  return this->getCriticalTimeStep(a_state, m_reactions);
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::getCriticalTimeStep(const State&        a_state,
                                            const ReactionList& a_criticalReactions) const noexcept
{
  CH_TIME("KMCSolver::getCriticalTimeStep");

  // TLDR: This computes the time until the firing of the next critical reaction.

  Real dt = std::numeric_limits<Real>::max();

  if (a_criticalReactions.size() > 0) {
    // Add numeric_limits<Real>::min to A and u to avoid division by zero.
    const Real A = std::numeric_limits<Real>::min() + this->totalPropensity(a_state, a_criticalReactions);

    dt = this->getCriticalTimeStep(A);
  }

  return dt;
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::getCriticalTimeStep(const std::vector<Real>& a_propensities) const noexcept
{
  CH_TIME("KineticMonteCarlo::getCriticalTimeStep(std::vector<Real>)");

  // To avoid division by zero later on.
  Real A = std::numeric_limits<Real>::min();

  for (const auto& p : a_propensities) {
    A += p;
  }

  return this->getCriticalTimeStep(A);
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::getCriticalTimeStep(const Real& a_totalPropensity) const noexcept
{
  CH_TIME("KineticMonteCarlo::getCriticalTimeStep(Real)");

  const Real u = std::numeric_limits<Real>::min() + Random::getUniformReal01();

  return log(1.0 / u) / a_totalPropensity;
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::getNonCriticalTimeStep(const State& a_state) const noexcept
{
  CH_TIME("KMCSolver::getNonCriticalTimeStep(State");

  const auto& partitionedReactions = this->partitionReactions(a_state, m_reactions);

  return this->getNonCriticalTimeStep(a_state, partitionedReactions.second);
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::getNonCriticalTimeStep(const State& a_state, const ReactionList& a_reactions) const noexcept
{
  CH_TIME("KMCSolver::getNonCriticalTimeStep(State, ReactionList");

  const std::vector<Real> propensities = this->propensities(a_state, a_reactions);

  return this->getNonCriticalTimeStep(a_state, a_reactions, propensities);
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::getNonCriticalTimeStep(const State&             a_state,
                                               const ReactionList&      a_nonCriticalReactions,
                                               const std::vector<Real>& a_nonCriticalPropensities) const noexcept
{
  CH_TIME("KMCSolver::getNonCriticalTimeStep(State, ReactionList, std::vector<Real>");

  CH_assert(a_nonCriticalReactions.size() == a_nonCriticalPropensities.size());

  constexpr Real one = 1.0;

  Real dt = std::numeric_limits<Real>::max();

  const size_t numReactions = a_nonCriticalReactions.size();

  if (numReactions > 0) {

    // 1. Gather a list of all reactants involved in the non-critical reactions.
    auto allReactants = a_nonCriticalReactions[0]->getReactants();
    for (size_t i = 1; i < numReactions; i++) {
      const auto& curReactants = a_nonCriticalReactions[i]->getReactants();

      allReactants.insert(allReactants.end(), curReactants.begin(), curReactants.end());
    }

    // Only do unique reactants.
    const auto ip = std::unique(allReactants.begin(), allReactants.end());
    allReactants.resize(std::distance(allReactants.begin(), ip));

    // 2. Iterate through all reactants and compute deviations.
    for (const auto& reactant : allReactants) {

      // Xi is the population of the current reactant. It might seem weird that we are indexing this
      // through the reactions rather then through the state. Which it is, but the reason for that is
      // that it is the REACTION that determines how we index the state population. This is simply a
      // design choice that permits the user to apply different type of reactions without changing
      // the underlying state.
      const T Xi = a_nonCriticalReactions[0]->population(reactant, a_state);

      if (Xi > (T)0) {
        Real mu     = 0.0;
        Real sigma2 = 0.0;

        // Set gi to 1 for now. A more complex version would parse this through an input parameter
        // where the user has inspected the highest-order-reaction.
        constexpr Real gi = 1.0;

        for (size_t i = 0; i < numReactions; i++) {
          const Real& p    = a_nonCriticalPropensities[i];
          const auto& muIJ = a_nonCriticalReactions[i]->getStateChange(reactant);

          mu += muIJ * p;
          sigma2 += muIJ * muIJ * p;
        }

        Real dt1 = std::numeric_limits<Real>::max();
        Real dt2 = std::numeric_limits<Real>::max();

        const Real f = std::max(m_eps * Xi / gi, one);

        mu     = std::abs(mu);
        sigma2 = std::abs(sigma2);

        if (mu > std::numeric_limits<Real>::min()) {
          dt1 = f / mu;
        }
        if (sigma2 > std::numeric_limits<Real>::min()) {
          dt2 = (f * f) / sigma2;
        }

        dt = std::min(dt, std::min(dt1, dt2));
      }
    }
  }

  return dt;
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::computeDt(const State&             a_state,
                                  const ReactionList&      a_nonCriticalReactions,
                                  const std::vector<Real>& a_nonCriticalPropensities,
                                  const Real               a_epsilon) const noexcept
{
  CH_TIME("KMCSolver::getNonCriticalTimeStep(State, ReactionList, std::vector<Real>, Real");

  CH_assert(a_nonCriticalReactions.size() == a_nonCriticalPropensities.size());

  constexpr Real one = 1.0;

  Real dt = std::numeric_limits<Real>::max();

  const size_t numReactions = a_nonCriticalReactions.size();

  if (numReactions > 0) {

    // 1. Gather a list of all reactants involved in the non-critical reactions.
    auto allReactants = a_nonCriticalReactions[0]->getReactants();
    for (size_t i = 1; i < numReactions; i++) {
      const auto& curReactants = a_nonCriticalReactions[i]->getReactants();

      allReactants.insert(allReactants.end(), curReactants.begin(), curReactants.end());
    }

    // Only do unique reactants.
    const auto ip = std::unique(allReactants.begin(), allReactants.end());
    allReactants.resize(std::distance(allReactants.begin(), ip));

    // 2. Iterate through all reactants and compute deviations.
    for (const auto& reactant : allReactants) {

      // Xi is the population of the current reactant. It might seem weird that we are indexing this
      // through the reactions rather then through the state. Which it is, but the reason for that is
      // that it is the REACTION that determines how we index the state population. This is simply a
      // design choice that permits the user to apply different type of reactions without changing
      // the underlying state.
      const T Xi = a_nonCriticalReactions[0]->population(reactant, a_state);

      if (Xi > (T)0) {
        Real mu = 0.0;

        // Set gi to 1 for now. A more complex version would parse this through an input parameter
        // where the user has inspected the highest-order-reaction.
        constexpr Real gi = 1.0;

        for (size_t i = 0; i < numReactions; i++) {
          const Real& p    = a_nonCriticalPropensities[i];
          const auto& muIJ = a_nonCriticalReactions[i]->getStateChange(reactant);

          mu += muIJ * p;
        }

        const Real f = std::max(a_epsilon * Xi / gi, one);

        mu = std::abs(mu);

        if (mu > std::numeric_limits<Real>::min()) {
          dt = std::min(dt, f / mu);
        }
      }
    }
  }

  return dt;
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepSSA(State& a_state) const noexcept
{
  CH_TIME("KMCSolver::stepSSA(State, Real)");

  this->stepSSA(a_state, m_reactions);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepSSA(State& a_state, const ReactionList& a_reactions) const noexcept
{
  CH_TIME("KMCSolver::stepSSA(State, ReactionList, Real)");

  if (a_reactions.size() > 0) {

    // Compute all propensities.
    const std::vector<Real> propensities = this->propensities(a_state, a_reactions);

    this->stepSSA(a_state, a_reactions, propensities);
  }
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepSSA(State&                   a_state,
                                const ReactionList&      a_reactions,
                                const std::vector<Real>& a_propensities) const noexcept
{
  CH_TIME("KMCSolver::stepSSA(State, ReactionList, std::vector<Real>)");

  CH_assert(a_reactions.size() == a_propensities.size());

  const size_t numReactions = a_reactions.size();

  if (numReactions > 0) {
    constexpr T one = (T)1;

    // Determine the reaction type as per Gillespie algorithm.
    Real A = 0.0;
    for (size_t i = 0; i < numReactions; i++) {
      A += a_propensities[i];
    }

    const Real u = Random::getUniformReal01();

    size_t r = std::numeric_limits<size_t>::max();

    Real sumProp = 0.0;
    for (size_t i = 0; i < a_reactions.size(); i++) {
      sumProp += a_propensities[i];

      if (sumProp >= u * A) {
        r = i;

        break;
      }
    }

    // Advance by one reaction.
    a_reactions[r]->advanceState(a_state, one);
  }
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::advanceSSA(State& a_state, const Real a_dt) const noexcept
{
  CH_TIME("KMCSolver::advanceSSA(State, Real)");

  this->advanceSSA(a_state, m_reactions, a_dt);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::advanceSSA(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept
{
  CH_TIME("KMCSolver::advanceSSA(State, ReactionList, Real)");

  const size_t numReactions = a_reactions.size();

  if (numReactions > 0) {

    // Simulated time within the SSA.
    Real curDt = 0.0;

    while (curDt <= a_dt) {

      // Compute the propensities and get the time to the next reaction.
      const std::vector<Real> propensities = this->propensities(a_state, a_reactions);

      const Real nextDt = this->getCriticalTimeStep(propensities);

      // Fire one reaction if occurs within a_dt.
      if (curDt + nextDt <= a_dt) {
        this->stepSSA(a_state, a_reactions, propensities);
      }

      curDt += nextDt;
    }
  }
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepExplicitEuler(State& a_state, const Real a_dt) const noexcept
{
  CH_TIME("KMCSolver::stepExplicitEuler(State, Real)");

  this->stepExplicitEuler(a_state, m_reactions, a_dt);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepExplicitEuler(State&              a_state,
                                          const ReactionList& a_reactions,
                                          const Real          a_dt) const noexcept
{
  CH_TIME("KMCSolver::stepExplicitEuler(State, ReactionList, std::vector<Real>, Real)");

  CH_assert(a_dt > 0.0);

  if (a_reactions.size() > 0) {
    const std::vector<Real> propensities = this->propensities(a_state, a_reactions);

    for (size_t i = 0; i < a_reactions.size(); i++) {

      // Number of reactions is always an integer -- draw from a Poisson distribution in long long. I'm just
      // using a large integer type to avoid potential overflows.
      const T numReactions = (T)Random::getPoisson<long long>(propensities[i] * a_dt);

      a_reactions[i]->advanceState(a_state, numReactions);
    }
  }
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepMidpoint(State& a_state, const Real a_dt) const noexcept
{
  CH_TIME("KMCSolver::stepMidpoint(State&, const Real)");

  this->stepMidpoint(a_state, m_reactions, a_dt);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepMidpoint(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept
{
  CH_TIMERS("KMCSolver::stepMidpoint(State&, const ReactionList&, const Real)");
  CH_TIMER("KMCSolver::stepTauMidPoint::first_loop", t1);
  CH_TIMER("KMCSolver::stepTauMidPoint::second_loop", t2);

  const int numReactions = a_reactions.size();

  if (numReactions > 0) {

    std::vector<Real> propensities = this->propensities(a_state, a_reactions);

    State Xdagger = a_state;

    CH_START(t1);
    for (size_t i = 0; i < numReactions; i++) {
      // TLDR: Predict a midpoint state -- unfortunately this means that as a_dt->0 we end up with plain
      //       tau-leaping. I don't know of a way to fix this without introducing double fluctuations (yet).
      a_reactions[i]->advanceState(Xdagger, (T)std::round(0.5 * propensities[i] * a_dt));
    }
    CH_STOP(t1);

    propensities = this->propensities(Xdagger, a_reactions);

    CH_START(t2);
    for (size_t i = 0; i < numReactions; i++) {
      const T curReactions = (T)Random::getPoisson<long long>(propensities[i] * a_dt);

      a_reactions[i]->advanceState(a_state, curReactions);
    }
    CH_STOP(t2);
  }
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepPRC(State& a_state, const Real a_dt) const noexcept
{
  CH_TIME("KMCSolver::stepPRC(State&, const Real)");

  this->stepPRC(a_state, m_reactions, a_dt);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepPRC(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept
{
  CH_TIMERS("KMCSolver::stepPRC(State&, const ReactionList&, const Real)");
  CH_TIMER("KMCSolver::stepPRC::first_loop", t1);
  CH_TIMER("KMCSolver::stepPRC::second_loop", t2);

  const int numReactions = a_reactions.size();

  if (numReactions > 0) {

    std::vector<Real> aj = this->propensities(a_state, a_reactions);

    const std::vector<Real> ak = aj;

    CH_START(t1);
    for (int j = 0; j < numReactions; j++) {
      for (int k = 0; k < numReactions; k++) {
        State x = a_state;

        a_reactions[k]->advanceState(x, (T)1);

        const Real etajk = a_reactions[j]->propensity(x) - ak[j];

        aj[j] += 0.5 * a_dt * ak[k] * etajk;
      }
    }
    CH_STOP(t1);

    CH_START(t2);
    for (size_t i = 0; i < numReactions; i++) {
      const T nr = (T)Random::getPoisson<long long>(aj[i] * a_dt);

      a_reactions[i]->advanceState(a_state, nr);
    }
    CH_STOP(t2);
  }
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepImplicitEuler(State& a_state, const Real a_dt) const noexcept
{
  CH_TIME("KMCSolver::stepImplicitEuler(State, a_dt)");

  this->stepImplicitEuler(a_state, m_reactions, a_dt);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepImplicitEuler(State&              a_state,
                                          const ReactionList& a_reactions,
                                          const Real          a_dt) const noexcept
{
  CH_TIME("KMCSolver::stepImplicitEuler(State, ReactionList, a_dt)");

  // TLDR: The implicit Euler tau-leaping scheme is equivalent to the solution of
  //
  // F(X) = X - (x + sum_j (nu_j * (P(a(x)*dt) - a(x)*dt)]) - dt*sum_j nu_j a_j(X)
  //      = X - c - dt*sum_j nu_j * a_j(X),
  //      = 0
  //
  // where c = (x + sum_j (nu_j * (P(a(x)*dt) - a(x)*dt)]) is a constant term throughout the Newton iterations. This term is presampled
  // before the Newton iterations begin.

  // Linearize input state onto something understandable to LAPACK
  const std::vector<T>              inputState = a_state.linearOut();
  const std::vector<std::vector<T>> nu         = this->getNu(a_state, a_reactions);
  const std::vector<Real>           ajX        = this->propensities(a_state, a_reactions);

  // Number of equations and number of reactions.
  const int N = inputState.size();
  const int M = a_reactions.size();

  // Lambda that computes the constant term c and the initial guess.
  auto compConstantTerm = [&](double* C, double* X, const State& state, const Real a_dt) -> void {
    // TLDR: To compute the constant term we perform a Poisson sampling and then linearize the output
    // state.

    State explicitEulerState = state;

    this->advanceTau(explicitEulerState, a_reactions, a_dt, KMCLeapPropagator::ExplicitEuler);

    const std::vector<T> eulerOut = explicitEulerState.linearOut();

    // Make c = (x + sum_j (nu_j * (P(a(x)*dt) - a(x)*dt)])
    for (int i = 0; i < N; i++) {
      C[i] = 1.0 * eulerOut[i];
      X[i] = 1.0 * eulerOut[i];
    }

    // Subtract the mean.
    for (int j = 0; j < M; j++) {
      const std::vector<T>& nuj = nu[j];

      for (int i = 0; i < N; i++) {
        C[i] -= nuj[i] * ajX[j] * a_dt;
      }
    }
  };

  // Lambda that computes the each equation.
  auto computeF = [&](double* F, const double* Xit, const double* C) -> void {
    // Compute the propensities for the Xit state. This requires us to round to the nearest integer.
    State stateXit = a_state;

    std::vector<T> linState(N);
    for (int i = 0; i < N; i++) {
      linState[i] = static_cast<T>(llround(Xit[i]));
    }

    stateXit.linearIn(linState);

    const std::vector<Real> ajXit = this->propensities(stateXit);

    // Compute Xit - c - dt * sum_j nu_j * aj(round(Xit))
    for (int i = 0; i < N; i++) {
      F[i] = Xit[i] - C[i];
    }

    for (int j = 0; j < M; j++) {
      const std::vector<T>& nuj = nu[j];

      for (int i = 0; i < N; i++) {
        F[i] -= a_dt * nuj[i] * ajXit[j];
      }
    }
  };

  // Compute the max-norm
  auto computeNorm = [&](double* F, const double* X, const double* C) -> Real {
    computeF(F, X, C);

    Real norm = 0.0;

    for (int i = 0; i < N; i++) {
      norm = std::max(norm, std::abs(F[i]));
    }

    return norm;
  };

  // Allocate memory for the Jacobian (J), the Newton increment (X), and F(X) (F). Also include the constant term from the Poisson sampling (c)
  double* J = new double[N * N];
  double* X = new double[N];
  double* F = new double[N];
  double* C = new double[N];

  // Temporary storage used for computing the Jacobian matrix.
  double* X1 = new double[N];
  double* X2 = new double[N];
  double* F1 = new double[N];
  double* F2 = new double[N];

  // Things that are required by LAPACK
  int  INFO = 0;
  int  NRHS = 1;
  int* IPIV = new int[N];

  // Compute the constant term.
  compConstantTerm(C, X, a_state, a_dt);

  // Compute the max-norm of F(0) to use as an exit criterion.
  for (int i = 0; i < N; i++) {
    X1[i] = 0.0;
  }

  const Real initNorm = computeNorm(F, X1, C);

  bool converged = true;

  for (int k = 0; k < m_maxIter; k++) {

    // Numerically computed the Jacobian using finite differences. The Jacobian is given by
    // Jij = dF_i/dx_j
    for (int j = 0; j < N; j++) {

      for (int s = 0; s < N; s++) {
        X1[s] = X[s];
        X2[s] = X[s];
      }

      X2[j] += std::max(0.01 * X[j], 1.0);

      computeF(F1, X1, C);
      computeF(F2, X2, C);

      for (int i = 0; i < N; i++) {
        J[i + j * N] = (F2[i] - F1[i]) / (X2[j] - X1[j]);
      }
    }

    // Compute F
    computeF(F, X, C);

    // Solve J*dX = F, but note that the true system is J*dX = -F, so we invert the
    // dX vector below.
    dgesv_((int*)&N, (int*)&NRHS, J, (int*)&N, IPIV, F, (int*)&N, &INFO);

    if (INFO != 0) {
#if 1 // Could not solve
      const std::string err = "KMCSolver<R, State, T>::stepImplicitEuler -- could not solve A*x = b";

      pout() << err << endl;
#endif
      converged = false;

      break;
    }
    else {
      // Increment and move on to next iteration if necessary. Note that F = -dX as per the comment above.
      for (int i = 0; i < N; i++) {
        X[i] = X[i] - F[i];
      }
    }

    // Recompute the norm and exit if necessary.
    const Real norm = computeNorm(F, X, C);

    if (norm / initNorm < m_exitTol) {
      break;
    }
  }

  // Turn X into an integer state.
  std::vector<T> outputState(N);

  if (converged) {
    for (int i = 0; i < N; i++) {
      outputState[i] = static_cast<T>(llround(X[i]));
    }
  }
  else {
    // Set X to (what is hopefully) an invalid state and rely on step rejection.
    for (int i = 0; i < N; i++) {
      outputState[i] = static_cast<T>(-1);
    }
  }

  a_state.linearIn(outputState);

  // Cleanup
  delete[] J;
  delete[] X;
  delete[] F;
  delete[] C;

  delete[] X1;
  delete[] X2;
  delete[] F1;
  delete[] F2;

  delete[] IPIV;

  return;
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::advanceTau(State&                   a_state,
                                   const Real&              a_dt,
                                   const KMCLeapPropagator& a_leapPropagator) const noexcept
{
  CH_TIME("KMCSolver::advanceTau(State, ReactionList, Real, KMCLeapPropagator)");

  this->advanceTau(a_state, m_reactions, a_dt, a_leapPropagator);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::advanceTau(State&                   a_state,
                                   const ReactionList&      a_reactions,
                                   const Real&              a_dt,
                                   const KMCLeapPropagator& a_leapPropagator) const noexcept
{
  CH_TIME("KMCSolver::advanceTau(State, ReactionList, Real, KMCLeapPropagator)");

  if (a_reactions.size() > 0) {
    Real curTime = 0.0;
    Real curDt   = a_dt;

    while (curTime < a_dt) {

      // Try big time step first.
      curDt = a_dt - curTime;

      bool valid = false;

      // Substepping so we end up with a valid state.
      while (!valid) {

        State state = a_state;

        // Do a tau-leaping step.
        switch (a_leapPropagator) {
        case KMCLeapPropagator::ExplicitEuler: {
          this->stepExplicitEuler(state, a_reactions, curDt);

          break;
        }
        case KMCLeapPropagator::Midpoint: {
          this->stepMidpoint(state, a_reactions, curDt);

          break;
        }
        case KMCLeapPropagator::PRC: {
          this->stepPRC(state, a_reactions, curDt);

          break;
        }
        case KMCLeapPropagator::ImplicitEuler: {
          this->stepImplicitEuler(state, a_reactions, curDt);

          break;
        }
        default: {
          break;
        }
        }

        // If this was a valid step, accept it. Else reduce dt.
        valid = state.isValidState();

        if (valid) {
          a_state = state;

          curTime += curDt;
        }
        else {
          curDt *= 0.5;
        }
      }
    }
  }
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::advanceHybrid(State&                   a_state,
                                      const Real               a_dt,
                                      const KMCLeapPropagator& a_leapPropagator) const noexcept
{
  CH_TIME("KMCSolver::advanceHybrid(State, Real, KMCLeapPropagator)");

  this->advanceHybrid(a_state, m_reactions, a_dt, a_leapPropagator);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::advanceHybrid(State&                   a_state,
                                      const ReactionList&      a_reactions,
                                      const Real               a_dt,
                                      const KMCLeapPropagator& a_leapPropagator) const noexcept
{
  CH_TIME("KMCSolver::advanceHybrid(State, ReactionList, Real, KMCLeapPropagator)");

  switch (a_leapPropagator) {
  case KMCLeapPropagator::ExplicitEuler: {
    this->advanceHybrid(a_state, a_reactions, a_dt, [this](State& s, const ReactionList& r, const Real dt) {
      this->stepExplicitEuler(s, r, dt);
    });

    break;
  }
  case KMCLeapPropagator::Midpoint: {
    this->advanceHybrid(a_state, a_reactions, a_dt, [this](State& s, const ReactionList& r, const Real dt) {
      this->stepMidpoint(s, r, dt);
    });

    break;
  }
  case KMCLeapPropagator::PRC: {
    this->advanceHybrid(a_state, a_reactions, a_dt, [this](State& s, const ReactionList& r, const Real dt) {
      this->stepPRC(s, r, dt);
    });

    break;
  }
  case KMCLeapPropagator::ImplicitEuler: {
    this->advanceHybrid(a_state, a_reactions, a_dt, [this](State& s, const ReactionList& r, const Real dt) {
      this->stepImplicitEuler(s, r, dt);
    });

    break;
  }
  default: {
    MayDay::Error("KMCSolver::advanceHybrid - unknown leap propagator requested");
  }
  }
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::advanceHybrid(
  State&                                                                               a_state,
  const ReactionList&                                                                  a_reactions,
  const Real                                                                           a_dt,
  const std::function<void(State&, const ReactionList& a_reactions, const Real a_dt)>& a_propagator) const noexcept
{
  CH_TIME("KMCSolver::advanceHybrid(State, ReactionList, Real, std::function)");

  constexpr T one = (T)1;

  // Simulated time within the advancement algorithm.
  Real curTime = 0.0;

  // Outer loop is for reactive substepping over a_dt.
  while (curTime < a_dt) {

    // Partition reactions into critical and non-critical reactions and compute the critical and non-critical time steps.
    const std::pair<ReactionList, ReactionList> partitionedReactions = this->partitionReactions(a_state, a_reactions);

    const ReactionList& criticalReactions    = partitionedReactions.first;
    const ReactionList& nonCriticalReactions = partitionedReactions.second;

    const std::vector<Real> propensitiesCrit    = this->propensities(a_state, criticalReactions);
    const std::vector<Real> propensitiesNonCrit = this->propensities(a_state, nonCriticalReactions);

    Real dtCrit    = this->getCriticalTimeStep(propensitiesCrit);
    Real dtNonCrit = this->getNonCriticalTimeStep(a_state, nonCriticalReactions, propensitiesNonCrit);

    // Try the various advancement algorithms; the loop is for step rejection in case we end up with an invalid state, e.g.
    // a state with a negative number of particles.
    bool validStep = false;

    while (!validStep) {

      // Do a backup of the advancement state to operate on. This is necessary because the tau-leaping
      // algorithms advance the state but we may need to reject those steps if we end up with a thermodynamically
      // invalid state.
      State state = a_state;

      // Compute the time step to be used.
      const Real curDt = std::min(a_dt - curTime, std::min(dtCrit, dtNonCrit));

      // Are we only doing non-critical reactions?
      const bool nonCriticalOnly = (dtNonCrit < dtCrit) || (criticalReactions.size() == 0) ||
                                   (dtCrit > (a_dt - curTime));

      if (nonCriticalOnly) {

        // Compute the total propensity because it might actually be better to run SSA steps.
        Real A = this->totalPropensity(state, a_reactions);

        if (A * curDt < m_SSAlim) {
          // TLDR: Tau-leaping is inefficient here so we switch to an SSA-based algorithm for the WHOLE reaction set.

          // Number of SSA steps taken and simulated time within the SSA. We will compute until either dtSSA < curDt or
          // we've exceeded the maximum number of SSA steps that the user has permitted (m_numSSA).
          Real dtSSA  = 0.0;
          T    numSSA = 0;

          while (dtSSA < curDt && numSSA < m_numSSA) {

            // Recompute propensities for the full reaction set and advance everything using the SSA.
            std::vector<Real> propensities = this->propensities(a_state, a_reactions);

            A = 0.0;
            for (size_t i = 0; i < propensities.size(); i++) {
              A += propensities[i];
            }

            // Compute the time to the next reaction.
            const Real dtReact = this->getCriticalTimeStep(A);

            if (dtSSA + dtReact < curDt) {
              this->stepSSA(a_state, a_reactions, propensities);

              dtSSA += dtReact;
              numSSA += one;
            }
            else {

              // Next reaction occured outside the substep -- break out of the loop.
              dtSSA = curDt;
            }
          }

          validStep = true;
          curTime += dtSSA;
        }
        else {
          // Do the tau-leaping.
          a_propagator(state, nonCriticalReactions, curDt);

          // Check if we need to reject the state and rather try again with a smaller non-critical time step.
          validStep = state.isValidState();

          if (validStep) {
            a_state = state;

            curTime += curDt;
          }
          else {
            dtNonCrit *= 0.5;
          }
        }
      }
      else {
        // TLDR: Here, one critical reaction fired. We advance the critical reaction using the SSA and the rest of the reactions
        //       using tau-leaping.
        this->stepSSA(state, criticalReactions, propensitiesCrit);

        // Do the tau-leaping over the non-critical reactions.
        a_propagator(state, nonCriticalReactions, curDt);

        // Check if we need to reject the state and rather try again with a smaller non-critical time step.
        validStep = state.isValidState();

        if (validStep) {
          a_state = state;

          curTime += curDt;
        }
        else {
          dtNonCrit *= 0.5;
        }
      }
    }
  }
}

#include <CD_NamespaceFooter.H>

#endif
