/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_KMCSolverImplem.H
 * @brief  Implementation of CD_KMCSolver.H
 * @author Robert Marskar
 */

#ifndef CD_KMCSOLVERIMPLEM_H
#define CD_KMCSOLVERIMPLEM_H

// Std includes
#include <limits>
#include <unordered_set>

// Chombo includes
#include <CH_assert.H>
#include <MayDay.H>
#include <parstream.H>

// Our includes
#include <CD_Random.H>
#include <CD_KMCSolver.H>
#include <CD_LaPackUtils.H>
#include <CD_NamespaceHeader.H>

template <typename R, typename State, typename T>
inline KMCSolver<R, State, T>::KMCSolver() noexcept
{
  this->setSolverParameters(0, 0, 100, std::numeric_limits<Real>::max(), 0.0, 1.E-6);
}

template <typename R, typename State, typename T>
inline KMCSolver<R, State, T>::KMCSolver(const ReactionList& a_reactions) noexcept
{
  this->define(a_reactions);
}

template <typename R, typename State, typename T>
inline KMCSolver<R, State, T>::~KMCSolver() noexcept
{}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::define(const ReactionList& a_reactions) noexcept
{
  m_reactions = a_reactions;

  // Default settings. These are equivalent to ALWAYS using tau-leaping.
  this->setSolverParameters(0, 0, 100, std::numeric_limits<Real>::max(), 0.0, 1.E-6);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::setSolverParameters(const T    a_numCrit,
                                            const T    a_numSSA,
                                            const T    a_maxIter,
                                            const Real a_eps,
                                            const Real a_SSAlim,
                                            const Real a_exitTol) noexcept
{
  m_Ncrit   = a_numCrit;
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
  return this->propensities(a_state, m_reactions);
}

template <typename R, typename State, typename T>
inline std::vector<Real>
KMCSolver<R, State, T>::propensities(const State& a_state, const ReactionList& a_reactions) const noexcept
{
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
  return this->totalPropensity(a_state, m_reactions);
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::totalPropensity(const State& a_state, const ReactionList& a_reactions) const noexcept
{
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
  return this->partitionReactions(a_state, m_reactions);
}

template <typename R, typename State, typename T>
inline std::pair<typename KMCSolver<R, State, T>::ReactionList, typename KMCSolver<R, State, T>::ReactionList>
KMCSolver<R, State, T>::partitionReactions(const State& a_state, const ReactionList& a_reactions) const noexcept
{
  ReactionList criticalReactions;
  ReactionList nonCriticalReactions;

  const size_t numReactions = a_reactions.size();

  // Each reaction lands in exactly one of the two lists, so reserving the full count up front avoids reallocations
  // during the emplace_back calls below.
  criticalReactions.reserve(numReactions);
  nonCriticalReactions.reserve(numReactions);

  for (size_t i = 0; i < numReactions; i++) {
    const T Lj = a_reactions[i]->computeCriticalNumberOfReactions(a_state);

    if (Lj < m_Ncrit) {
      criticalReactions.emplace_back(a_reactions[i]);
    }
    else {
      nonCriticalReactions.emplace_back(a_reactions[i]);
    }
  }

  // Move rather than copy: this runs once per substep per grid cell, and copying the two lists would
  // duplicate both allocations and bump the refcount of every reaction twice more.
  return std::make_pair(std::move(criticalReactions), std::move(nonCriticalReactions));
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
  const Real u = std::numeric_limits<Real>::min() + Random::getUniformReal01();

  return log(1.0 / u) / a_totalPropensity;
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::getNonCriticalTimeStep(const State& a_state) const noexcept
{
  const auto& partitionedReactions = this->partitionReactions(a_state, m_reactions);

  const std::vector<Real> propensities = this->propensities(a_state, partitionedReactions.second);

  return this->getNonCriticalTimeStep(a_state, m_reactions, partitionedReactions.second, propensities);
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::getNonCriticalTimeStep(const State& a_state, const ReactionList& a_reactions) const noexcept
{
  const std::vector<Real> propensities = this->propensities(a_state, a_reactions);

  return this->getNonCriticalTimeStep(a_state, a_reactions, a_reactions, propensities);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::gatherDistinctReactants(const ReactionList& a_reactions) const noexcept
{
  // Size the scratch from the largest reactant index the reactions refer to. Species above that are
  // not reactants of any of these reactions and can never be reduced over below.
  size_t numSpecies = 0;

  for (const auto& reaction : a_reactions) {
    for (const size_t reactant : reaction->getReactants()) {
      numSpecies = std::max(numSpecies, reactant + 1);
    }
  }

  if (m_seenScratch.size() < numSpecies) {
    m_seenScratch.resize(numSpecies, 0);
    m_muScratch.resize(numSpecies, 0.0);
    m_sigmaScratch.resize(numSpecies, 0.0);
    m_lossScratch.resize(numSpecies, 0.0);
  }

  // Collect the distinct reactants. The stamps make this a single pass with no sort: a reactant is
  // appended the first time it is seen and skipped afterwards.
  m_reactantScratch.clear();

  for (const auto& reaction : a_reactions) {
    for (const size_t reactant : reaction->getReactants()) {
      if (!m_seenScratch[reactant]) {
        m_seenScratch[reactant] = 1;

        m_reactantScratch.push_back(reactant);
      }
    }
  }

  // Leave the stamps clear for the next call.
  for (const size_t reactant : m_reactantScratch) {
    m_seenScratch[reactant] = 0;
  }
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::getNonCriticalTimeStep(const State&             a_state,
                                               const ReactionList&      a_reactions,
                                               const ReactionList&      a_nonCriticalReactions,
                                               const std::vector<Real>& a_nonCriticalPropensities) const noexcept
{
  CH_assert(a_nonCriticalReactions.size() == a_nonCriticalPropensities.size());
  CH_assert(a_nonCriticalReactions.size() <= a_reactions.size());

  constexpr Real one = 1.0;

  Real dt = std::numeric_limits<Real>::max();

  const size_t numReactions = a_nonCriticalReactions.size();

  if (numReactions > 0) {

    // 1. Gather the distinct reactants of ALL reactions. A reactant of a critical reaction that the non-critical
    // reactions produce must be bounded too, since the critical propensities are frozen during the leap and must
    // not drift more than the leap condition allows.
    this->gatherDistinctReactants(a_reactions);

    for (const size_t reactant : m_reactantScratch) {
      m_muScratch[reactant]    = 0.0;
      m_sigmaScratch[reactant] = 0.0;
      m_lossScratch[reactant]  = 0.0;
    }

    // 2. Accumulate the expected net change, its variance, and the expected consumption per reactant.
    // The consumption only collects the reactions that remove the reactant, so it is the loss side
    // of the net change without the production that may cancel it. The loops are ordered reaction-
    // major so that each reaction's state change is fetched once per reactant from its dense table
    // rather than through a map lookup; the sums are the same as reactant-major.
    for (size_t i = 0; i < numReactions; i++) {
      const Real& p = a_nonCriticalPropensities[i];

      for (const size_t reactant : m_reactantScratch) {
        const auto muIJ = a_nonCriticalReactions[i]->getStateChange(reactant);

        m_muScratch[reactant] += muIJ * p;
        m_sigmaScratch[reactant] += muIJ * muIJ * p;

        if (muIJ < 0) {
          m_lossScratch[reactant] -= muIJ * p;
        }
      }
    }

    // 3. Reduce over the reactants.
    for (const size_t reactant : m_reactantScratch) {

      // Xi is the population of the current reactant. It might seem weird that we are indexing this
      // through the reactions rather then through the state. Which it is, but the reason for that is
      // that it is the REACTION that determines how we index the state population. This is simply a
      // design choice that permits the user to apply different type of reactions without changing
      // the underlying state.
      const T Xi = a_nonCriticalReactions[0]->population(reactant, a_state);

      // Set gi to 1 for now. A more complex version would parse this through an input parameter
      // where the user has inspected the highest-order-reaction.
      constexpr Real gi = 1.0;

      Real dt1 = std::numeric_limits<Real>::max();
      Real dt2 = std::numeric_limits<Real>::max();
      Real dt3 = std::numeric_limits<Real>::max();

      // The floor of one applies at every population, so a reactant at zero is bounded to one
      // expected net change rather than left unbounded. This is what lets the critical reactions
      // respond to a species that the leap produces from nothing.
      const Real f = std::max(m_eps * Xi / gi, one);

      const Real mu     = std::abs(m_muScratch[reactant]);
      const Real sigma2 = std::abs(m_sigmaScratch[reactant]);
      const Real loss   = m_lossScratch[reactant];

      // Leap condition on the mean and variance of the net change.
      if (mu > std::numeric_limits<Real>::min()) {
        dt1 = f / mu;
      }
      if (sigma2 > std::numeric_limits<Real>::min()) {
        dt2 = (f * f) / sigma2;
      }

      // Bound on the gross consumption: the reactions that remove this reactant may not, between
      // them, be expected to fire more than f times. This is what keeps the Poisson draws of a fast
      // reaction below the population it draws on when a producing reaction hides its loss from mu.
      if (loss > std::numeric_limits<Real>::min()) {
        dt3 = f / loss;
      }

      dt = std::min(dt, std::min(dt1, std::min(dt2, dt3)));
    }
  }

  return dt;
}

template <typename R, typename State, typename T>
inline Real
KMCSolver<R, State, T>::computeDt(const State&             a_state,
                                  const ReactionList&      a_reactions,
                                  const std::vector<Real>& a_propensities,
                                  const Real               a_epsilon) const noexcept
{
  CH_assert(a_reactions.size() == a_propensities.size());

  constexpr Real one = 1.0;

  Real dt = std::numeric_limits<Real>::max();

  const size_t numReactions = a_reactions.size();

  if (numReactions > 0) {

    // 1. Gather the distinct reactants involved in the input reactions.
    this->gatherDistinctReactants(a_reactions);

    for (const size_t reactant : m_reactantScratch) {
      m_muScratch[reactant] = 0.0;
    }

    // 2. Accumulate the expected change per reactant. The loops are ordered reaction-major so that
    // each reaction's state change is fetched once per reactant from its dense table rather than
    // through a map lookup; the sums are the same as reactant-major.
    for (size_t i = 0; i < numReactions; i++) {
      const Real& p = a_propensities[i];

      for (const size_t reactant : m_reactantScratch) {
        const auto muIJ = a_reactions[i]->getStateChange(reactant);

        m_muScratch[reactant] += muIJ * p;
      }
    }

    // 3. Reduce over the reactants.
    for (const size_t reactant : m_reactantScratch) {

      // Xi is the population of the current reactant. It might seem weird that we are indexing this
      // through the reactions rather then through the state. Which it is, but the reason for that is
      // that it is the REACTION that determines how we index the state population. This is simply a
      // design choice that permits the user to apply different type of reactions without changing
      // the underlying state.
      const T Xi = R::population(reactant, a_state);

      if (Xi > (T)0) {

        // Set gi to 1 for now. A more complex version would parse this through an input parameter
        // where the user has inspected the highest-order-reaction.
        constexpr Real gi = 1.0;

        const Real f  = std::max(a_epsilon * Xi / gi, one);
        const Real mu = std::abs(m_muScratch[reactant]);

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
  this->stepSSA(a_state, m_reactions);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepSSA(State& a_state, const ReactionList& a_reactions) const noexcept
{
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

    size_t r = numReactions - 1;

    Real sumProp = 0.0;
    for (size_t i = 0; i + 1 < numReactions; i++) {
      sumProp += a_propensities[i];

      if (sumProp >= u * A) {
        r = i;

        break;
      }
    }

    CH_assert(r < a_reactions.size());

    // Advance by one reaction.
    a_reactions[r]->advanceState(a_state, one);
  }
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::advanceSSA(State& a_state, const Real a_dt) const noexcept
{
  this->advanceSSA(a_state, m_reactions, a_dt);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::advanceSSA(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept
{
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
  this->stepExplicitEuler(a_state, m_reactions, a_dt);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepExplicitEuler(State&              a_state,
                                          const ReactionList& a_reactions,
                                          const Real          a_dt) const noexcept
{
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
  this->stepMidpoint(a_state, m_reactions, a_dt);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepMidpoint(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept
{
  const int numReactions = a_reactions.size();

  if (numReactions > 0) {

    std::vector<Real> propensities = this->propensities(a_state, a_reactions);

    State Xdagger = a_state;

    for (size_t i = 0; i < numReactions; i++) {
      // TLDR: Predict a midpoint state -- unfortunately this means that as a_dt->0 we end up with plain
      //       tau-leaping. I don't know of a way to fix this without introducing double fluctuations (yet).
      a_reactions[i]->advanceState(Xdagger, (T)std::round(0.5 * propensities[i] * a_dt));
    }

    propensities = this->propensities(Xdagger, a_reactions);

    for (size_t i = 0; i < numReactions; i++) {
      const T curReactions = (T)Random::getPoisson<long long>(propensities[i] * a_dt);

      a_reactions[i]->advanceState(a_state, curReactions);
    }
  }
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepPRC(State& a_state, const Real a_dt) const noexcept
{
  this->stepPRC(a_state, m_reactions, a_dt);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepPRC(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept
{
  const int numReactions = a_reactions.size();

  if (numReactions > 0) {

    std::vector<Real> aj = this->propensities(a_state, a_reactions);

    const std::vector<Real> ak = aj;

    for (int j = 0; j < numReactions; j++) {
      for (int k = 0; k < numReactions; k++) {
        State x = a_state;

        a_reactions[k]->advanceState(x, (T)1);

        const Real etajk = a_reactions[j]->propensity(x) - ak[j];

        aj[j] += 0.5 * a_dt * ak[k] * etajk;
      }
    }

    for (size_t i = 0; i < numReactions; i++) {
      const T nr = (T)Random::getPoisson<long long>(aj[i] * a_dt);

      a_reactions[i]->advanceState(a_state, nr);
    }
  }
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepImplicitEuler(State& a_state, const Real a_dt) const noexcept
{
  this->stepImplicitEuler(a_state, m_reactions, a_dt);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::stepImplicitEuler(State&              a_state,
                                          const ReactionList& a_reactions,
                                          const Real          a_dt) const noexcept
{
  // TLDR: The implicit Euler tau-leaping scheme is equivalent to the solution of
  //
  // F(X) = X - (x + sum_j (nu_j * (P(a(x)*dt) - a(x)*dt)]) - dt*sum_j nu_j a_j(X)
  //      = X - c - dt*sum_j nu_j * a_j(X),
  //      = 0
  //
  // where c = (x + sum_j (nu_j * (P(a(x)*dt) - a(x)*dt)]) is a constant term throughout the Newton iterations. This
  // term is presampled before the Newton iterations begin.

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

  // Allocate memory for the Jacobian (J), the Newton increment (X), and F(X) (F). Also include the constant term from
  // the Poisson sampling (c)
  std::vector<double> J(static_cast<size_t>(N * N));
  std::vector<double> X(static_cast<size_t>(N));
  std::vector<double> F(static_cast<size_t>(N));
  std::vector<double> C(static_cast<size_t>(N));

  // Temporary storage used for computing the Jacobian matrix.
  std::vector<double> X1(static_cast<size_t>(N));
  std::vector<double> X2(static_cast<size_t>(N));
  std::vector<double> F2(static_cast<size_t>(N));

  // Things that are required by LAPACK
  int              INFO = 0;
  int              NRHS = 1;
  std::vector<int> IPIV(static_cast<size_t>(N));

  // Compute the constant term.
  compConstantTerm(C.data(), X.data(), a_state, a_dt);

  // Compute the max-norm of F(0) to use as an exit criterion.
  for (int i = 0; i < N; i++) {
    X1[i] = 0.0;
  }

  const Real initNorm = computeNorm(F.data(), X1.data(), C.data());

  bool converged = true;

  for (int k = 0; k < m_maxIter; k++) {

    // Compute the residual F(X) at the current iterate. It is the unperturbed column used in every finite-difference
    // Jacobian column (it does not depend on the perturbation direction j) and is also the right-hand side for the
    // linear solve below, so it is computed once per Newton iteration rather than N+1 times.
    computeF(F.data(), X.data(), C.data());

    // Numerically computed the Jacobian using finite differences. The Jacobian is given by
    // Jij = dF_i/dx_j
    for (int j = 0; j < N; j++) {

      for (int s = 0; s < N; s++) {
        X2[s] = X[s];
      }

      X2[j] += std::max(0.01 * X[j], 1.0);

      computeF(F2.data(), X2.data(), C.data());

      for (int i = 0; i < N; i++) {
        J[i + j * N] = (F2[i] - F[i]) / (X2[j] - X[j]);
      }
    }

    // Solve J*dX = F, but note that the true system is J*dX = -F, so we invert the
    // dX vector below.
    dgesv_((int*)&N, &NRHS, J.data(), (int*)&N, IPIV.data(), F.data(), (int*)&N, &INFO);

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
    const Real norm = computeNorm(F.data(), X.data(), C.data());

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
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::advanceTau(State&                   a_state,
                                   const Real&              a_dt,
                                   const KMCLeapPropagator& a_leapPropagator) const noexcept
{
  this->advanceTau(a_state, m_reactions, a_dt, a_leapPropagator);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::advanceTau(State&                   a_state,
                                   const ReactionList&      a_reactions,
                                   const Real&              a_dt,
                                   const KMCLeapPropagator& a_leapPropagator) const noexcept
{
  if (a_reactions.size() > 0) {
    Real curTime = 0.0;

    while (curTime < a_dt) {

      // Bound the substep by the tau-leaping condition tau = eps * X_i / |sum_j nu_ij * a_j|, as advanceHybrid
      // does for its non-critical reactions. The propensities are recomputed each substep because the state
      // changes.
      const std::vector<Real> propensities = this->propensities(a_state, a_reactions);

      const Real dtLeap = this->getNonCriticalTimeStep(a_state, a_reactions, a_reactions, propensities);

      Real curDt = std::min(a_dt - curTime, dtLeap);

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
  this->advanceHybrid(a_state, m_reactions, a_dt, a_leapPropagator);
}

template <typename R, typename State, typename T>
inline void
KMCSolver<R, State, T>::advanceHybrid(State&                   a_state,
                                      const ReactionList&      a_reactions,
                                      const Real               a_dt,
                                      const KMCLeapPropagator& a_leapPropagator) const noexcept
{
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
  constexpr T one = (T)1;

  // Simulated time within the advancement algorithm.
  Real curTime = 0.0;

  // Outer loop is for reactive substepping over a_dt.
  while (curTime < a_dt) {

    // Partition reactions into critical and non-critical reactions. The propensities are those of the state at the
    // start of the substep: the critical firing (if any) is selected from them, and the non-critical leap starts from
    // the same state.
    const std::pair<ReactionList, ReactionList> partitionedReactions = this->partitionReactions(a_state, a_reactions);

    const ReactionList& criticalReactions    = partitionedReactions.first;
    const ReactionList& nonCriticalReactions = partitionedReactions.second;

    const std::vector<Real> propensitiesCrit    = this->propensities(a_state, criticalReactions);
    const std::vector<Real> propensitiesNonCrit = this->propensities(a_state, nonCriticalReactions);

    // Leap candidate from the tau-leaping condition on the non-critical reactions. Halved on every rejected step.
    Real dtNonCrit = this->getNonCriticalTimeStep(a_state, a_reactions, nonCriticalReactions, propensitiesNonCrit);

    // Total propensity of the whole reaction set, which decides between tau-leaping and the SSA below.
    const Real A = this->totalPropensity(a_state, a_reactions);

    // The loop is for step rejection in case we end up with an invalid state, e.g. a state with a negative number of
    // particles.
    bool validStep = false;

    while (!validStep) {

      // The leap candidate, bounded by the remaining time. It is a function of the state only, so the choice between
      // the SSA and tau-leaping below is made without reference to any random draw.
      const Real dtLeap = std::min(a_dt - curTime, dtNonCrit);

      // Fall back to an exact SSA advancement of the WHOLE reaction set when the leap is expected to carry fewer than
      // m_SSAlim firings; in that limit the SSA is both cheaper and more accurate than tau-leaping. This is only done
      // when the user has permitted at least one SSA step (m_numSSA > 0); otherwise we always tau-leap, which also
      // avoids a no-progress infinite loop when m_numSSA == 0.
      const bool useSSA = (m_numSSA >= one) && (A * dtLeap < m_SSAlim);

      if (useSSA) {
        // Advance with the SSA until either dtLeap has elapsed or m_numSSA reactions have fired. A firing whose
        // waiting time reaches beyond dtLeap is not applied; the state is advanced to dtLeap and the waiting time is
        // drawn again in the next substep, which is exact because the waiting time is memoryless.
        Real dtSSA  = 0.0;
        T    numSSA = 0;

        while (dtSSA < dtLeap && numSSA < m_numSSA) {

          // Recompute propensities for the full reaction set and advance everything using the SSA.
          const std::vector<Real> propensities = this->propensities(a_state, a_reactions);

          const Real dtReact = this->getCriticalTimeStep(propensities);

          if (dtSSA + dtReact < dtLeap) {
            this->stepSSA(a_state, a_reactions, propensities);

            dtSSA += dtReact;
            numSSA += one;
          }
          else {
            dtSSA = dtLeap;
          }
        }

        validStep = true;
        curTime += dtSSA;
      }
      else {
        // Waiting time to the next critical firing. Drawn only now, after the SSA decision, and drawn again on every
        // retry, so that every draw is either used or discarded together with the step it belonged to.
        const Real dtCrit = (criticalReactions.size() > 0) ? this->getCriticalTimeStep(propensitiesCrit)
                                                           : std::numeric_limits<Real>::max();

        // One critical reaction fires if it is due before the leap candidate runs out; the substep then ends at the
        // firing.
        const bool fireCritical = (dtCrit < dtLeap);

        const Real curDt = fireCritical ? dtCrit : dtLeap;

        // Operate on a copy so that the step can be rejected.
        State state = a_state;

        // Tau-leap the non-critical reactions over curDt. This is done FIRST, while the state is still at the start of
        // the substep, so that the leap propensities are evaluated at the start-of-substep state. The single critical
        // firing is then applied as a fixed stoichiometric increment; its reaction is selected from the
        // start-of-substep propensities (propensitiesCrit), so applying it last is exact. Both contributions are thus
        // consistent with the start-of-substep state, as required by the Cao-Gillespie-Petzold scheme.
        a_propagator(state, nonCriticalReactions, curDt);

        if (fireCritical) {
          this->stepSSA(state, criticalReactions, propensitiesCrit);
        }

        // Check if we need to reject the state and rather try again with a smaller leap candidate.
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
