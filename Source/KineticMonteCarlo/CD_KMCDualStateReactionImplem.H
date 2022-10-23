/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCDualStateReactionImplem.H
  @brief  Implementation of CD_KMCDualStateReaction.H
  @author Robert Marskar
*/

#ifndef CD_KMCDualStateReactionImplem_H
#define CD_KMCDualStateReactionImplem_H

// Our includes
#include <CD_KMCDualStateReaction.H>
#include <CD_NamespaceHeader.H>

template <typename T, typename State>
inline KMCDualStateReaction<T, State>::KMCDualStateReaction(const std::list<size_t>& a_lhsReactives,
                                                            const std::list<size_t>& a_rhsReactives,
                                                            const std::list<size_t>& a_rhsNonReactives) noexcept
{
  m_lhsReactives    = a_lhsReactives;
  m_rhsReactives    = a_rhsReactives;
  m_rhsNonReactives = a_rhsNonReactives;

  CH_assert(a_lhsReactives > 0);

  this->computeStateChanges();
}

template <typename T, typename State>
inline KMCDualStateReaction<T, State>::~KMCDualStateReaction()
{}

template <typename T, typename State>
inline void
KMCDualStateReaction<T, State>::computeStateChanges() noexcept
{
  // Consumed species.
  for (const auto& r : m_lhsReactives) {
    if (m_reactiveStateChange.find(r) == m_reactiveStateChange.end()) {
      m_reactiveStateChange.emplace(r, -1);
    }
    else {
      m_reactiveStateChange[r]--;
    }
  }

  // Produced species.
  for (const auto& p : m_rhsReactives) {
    if (m_reactiveStateChange.find(p) == m_reactiveStateChange.end()) {
      m_reactiveStateChange.emplace(p, +1);
    }
    else {
      m_reactiveStateChange[p]++;
    }
  }

  // Produced photons.
  for (const auto& p : m_rhsNonReactives) {
    if (m_nonReactiveStateChange.find(p) == m_nonReactiveStateChange.end()) {
      m_nonReactiveStateChange.emplace(p, +1);
    }
    else {
      m_nonReactiveStateChange[p]++;
    }
  }

  // Compute the propensity factor that we need when there are bi- or tri-particle reactions involving the same species. We need to do this
  // because for biparticle reactions of N particles of the same type, there are 0.5 * N * (N-1) unique pairs of particles. Or in general when
  // we have a k-order reaction we have 'N choose k' = N!/(k! * (N-k)!) combinations. In the propensity function we compute the propensity as
  //
  //    a = rate * factor * N * (N-1) * (N-2) ...
  //
  // For this to make sense we have
  //
  //    a = rate * N * (N-1) * (N-2) ... (N-k+1) * (N-k)!/(k! * (N-k)!)
  //      = rate * N * (N-1) * (N-2) ... * 1/k!
  //
  // so the factor is just 1/k! (for each species).

  m_propensityFactor = 1.0;

  std::map<size_t, size_t> reactantNumbers;
  for (const auto& r : m_lhsReactives) {
    if (reactantNumbers.find(r) == reactantNumbers.end()) {
      reactantNumbers.emplace(r, 1);
    }
    else {
      reactantNumbers[r]++;
    }
  }

  // Factorial function.
  T factorial = [](const T& N) -> T {
    T fac = 1;
    for (T i = 2; i <= N; i++) {
      fac *= i;
    }

    return fac;
  };

  for (const auto& rn : reactantNumbers) {
    m_propensityFactor *= 1.0 / factorial(rn.second);
  }
}

template <typename T, typename State>
inline Real&
KMCDualStateReaction<T, State>::rate() const noexcept
{
  return m_rate;
}

template <typename T, typename State>
inline T
KMCDualStateReaction<T, State>::population(const size_t& a_reactant, const State& a_state) const noexcept
{
  const auto& reactiveState = a_state.getReactiveState();

  CH_assert(reactiveState.size() > a_reactant);

  return reactiveState[a_reactant];
}

template <typename T, typename State>
inline Real
KMCDualStateReaction<T, State>::propensity(const State& a_state) const noexcept
{
#ifndef NDEBUG
  this->sanityCheck(a_state);
#endif

  Real A = m_rate * m_propensityFactor;

  auto reactiveState = a_state.getReactiveState();

  for (const auto& r : m_lhsReactives) {
    A *= reactiveState[r];

    reactiveState[r]--;
  }

  return A;
}

template <typename T, typename State>
inline T
KMCDualStateReaction<T, State>::computeCriticalNumberOfReactions(const State& a_state) const noexcept
{
#ifndef NDEBUG
  this->sanityCheck(a_state);
#endif

  T Lj = std::numeric_limits<long long>::max();

  const auto& reactiveState = a_state.getReactiveState();

  for (const auto& s : m_reactiveStateChange) {

    const size_t rI   = s.first;
    const T      nuIJ = s.second;

    if (nuIJ < 0) {
      Lj = std::min(Lj, reactiveState[rI] / std::abs(nuIJ));
    }
  }

  return Lj;
}

template <typename T, typename State>
inline std::list<size_t>
KMCDualStateReaction<T, State>::getReactants() const noexcept
{
  return m_lhsReactives;
}

template <typename T, typename State>
inline T
KMCDualStateReaction<T, State>::getStateChange(const size_t a_particleReactant) const noexcept
{
  T nuIJ = 0;

  if (m_reactiveStateChange.find(a_particleReactant) != m_reactiveStateChange.end()) {
    nuIJ = m_reactiveStateChange.at(a_particleReactant);
  }

  return nuIJ;
}

template <typename T, typename State>
inline void
KMCDualStateReaction<T, State>::advanceState(State& a_state, const T& a_numReactions) const noexcept
{
#ifndef NDEBUG
  this->sanityCheck(a_state);
#endif

  auto& reactiveState = a_state.getReactiveState();
  auto& photonState   = a_state.getNonReactiveState();

  for (const auto& s : m_reactiveStateChange) {
    reactiveState[s.first] += a_numReactions * s.second;
  }

  for (const auto& s : m_nonReactiveStateChange) {
    photonState[s.first] += a_numReactions * s.second;
  }
}

template <typename T, typename State>
inline void
KMCDualStateReaction<T, State>::sanityCheck(const State& a_state) const noexcept
{
  const auto& reactiveState = a_state.getReactiveState();
  const auto& photonState   = a_state.getNonReactiveState();

  for (const auto& idx : m_lhsReactives) {
    CH_assert(reactiveState.size() > idx);
  }

  for (const auto& idx : m_rhsReactives) {
    CH_assert(reactiveState.size() > idx);
  }

  for (const auto& idx : m_rhsNonReactives) {
    CH_assert(photonState.size() > idx);
  }
}

#include <CD_NamespaceFooter.H>

#endif
