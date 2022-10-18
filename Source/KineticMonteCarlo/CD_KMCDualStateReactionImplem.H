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
#include <CD_KMCDualStateReactionState.H>
#include <CD_NamespaceHeader.H>

template <typename T, typename State>
inline KMCDualStateReaction<T, State>::KMCDualStateReaction(const std::list<size_t>& a_particleReactants,
                                                            const std::list<size_t>& a_particleProducts,
                                                            const std::list<size_t>& a_photonProducts) noexcept
{
  m_particleReactants = a_particleReactants;
  m_particleProducts  = a_particleProducts;
  m_photonProducts    = a_photonProducts;

  CH_assert(a_particleReactants > 0);

  this->computeStateChanges();
}

template <typename T, typename State>
inline KMCDualStateReaction<T, State>::~KMCDualStateReaction()
{}

template <typename T, typename State>
inline void
KMCDualStateReaction<T, State>::computeStateChanges() noexcept
{
  // Consumed particles.
  for (const auto& r : m_particleReactants) {
    if (m_particleStateChange.find(r) == m_particleStateChange.end()) {
      m_particleStateChange.emplace(r, -1);
    }
    else {
      m_particleStateChange[r]--;
    }
  }

  // Produced particles.
  for (const auto& p : m_particleProducts) {
    if (m_particleStateChange.find(p) == m_particleStateChange.end()) {
      m_particleStateChange.emplace(p, +1);
    }
    else {
      m_particleStateChange[p]++;
    }
  }

  // Produced photons.
  for (const auto& p : m_photonProducts) {
    if (m_photonStateChange.find(p) == m_photonStateChange.end()) {
      m_photonStateChange.emplace(p, +1);
    }
    else {
      m_photonStateChange[p]++;
    }
  }

  // Compute the reaction order.
  m_reactionOrder = (T)a_particleReactants.size();
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
  const auto& particleState = a_state.getReactiveState();

  CH_assert(particleState.size() > a_reactant);

  return particleState[a_reactant];
}

template <typename T, typename State>
inline Real
KMCDualStateReaction<T, State>::propensity(const State& a_state) const noexcept
{
#ifndef NDEBUG
  this->sanityCheck(a_state);
#endif

  Real A = m_rate;

  const auto& particleState = a_state.getReactiveState();

  for (const auto& r : m_particleReactants) {
    A *= particleState[r];
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

  const auto& particleState = a_state.getReactiveState();

  for (const auto& s : m_particleStateChange) {

    const size_t rI   = s.first;
    const T      nuIJ = s.second;

    if (nuIJ < 0) {
      Lj = std::min(Lj, particleState[rI] / std::abs(nuIJ));
    }
  }

  return Lj;
}

template <typename T, typename State>
inline std::list<size_t>
KMCDualStateReaction<T, State>::getReactants() const noexcept
{
  return m_particleReactants;
}

template <typename T, typename State>
inline T
KMCDualStateReaction<T>::getStateChange(const size_t a_particleReactant) const noexcept
{
  T nuIJ = 0;

  if (m_particleStateChange.find(a_particleReactant) != m_particleStateChange.end()) {
    nuIJ = m_particleStateChange.at(a_particleReactant);
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

  auto& particleState = a_state.getReactiveState();
  auto& photonState   = a_state.getNonReactiveState();

  for (const auto& s : m_particleStateChange) {
    particleState[s.first] += a_numReactions * s.second;
  }

  for (const auto& s : m_photonStateChange) {
    photonState[s.first] += a_numReactions * s.second;
  }
}

template <typename T, typename State>
inline void
KMCDualStateReaction<T, State>::sanityCheck(const State& a_state) const noexcept
{
  const auto& particleState = a_state.getReactiveState();
  const auto& photonState   = a_state.getNonReactiveState();

  for (const auto& idx : m_particleReactants) {
    CH_assert(particleState.size() > idx);
  }

  for (const auto& idx : m_particleProducts) {
    CH_assert(particleState.size() > idx);
  }

  for (const auto& idx : m_photonProducts) {
    CH_assert(photonState.size() > idx);
  }
}

#include <CD_NamespaceFooter.H>

#endif
