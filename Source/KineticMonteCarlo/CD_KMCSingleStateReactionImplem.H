/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCSingleStateReactionImplem.H
  @brief  Implementation of CD_KMCSingleStateReaction.H
  @author Robert Marskar
*/

#ifndef CD_KMCSingleStateReactionImplem_H
#define CD_KMCSingleStateReactionImplem_H

// Our includes
#include <CD_KMCSingleStateReaction.H>
#include <CD_NamespaceHeader.H>

template <typename T, typename State>
inline KMCSingleStateReaction<T, State>::KMCSingleStateReaction(const std::list<size_t>& a_reactants,
                                                                const std::list<size_t>& a_products) noexcept
{
  m_reactants = a_reactants;
  m_products  = a_products;

  CH_assert(m_reactants.size() > 0);

  this->computeStateChanges();
}

template <typename T, typename State>
inline KMCSingleStateReaction<T, State>::~KMCSingleStateReaction()
{}

template <typename T, typename State>
inline void
KMCSingleStateReaction<T, State>::computeStateChanges() noexcept
{
  // Consumed particles.
  for (const auto& r : m_reactants) {
    if (m_stateChange.find(r) == m_stateChange.end()) {
      m_stateChange.emplace(r, -1);
    }
    else {
      m_stateChange[r]--;
    }
  }

  // Produced particles.
  for (const auto& p : m_products) {
    if (m_stateChange.find(p) == m_stateChange.end()) {
      m_stateChange.emplace(p, +1);
    }
    else {
      m_stateChange[p]++;
    }
  }
}

template <typename T, typename State>
inline Real&
KMCSingleStateReaction<T, State>::rate() const noexcept
{
  return m_rate;
}

template <typename T, typename State>
inline T
KMCSingleStateReaction<T, State>::population(const size_t& a_reactant, const State& a_state) const noexcept
{
  return a_state[a_reactant];
}

template <typename T, typename State>
inline Real
KMCSingleStateReaction<T, State>::propensity(const State& a_state) const noexcept
{
  Real A = m_rate;

  for (const auto& r : m_reactants) {
    A *= a_state[r];
  }

  return A;
}

template <typename T, typename State>
inline T
KMCSingleStateReaction<T, State>::computeCriticalNumberOfReactions(const State& a_state) const noexcept
{
  T Lj = std::numeric_limits<long long>::max();

  for (const auto& s : m_stateChange) {

    const size_t rI   = s.first;
    const T      nuIJ = s.second;

    if (nuIJ < 0) {
      Lj = std::min(Lj, a_state[rI] / std::abs(nuIJ));
    }
  }

  return Lj;
}

template <typename T, typename State>
inline std::list<size_t>
KMCSingleStateReaction<T, State>::getReactants() const noexcept
{
  return m_reactants;
}

template <typename T, typename State>
inline T
KMCSingleStateReaction<T, State>::getStateChange(const size_t a_particleReactant) const noexcept
{
  T nuIJ = 0;

  if (m_stateChange.find(a_particleReactant) != m_stateChange.end()) {
    nuIJ = m_stateChange.at(a_particleReactant);
  }

  return nuIJ;
}

template <typename T, typename State>
inline void
KMCSingleStateReaction<T, State>::advanceState(State& a_state, const T& a_numReactions) const noexcept
{
  for (const auto& s : m_stateChange) {
    a_state[s.first] += a_numReactions * s.second;
  }
}

#include <CD_NamespaceFooter.H>

#endif
