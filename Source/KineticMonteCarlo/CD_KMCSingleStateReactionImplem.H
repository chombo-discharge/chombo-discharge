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

template <typename State, typename T>
inline KMCSingleStateReaction<State, T>::KMCSingleStateReaction(const std::list<size_t>& a_reactants,
                                                                const std::list<size_t>& a_products) noexcept
{
  m_reactants = a_reactants;
  m_products  = a_products;

  this->computeStateChanges();
}

template <typename State, typename T>
inline KMCSingleStateReaction<State, T>::~KMCSingleStateReaction()
{}

template <typename State, typename T>
inline void
KMCSingleStateReaction<State, T>::computeStateChanges() noexcept
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
  for (const auto& r : m_reactants) {
    if (reactantNumbers.find(r) == reactantNumbers.end()) {
      reactantNumbers.emplace(r, 1);
    }
    else {
      reactantNumbers[r]++;
    }
  }

  for (const auto& rn : reactantNumbers) {
    m_propensityFactor *= 1.0 / factorial(rn.second);
  }
}

template <typename State, typename T>
inline Real&
KMCSingleStateReaction<State, T>::rate() const noexcept
{
  return m_rate;
}

template <typename State, typename T>
inline T
KMCSingleStateReaction<State, T>::population(const size_t& a_reactant, const State& a_state) const noexcept
{
  return a_state[a_reactant];
}

template <typename State, typename T>
inline Real
KMCSingleStateReaction<State, T>::propensity(const State& a_state) const noexcept
{
  Real A = m_rate * m_propensityFactor;

  State tmp = a_state;

  for (const auto& r : m_reactants) {
    A *= tmp[r];

    tmp[r]--;
  }

  return A;
}

template <typename State, typename T>
inline T
KMCSingleStateReaction<State, T>::computeCriticalNumberOfReactions(const State& a_state) const noexcept
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

template <typename State, typename T>
inline std::list<size_t>
KMCSingleStateReaction<State, T>::getReactants() const noexcept
{
  return m_reactants;
}

template <typename State, typename T>
inline T
KMCSingleStateReaction<State, T>::getStateChange(const size_t a_particleReactant) const noexcept
{
  T nuIJ = 0;

  if (m_stateChange.find(a_particleReactant) != m_stateChange.end()) {
    nuIJ = m_stateChange.at(a_particleReactant);
  }

  return nuIJ;
}

template <typename State, typename T>
inline void
KMCSingleStateReaction<State, T>::advanceState(State& a_state, const T& a_numReactions) const noexcept
{
  for (const auto& s : m_stateChange) {
    a_state[s.first] += a_numReactions * s.second;
  }
}

#include <CD_NamespaceFooter.H>

#endif
