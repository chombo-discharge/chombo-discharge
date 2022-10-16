/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PlasmaReactionImplem.H
  @brief  Implementation of CD_PlasmaReaction.H
  @author Robert Marskar
*/

#ifndef CD_PlasmaReactionImplem_H
#define CD_PlasmaReactionImplem_H

// Our includes
#include <CD_PlasmaReaction.H>
#include <CD_NamespaceHeader.H>

template <typename T, typename State>
inline PlasmaReaction<T, State>::PlasmaReaction(const std::list<size_t>& a_particleReactants,
                                                const std::list<size_t>& a_particleProducts,
                                                const std::list<size_t>& a_photonProducts) noexcept
{
  m_particleReactants = a_particleReactants;
  m_particleProducts  = a_particleProducts;
  m_photonProducts    = a_photonProducts;

  this->computeStateChanges();
}

template <typename T, typename State>
inline PlasmaReaction<T, State>::~PlasmaReaction()
{}

template <typename T, typename State>
inline void
PlasmaReaction<T, State>::computeStateChanges() noexcept
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
}

template <typename T, typename State>
inline Real&
PlasmaReaction<T, State>::rate() const noexcept
{
  return m_rate;
}

template <typename T, typename State>
inline Real
PlasmaReaction<T, State>::propensity(const State& a_particleState) const noexcept
{
  this->sanityCheckParticles(a_particleState);

  Real prop = m_rate;

  for (const auto& idx : m_particleReactants) {
    prop *= a_particleState[idx];
  }

  return prop;
}

template <typename T, typename State>
inline void
PlasmaReaction<T, State>::advanceState(State&   a_particleState,
                                       State&   a_photonState,
                                       const T& a_numReactions) const noexcept
{
  this->sanityCheckParticles(a_particleState);
  this->sanityCheckPhotons(a_photonState);

  // Particle state consumption/gain.
  for (const auto& change : m_particleStateChange) {
    a_particleState[change.first] += a_numReactions * change.second;
  }

  // Photon gain.
  for (const auto& change : m_photonStateChange) {
    a_photonState[change.first] += a_numReactions * change.second;
  }
}

template <typename T, typename State>
inline const std::map<size_t, int>&
PlasmaReaction<T, State>::getParticleStateChange() const noexcept
{
  return m_particleStateChange;
}

template <typename T, typename State>
inline const std::map<size_t, int>&
PlasmaReaction<T, State>::getPhotonStateChange() const noexcept
{
  return m_photonStateChange;
}

template <typename T, typename State>
inline const std::list<size_t>&
PlasmaReaction<T, State>::getParticleReactants() const noexcept
{
  return m_particleReactants;
}

template <typename T, typename State>
inline const std::list<size_t>&
PlasmaReaction<T, State>::getParticleProducts() const noexcept
{
  return m_particleProducts;
}

template <typename T, typename State>
inline const std::list<size_t>&
PlasmaReaction<T, State>::getPhotonProducts() const noexcept
{
  return m_photonProducts;
}

template <typename T, typename State>
inline void
PlasmaReaction<T, State>::sanityCheckParticles(const State& a_particleState) const noexcept
{
  for (const auto& idx : m_particleReactants) {
    CH_assert(a_particleState.size() > idx);
  }

  for (const auto& idx : m_particleProducts) {
    CH_assert(a_particleState.size() > idx);
  }
}

template <typename T, typename State>
inline void
PlasmaReaction<T, State>::sanityCheckPhotons(const State& a_photonState) const noexcept
{
  for (const auto& idx : m_photonProducts) {
    CH_assert(a_photonState.size() > idx);
  }
}

#include <CD_NamespaceFooter.H>

#endif
