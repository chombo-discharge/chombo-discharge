/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCPlasmaReactionStateImplem.H
  @brief  Implementation of CD_KMCPlasmaReactionState.H
  @author Robert Marskar
*/

#ifndef CD_KMCPlasmaReactionStateImplem_H
#define CD_KMCPlasmaReactionStateImplem_H

// Our includes
#include <CD_KMCPlasmaReactionState.H>
#include <CD_NamespaceHeader.H>

template <typename T>
KMCPlasmaReactionState<T>::KMCPlasmaReactionState(const size_t a_numParticleSpecies,
                                                  const size_t a_numPhotonSpecies) noexcept
{
  m_particleState.resize(a_numParticleSpecies);
  m_photonState.resize(a_numPhotonSpecies);
}

template <typename T>
inline KMCPlasmaReactionState<T>::~KMCPlasmaReactionState()
{}

template <typename T>
inline bool
KMCPlasmaReactionState<T>::isValidState() const noexcept
{
  bool isValid = true;

  for (const auto& p : m_particleState) {
    if (p < 0) {
      isValid = false;

      break;
    }
  }

  for (const auto& p : m_photonState) {
    if (p < 0) {
      isValid = false;

      break;
    }
  }

  return isValid;
}

template <typename T>
inline KMCPlasmaReactionState<T>&
KMCPlasmaReactionState<T>::operator+=(const KMCPlasmaReactionState<T>& a_other) noexcept
{
  const auto& otherParticleState = a_other.getParticleState();
  const auto& otherPhotonState   = a_other.getPhotonState();

  const size_t numParticleSpecies = m_particleState.size();
  const size_t numPhotonSpecies   = m_photonState.size();

  CH_assert(otherParticleState.size() == m_particleState.size());
  CH_assert(otherPhotonState.size() == m_photonState.size());

  for (size_t i = 0; i < numParticleSpecies; i++) {
    m_particleState[i] += otherParticleState[i];
  }

  for (size_t i = 0; i < numPhotonSpecies; i++) {
    m_photonState[i] += otherPhotonState[i];
  }

  return *this;
}

template <typename T>
inline typename KMCPlasmaReactionState<T>::State&
KMCPlasmaReactionState<T>::getParticleState() noexcept
{
  return m_particleState;
}

template <typename T>
inline const typename KMCPlasmaReactionState<T>::State&
KMCPlasmaReactionState<T>::getParticleState() const noexcept
{
  return m_particleState;
}

template <typename T>
inline typename KMCPlasmaReactionState<T>::State&
KMCPlasmaReactionState<T>::getPhotonState() noexcept
{
  return m_photonState;
}

template <typename T>
inline const typename KMCPlasmaReactionState<T>::State&
KMCPlasmaReactionState<T>::getPhotonState() const noexcept
{
  return m_photonState;
}

#include <CD_NamespaceFooter.H>

#endif
