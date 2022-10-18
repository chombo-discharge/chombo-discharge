/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCSingleStateImplem.H
  @brief  Implementation of CD_KMCSingleState.H
  @author Robert Marskar
*/

#ifndef CD_KMCSingleStateImplem_H
#define CD_KMCSingleStateImplem_H

// Our includes
#include <CD_KMCSingleState.H>
#include <CD_NamespaceHeader.H>

template <typename T>
KMCSingleState<T>::KMCSingleState(const size_t a_numSpecies) noexcept
{
  m_state.resize(a_numSpecies);
}

template <typename T>
inline KMCSingleState<T>::~KMCSingleState()
{}

template <typename T>
inline bool
KMCSingleState<T>::isValidState() const noexcept
{
  bool isValid = true;

  for (const auto& p : m_state) {
    if (p < 0) {
      isValid = false;

      break;
    }
  }

  return isValid;
}

template <typename T>
inline T&
KMCSingleState<T>::operator[](const size_t a_idx) noexcept
{
  CH_assert(a_idx < m_state.size());

  return m_state[a_idx];
}

template <typename T>
inline const T&
KMCSingleState<T>::operator[](const size_t a_idx) const noexcept
{
  CH_assert(a_idx < m_state.size());

  return m_state[a_idx];
}

template <typename T>
inline typename KMCSingleState<T>::State&
KMCSingleState<T>::getState() noexcept
{
  return m_state;
}

template <typename T>
inline const typename KMCSingleState<T>::State&
KMCSingleState<T>::getState() const noexcept
{
  return m_state;
}

#include <CD_NamespaceFooter.H>

#endif
