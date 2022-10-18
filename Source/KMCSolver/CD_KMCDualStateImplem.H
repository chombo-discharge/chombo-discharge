/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCDualStateImplem.H
  @brief  Implementation of CD_KMCDualState.H
  @author Robert Marskar
*/

#ifndef CD_KMCDualStateImplem_H
#define CD_KMCDualStateImplem_H

// Our includes
#include <CD_KMCDualState.H>
#include <CD_NamespaceHeader.H>

template <typename T>
KMCDualState<T>::KMCDualState(const size_t a_numReactiveSpecies, const size_t a_numNonReactiveSpecies) noexcept
{
  m_reactiveState.resize(a_numReactiveSpecies);
  m_nonReactiveState.resize(a_numNonReactiveSpecies);
}

template <typename T>
inline KMCDualState<T>::~KMCDualState()
{}

template <typename T>
inline bool
KMCDualState<T>::isValidState() const noexcept
{
  bool isValid = true;

  for (const auto& p : m_reactiveState) {
    if (p < 0) {
      isValid = false;

      break;
    }
  }

  for (const auto& p : m_nonReactiveState) {
    if (p < 0) {
      isValid = false;

      break;
    }
  }

  return isValid;
}

template <typename T>
inline typename KMCDualState<T>::State&
KMCDualState<T>::getReactiveState() noexcept
{
  return m_reactiveState;
}

template <typename T>
inline const typename KMCDualState<T>::State&
KMCDualState<T>::getReactiveState() const noexcept
{
  return m_reactiveState;
}

template <typename T>
inline typename KMCDualState<T>::State&
KMCDualState<T>::getNonReactiveState() noexcept
{
  return m_nonReactiveState;
}

template <typename T>
inline const typename KMCDualState<T>::State&
KMCDualState<T>::getNonReactiveState() const noexcept
{
  return m_nonReactiveState;
}

#include <CD_NamespaceFooter.H>

#endif
