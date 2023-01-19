/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrIteratorImplem.H
  @brief  Implementation of CdrIterator.H
  @author Robert Marskar
*/

#ifndef CD_CdrIteratorImplem_H
#define CD_CdrIteratorImplem_H

// Our includes
#include <CD_CdrIterator.H>
#include <CD_NamespaceHeader.H>

template <class T>
CdrIterator<T>::CdrIterator(CdrLayout<T>& a_layout, const SpeciesIteration a_mode)
{
  CH_TIME("CdrIterator<T>::CdrIterator(CdrLayout<T>, SpeciesIteration)");

  // Default settings
  m_solvers    = a_layout.getSolvers();
  m_species    = a_layout.getSpecies();
  m_numSolvers = m_solvers.size();
  m_mode       = a_mode;

  m_isDefined = true;

  reset();
}

template <class T>
CdrIterator<T>::~CdrIterator()
{
  CH_TIME("CdrIterator<T>::~CdrIterator()");
}

template <class T>
int
CdrIterator<T>::getNumberOfSolvers()
{
  CH_TIME("CdrIterator<T>::getNumberOfSolvers()");

  CH_assert(m_isDefined);

  return m_numSolvers;
}

template <class T>
int
CdrIterator<T>::index() const
{
  CH_TIME("CdrIterator<T>::index()");

  CH_assert(m_isDefined);

  return m_isolver;
}

template <class T>
void
CdrIterator<T>::reset()
{
  CH_TIME("CdrIterator<T>::reset()");

  CH_assert(m_isDefined);

  m_isolver = 0;
}

template <class T>
bool
CdrIterator<T>::ok() const
{
  CH_TIME("CdrIterator<T>::ok()");

  CH_assert(m_isDefined);

  return (m_isolver < m_numSolvers);
}

template <class T>
void
CdrIterator<T>::operator++()
{
  CH_TIME("CdrIterator<T>::operator++()");

  CH_assert(m_isDefined);

  if (m_mode == SpeciesIteration::All) {
    m_isolver++;
  }
  else {
    m_isolver++;

    for (int isolver = m_isolver; this->ok(); ++isolver) {
      const RefCountedPtr<CdrSpecies>& species = m_species[isolver];

      if (m_mode == SpeciesIteration::Charged) {
        if (species->getChargeNumber() != 0) {
          m_isolver = isolver;
          break;
        }
      }
      else if (m_mode == SpeciesIteration::Negative) {
        if (species->getChargeNumber() < 0) {
          m_isolver = isolver;
          break;
        }
      }
      else if (m_mode == SpeciesIteration::Positive) {
        if (species->getChargeNumber() > 0) {
          m_isolver = isolver;
          break;
        }
      }
      else if (m_mode == SpeciesIteration::Neutral) {
        if (species->getChargeNumber() == 0) {
          m_isolver = isolver;
          break;
        }
      }
    }
  }
}

template <class T>
RefCountedPtr<T>&
CdrIterator<T>::operator()()
{
  CH_TIME("CdrIterator<T>::operator()");

  CH_assert(m_isDefined);
  CH_assert(m_isolver < m_numSolvers);

  return m_solvers[m_isolver];
}

template <class T>
RefCountedPtr<CdrSpecies>&
CdrIterator<T>::getSpecies()
{
  CH_TIME("CdrIterator<T>::getSpecies()");

  CH_assert(m_isDefined);
  CH_assert(m_isolver < m_numSolvers);

  return m_species[m_isolver];
}

#include <CD_NamespaceFooter.H>

#endif
