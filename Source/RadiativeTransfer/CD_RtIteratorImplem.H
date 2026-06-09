/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RtIteratorImplem.H
  @brief  Implementation of RtIterator.H
  @author Robert Marskar
*/

#ifndef CD_RTITERATORIMPLEM_H
#define CD_RTITERATORIMPLEM_H

// Our includes
#include <CD_RtIterator.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Iterator class for RtLayout
  @param[in,out] a_layout Description
*/
template <class T>
RtIterator<T>::RtIterator(RtLayout<T>& a_layout)
{
  CH_TIME("RtIterator::RtIterator");

  m_solvers = a_layout.getSolvers();
  m_species = a_layout.getSpecies();
  m_num     = m_solvers.size();

  this->reset();
}

template <class T>
RtIterator<T>::~RtIterator()
{
  CH_TIME("RtIterator::~RtIterator");
}

template <class T>
int
RtIterator<T>::getNumberOfSolvers()
{
  return m_num;
}

template <class T>
int
RtIterator<T>::index() const
{
  return m_isolver;
}

template <class T>
void
RtIterator<T>::reset()
{
  m_isolver = 0;
}

template <class T>
bool
RtIterator<T>::ok()
{
  return (m_isolver < m_num);
}

template <class T>
void
RtIterator<T>::operator++()
{
  m_isolver++;
}

template <class T>
RefCountedPtr<T>&
RtIterator<T>::operator()()
{
  CH_assert(m_isolver < m_num);
  return m_solvers[m_isolver];
}
template <class T>
RefCountedPtr<RtSpecies>&
RtIterator<T>::getSpecies()
{
  CH_assert(m_isolver < m_num);
  return m_species[m_isolver];
}

#include <CD_NamespaceFooter.H>

#endif
