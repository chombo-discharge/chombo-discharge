/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_NonCommParticleImplem.H
  @brief  Implementation of CD_NonCommParticle.H
  @author Robert Marskar
*/

#ifndef CD_NONCOMMPARTICLEIMPLEM_H
#define CD_NONCOMMPARTICLEIMPLEM_H

// Our includes
#include <CD_NonCommParticle.H>
#include <CD_NamespaceHeader.H>

template <size_t M, size_t N>
inline NonCommParticle<M, N>::NonCommParticle()
{
  for (auto& s : m_scalars) {
    s = 0.0;
  }

  for (auto& v : m_vectors) {
    v = RealVect::Zero;
  }
}

template <size_t M, size_t N>
inline NonCommParticle<M, N>::NonCommParticle(const NonCommParticle<M, N>& a_other)
{
  this->m_scalars = a_other.m_scalars;
  this->m_vectors = a_other.m_vectors;
}

template <size_t M, size_t N>
inline NonCommParticle<M, N>::~NonCommParticle()
{}

template <size_t M, size_t N>
template <size_t K>
inline Real&
NonCommParticle<M, N>::real()
{
  return std::get<K>(m_scalars);
}

template <size_t M, size_t N>
template <size_t K>
inline const Real&
NonCommParticle<M, N>::real() const
{
  return std::get<K>(m_scalars);
}

template <size_t M, size_t N>
template <size_t K>
inline RealVect&
NonCommParticle<M, N>::vect()
{
  return std::get<K>(m_vectors);
}

template <size_t M, size_t N>
template <size_t K>
inline const RealVect&
NonCommParticle<M, N>::vect() const
{
  return std::get<K>(m_vectors);
}

template <size_t M, size_t N>
inline bool
NonCommParticle<M, N>::operator==(const NonCommParticle<M, N>& a_other) const
{
  return (m_scalars = a_other.m_scalars) && (m_vectors == a_other.m_vectors);
}

template <size_t M, size_t N>
inline bool
NonCommParticle<M, N>::operator!=(const NonCommParticle<M, N>& a_other) const
{
  return !(*this == a_other);
}

#include <CD_NamespaceFooter.H>

#endif
