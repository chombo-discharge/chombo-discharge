/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_GenericParticleImplem.H
  @brief  Implementation of CD_GenericParticle.H
  @author Robert Marskar
*/

#ifndef CD_GenericParticleImplem_H
#define CD_GenericParticleImplem_H

// Our includes
#include <CD_GenericParticle.H>
#include <CD_NamespaceHeader.H>

template <size_t M, size_t N>
inline GenericParticle<M, N>::GenericParticle()
{
  m_position = RealVect::Zero;

  for (auto& s : m_scalars) {
    s = 0.0;
  }

  for (auto& v : m_vectors) {
    v = RealVect::Zero;
  }
}

template <size_t M, size_t N>
inline GenericParticle<M, N>::GenericParticle(const GenericParticle<M, N>& a_other)
{
  this->m_position = a_other.m_position;
  this->m_scalars  = a_other.m_scalars;
  this->m_vectors  = a_other.m_vectors;
}

template <size_t M, size_t N>
inline GenericParticle<M, N>::~GenericParticle()
{}

template <size_t M, size_t N>
inline RealVect&
GenericParticle<M, N>::position()
{
  return m_position;
}

template <size_t M, size_t N>
inline const RealVect&
GenericParticle<M, N>::position() const
{
  return m_position;
}

template <size_t M, size_t N>
template <size_t K>
inline Real&
GenericParticle<M, N>::real()
{
  return std::get<K>(m_scalars);
}

template <size_t M, size_t N>
template <size_t K>
inline const Real&
GenericParticle<M, N>::real() const
{
  return std::get<K>(m_scalars);
}

template <size_t M, size_t N>
template <size_t K>
inline RealVect&
GenericParticle<M, N>::vect()
{
  return std::get<K>(m_vectors);
}

template <size_t M, size_t N>
template <size_t K>
inline const RealVect&
GenericParticle<M, N>::vect() const
{
  return std::get<K>(m_vectors);
}

template <size_t M, size_t N>
inline bool
GenericParticle<M, N>::operator==(const GenericParticle<M, N>& a_p) const
{
  return (this->m_position == a_p.m_position && this->m_scalars == a_p.m_scalars && this->m_vectors == a_p.m_vectors);
}

template <size_t M, size_t N>
inline bool
GenericParticle<M, N>::operator!=(const GenericParticle<M, N>& a_p) const
{
  return !(*this == a_p);
}

template <size_t M, size_t N>
inline int
GenericParticle<M, N>::size() const
{
  // Size is the position plus M Reals and N RealVects
  return SpaceDim * sizeof(Real) + M * sizeof(Real) + N * SpaceDim * sizeof(Real);
}

template <size_t M, size_t N>
inline void
GenericParticle<M, N>::linearOut(void* buf) const
{
  Real* buffer = (Real*)buf;

  // Linearize m_position onto buffer.
  *buffer++ = m_position[0];
  *buffer++ = m_position[1];
#if CH_SPACEDIM == 3
  *buffer++ = m_position[2];
#endif

  // Linearize m_scalars onto the buffer
  for (size_t i = 0; i < M; i++) {
    *buffer++ = m_scalars[i];
  }

  // Linearize vectors onto the buffer
  for (size_t i = 0; i < N; i++) {
    const RealVect& v = m_vectors[i];

    *buffer++ = v[0];
    *buffer++ = v[1];
#if CH_SPACEDIM == 3
    *buffer++ = v[2];
#endif
  }
}

template <size_t M, size_t N>
inline void
GenericParticle<M, N>::linearIn(void* buf)
{
  Real* buffer = (Real*)buf;

  m_position[0] = *buffer++;
  m_position[1] = *buffer++;
#if CH_SPACEDIM == 3
  m_position[2] = *buffer++;
#endif

  // Linearize buffer onto scalars
  for (size_t i = 0; i < M; i++) {
    m_scalars[i] = *buffer++;
  }

  // Linearize buffer onto vectors
  for (size_t i = 0; i < N; i++) {
    RealVect& v = m_vectors[i];

    v[0] = *buffer++;
    v[1] = *buffer++;
#if CH_SPACEDIM == 3
    v[2] = *buffer++;
#endif
  }
}

template <size_t M, size_t N>
inline std::ostream&
operator<<(std::ostream& ostr, const GenericParticle<M, N>& p)
{
  ostr << "GenericParticle : \n";

  // Print position.
  ostr << "\tPosition  = " << p.position() << "\n";

  return ostr;
}

#include <CD_NamespaceFooter.H>

#endif
