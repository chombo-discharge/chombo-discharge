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
  m_particleID = -1;
  m_rankID     = -1;
  m_position   = RealVect::Zero;

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
inline const int32_t&
GenericParticle<M, N>::particleID() const noexcept
{
  return m_particleID;
}

template <size_t M, size_t N>
inline int32_t&
GenericParticle<M, N>::particleID() noexcept
{
  return m_particleID;
}

template <size_t M, size_t N>
inline const int32_t&
GenericParticle<M, N>::rankID() const noexcept
{
  return m_rankID;
}

template <size_t M, size_t N>
inline int32_t&
GenericParticle<M, N>::rankID() noexcept
{
  return m_rankID;
}

template <size_t M, size_t N>
inline const std::array<Real, M>&
GenericParticle<M, N>::getReals() const noexcept
{
  return m_scalars;
}

template <size_t M, size_t N>
inline std::array<Real, M>&
GenericParticle<M, N>::getReals() noexcept
{
  return m_scalars;
}

template <size_t M, size_t N>
inline const std::array<RealVect, N>&
GenericParticle<M, N>::getVects() const noexcept
{
  return m_vectors;
}

template <size_t M, size_t N>
inline std::array<RealVect, N>&
GenericParticle<M, N>::getVects() noexcept
{
  return m_vectors;
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
inline bool
GenericParticle<M, N>::operator<(const GenericParticle<M, N>& a_other) const noexcept
{
  const RealVect& x = m_position;
  const RealVect& y = a_other.m_position;

  if (x[0] < y[0]) {
    return true;
  }
#if CH_SPACEDIM > 1
  if (x[0] > y[0]) {
    return false;
  }
  if (x[1] < y[1]) {
    return true;
  }
#endif
#if CH_SPACEDIM > 2
  if (x[1] > y[1]) {
    return false;
  }
  if (x[2] < y[2]) {
    return true;
  }
#endif

  return false;
}

template <size_t M, size_t N>
inline int
GenericParticle<M, N>::size() const
{
  // Size is the sum of the particle ID, rank ID, position, plus M Reals and N RealVects
  int size = 0;

  size += sizeof(int32_t);             // Particle ID
  size += sizeof(int32_t);             // Rank ID
  size += SpaceDim * sizeof(Real);     // Position
  size += M * sizeof(Real);            // m_scalars
  size += N * SpaceDim * sizeof(Real); // m_vectors

  return size;
}

template <size_t M, size_t N>
inline void
GenericParticle<M, N>::linearOut(void* const a_buffer) const
{
  uint8_t* p = static_cast<uint8_t*>(a_buffer);

  detail::pushParticleProperty(p, m_particleID);
  detail::pushParticleProperty(p, m_rankID);

  detail::pushParticleProperty(p, m_position[0]);
  detail::pushParticleProperty(p, m_position[1]);
#if CH_SPACEDIM == 3
  detail::pushParticleProperty(p, m_position[2]);
#endif

  for (size_t i = 0; i < M; i++) {
    detail::pushParticleProperty(p, m_scalars[i]);
  }

  for (size_t i = 0; i < N; i++) {
    detail::pushParticleProperty(p, m_vectors[i][0]);
    detail::pushParticleProperty(p, m_vectors[i][1]);
#if CH_SPACEDIM == 3
    detail::pushParticleProperty(p, m_vectors[i][2]);
#endif
  }
}

template <size_t M, size_t N>
inline void
GenericParticle<M, N>::linearIn(const void* const a_buffer)
{
  const uint8_t* p = static_cast<const uint8_t*>(a_buffer);

  detail::pullParticleProperty(p, m_particleID);
  detail::pullParticleProperty(p, m_rankID);

  detail::pullParticleProperty(p, m_position[0]);
  detail::pullParticleProperty(p, m_position[1]);
#if CH_SPACEDIM == 3
  detail::pullParticleProperty(p, m_position[2]);
#endif

  for (size_t i = 0; i < M; i++) {
    detail::pullParticleProperty(p, m_scalars[i]);
  }

  for (size_t i = 0; i < N; i++) {
    detail::pullParticleProperty(p, m_vectors[i][0]);
    detail::pullParticleProperty(p, m_vectors[i][1]);
#if CH_SPACEDIM == 3
    detail::pullParticleProperty(p, m_vectors[i][2]);
#endif
  }
}

template <size_t M, size_t N>
inline int
GenericParticle<M, N>::H5size() const
{
  int size = 0;

  size += SpaceDim * sizeof(Real);     // m_position
  size += M * sizeof(Real);            // m_scalars
  size += N * SpaceDim * sizeof(Real); // m_vectors

  return size;
}

template <size_t M, size_t N>
inline void
GenericParticle<M, N>::H5linearOut(void* const a_buffer) const
{
  uint8_t* p = static_cast<uint8_t*>(a_buffer);

  detail::pushParticleProperty(p, m_position[0]);
  detail::pushParticleProperty(p, m_position[1]);
#if CH_SPACEDIM == 3
  detail::pushParticleProperty(p, m_position[2]);
#endif

  for (size_t i = 0; i < M; i++) {
    detail::pushParticleProperty(p, m_scalars[i]);
  }

  for (size_t i = 0; i < N; i++) {
    detail::pushParticleProperty(p, m_vectors[i][0]);
    detail::pushParticleProperty(p, m_vectors[i][1]);
#if CH_SPACEDIM == 3
    detail::pushParticleProperty(p, m_vectors[i][2]);
#endif
  }
}

template <size_t M, size_t N>
inline void
GenericParticle<M, N>::H5linearIn(const void* const a_buffer)
{
  const uint8_t* p = static_cast<const uint8_t*>(a_buffer);

  detail::pullParticleProperty(p, m_position[0]);
  detail::pullParticleProperty(p, m_position[1]);
#if CH_SPACEDIM == 3
  detail::pullParticleProperty(p, m_position[2]);
#endif

  for (size_t i = 0; i < M; i++) {
    detail::pullParticleProperty(p, m_scalars[i]);
  }

  for (size_t i = 0; i < N; i++) {
    detail::pullParticleProperty(p, m_vectors[i][0]);
    detail::pullParticleProperty(p, m_vectors[i][1]);
#if CH_SPACEDIM == 3
    detail::pullParticleProperty(p, m_vectors[i][2]);
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
