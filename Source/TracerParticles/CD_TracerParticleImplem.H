/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TracerParticleImplem.H
  @brief  Implementation of CD_TracerParticle.H
  @author Robert Marskar
*/

#ifndef CD_TracerParticleImplem_H
#define CD_TracerParticleImplem_H

// Our includes
#include <CD_TracerParticle.H>
#include <CD_NamespaceHeader.H>

template <size_t M, size_t N>
inline TracerParticle<M, N>::TracerParticle()
{

  m_position = RealVect::Zero;
  m_weight   = 0.0;
  m_velocity = RealVect::Zero;

  for (auto& s : m_scalars) {
    s = 0.0;
  }

  for (auto& v : m_vectors) {
    v = RealVect::Zero;
  }
}

template <size_t M, size_t N>
inline TracerParticle<M, N>::TracerParticle(const TracerParticle<M, N>& a_other)
{
  this->m_position = a_other.m_position;
  this->m_weight   = a_other.m_weight;
  this->m_velocity = a_other.m_velocity;
  this->m_scalars  = a_other.m_scalars;
  this->m_vectors  = a_other.m_vectors;
}

template <size_t M, size_t N>
inline TracerParticle<M, N>::TracerParticle(const RealVect& a_position) : TracerParticle<M, N>()
{
  this->m_position = a_position;
}

template <size_t M, size_t N>
inline TracerParticle<M, N>::~TracerParticle()
{}

template <size_t M, size_t N>
inline Real&
TracerParticle<M, N>::weight()
{
  return m_weight;
}

template <size_t M, size_t N>
inline const Real&
TracerParticle<M, N>::weight() const
{
  return m_weight;
}

template <size_t M, size_t N>
inline RealVect&
TracerParticle<M, N>::position()
{
  return m_position;
}

template <size_t M, size_t N>
inline const RealVect&
TracerParticle<M, N>::position() const
{
  return m_position;
}

template <size_t M, size_t N>
inline RealVect&
TracerParticle<M, N>::velocity()
{
  return m_velocity;
}

template <size_t M, size_t N>
inline const RealVect&
TracerParticle<M, N>::velocity() const
{
  return m_velocity;
}

template <size_t M, size_t N>
template <size_t K>
inline Real&
TracerParticle<M, N>::scalar()
{
  return std::get<K>(m_scalars);
}

template <size_t M, size_t N>
template <size_t K>
inline const Real&
TracerParticle<M, N>::scalar() const
{
  return std::get<K>(m_scalars);
}

template <size_t M, size_t N>
template <size_t K>
inline RealVect&
TracerParticle<M, N>::vect()
{
  return std::get<K>(m_vectors);
}

template <size_t M, size_t N>
template <size_t K>
inline const RealVect&
TracerParticle<M, N>::vect() const
{
  return std::get<K>(m_vectors);
}

template <size_t M, size_t N>
inline bool
TracerParticle<M, N>::operator==(const TracerParticle<M, N>& a_p) const
{
  return (this->m_position == a_p.m_position) && (this->m_weight == a_p.m_weight)
         && (this->m_velocity == a_p.m_velocity);
}

template <size_t M, size_t N>
inline bool
TracerParticle<M, N>::operator!=(const TracerParticle<M, N>& a_p) const
{
  return !(*this == a_p);
}

template <size_t M, size_t N>
inline size_t
TracerParticle<M, N>::size() const
{
  // clang-format off
  //     position..............    weight.....    velocity..............    scalars........    vectors...................
  return SpaceDim * sizeof(Real) + sizeof(Real) + SpaceDim * sizeof(Real) + M * sizeof(Real) + N * SpaceDim * sizeof(Real);
  // clang-format on
}

template <size_t M, size_t N>
inline void
TracerParticle<M, N>::linearOut(void* buf) const
{
  Real* buffer = (Real*)buf;

  // Linearize m_position onto buffer.
  *buffer++ = m_position[0];
  *buffer++ = m_position[1];
#if CH_SPACEDIM == 3
  *buffer++ = m_position[2];
#endif

  // Linearize m_velocity onto buffer
  *buffer++ = m_velocity[0];
  *buffer++ = m_velocity[1];
#if CH_SPACEDIM == 3
  *buffer++ = m_velocity[2];
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

  // Linearize weight onto buffer
  *buffer = m_weight;
}

template <size_t M, size_t N>
inline void
TracerParticle<M, N>::linearIn(void* buf)
{
  Real* buffer = (Real*)buf;

  m_position[0] = *buffer++;
  m_position[1] = *buffer++;
#if CH_SPACEDIM == 3
  m_position[2] = *buffer++;
#endif

  // Put velocity onto buffer.
  m_velocity[0] = *buffer++;
  m_velocity[1] = *buffer++;
#if CH_SPACEDIM == 3
  m_velocity[2] = *buffer++;
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

  // Put weight onto buffer
  m_weight = *buffer;
}

template <size_t M, size_t N>
inline std::ostream&
operator<<(std::ostream& ostr, const TracerParticle<M, N>& p)
{
  ostr << "TracerParticle : \n";

  // Print position.
  ostr << "\tPosition  = " << p.position() << "\n";
  ostr << "\tWeight    = " << p.weight() << "\n";
  ostr << "\tVelocity  = " << p.velocity() << "\n";

  return ostr;
}

#include <CD_NamespaceFooter.H>

#endif
