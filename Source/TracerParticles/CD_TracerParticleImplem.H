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
  this->m_position = RealVect::Zero;
  this->m_weight   = 0.0;
  this->m_velocity = RealVect::Zero;

  for (auto& s : this->m_scalars) {
    s = 0.0;
  }

  for (auto& v : this->m_vectors) {
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
inline int
TracerParticle<M, N>::size() const
{
  return (GenericParticle<M, N>::size() + sizeof(Real) + SpaceDim * sizeof(Real));
}

template <size_t M, size_t N>
inline void
TracerParticle<M, N>::linearOut(void* buf) const
{
  // clang-format off
  Real* buffer = (Real*) buf;

  // Linearize m_position onto buffer.
  D_TERM(*buffer++ = this->m_position[0];,
	 *buffer++ = this->m_position[1];,
	 *buffer++ = this->m_position[2];);

  // Linearize m_scalars onto the buffer
  for (size_t i = 0; i < M; i++) {
    *buffer++ = this->m_scalars[i];
  }

  // Linearize vectors onto the buffer
  for (size_t i = 0; i < N; i++) {
    const RealVect& v = this->m_vectors[i];

    D_TERM(*buffer++ = v[0];, *buffer++ = v[1];, *buffer++ = v[2];);
  }

  // Linearize m_velocity onto buffer
  D_TERM(*buffer++ = m_velocity[0];,
	 *buffer++ = m_velocity[1];,
	 *buffer++ = m_velocity[2];);

  // // Linearize weight onto buffer
  *buffer = m_weight;

  // clang-format on
}

template <size_t M, size_t N>
inline void
TracerParticle<M, N>::linearIn(void* buf)
{
  // clang-format off
  
  Real* buffer = (Real*)buf;
  
  D_TERM(this->m_position[0] = *buffer++;,
	 this->m_position[1] = *buffer++;,
	 this->m_position[2] = *buffer++;);

  // Linearize buffer onto scalars
  for (size_t i = 0; i < M; i++) {
    this->m_scalars[i] = *buffer++;
  }

  // Linearize buffer onto vectors
  for (size_t i = 0; i < N; i++) {
    RealVect& v = this->m_vectors[i];

    D_TERM(v[0] = *buffer++;,
	   v[1] = *buffer++;,
	   v[2] = *buffer++;);
  }

  // Put velocity onto buffer.
  D_TERM(m_velocity[0] = *buffer++;,
	 m_velocity[1] = *buffer++;,
	 m_velocity[2] = *buffer++;);

  // Put weight onto buffer
  m_weight = *buffer;

  // clang-format on
}

template <size_t M, size_t N>
template <size_t K>
inline Real&
TracerParticle<M, N>::real()
{
  return GenericParticle<M, N>::template real<K>();
}

template <size_t M, size_t N>
template <size_t K>
inline const Real&
TracerParticle<M, N>::real() const
{
  return GenericParticle<M, N>::template real<K>();
}

template <size_t M, size_t N>
template <size_t K>
inline RealVect&
TracerParticle<M, N>::vect()
{
  return GenericParticle<M, N>::template vect<K>();
}

template <size_t M, size_t N>
template <size_t K>
inline const RealVect&
TracerParticle<M, N>::vect() const
{
  return GenericParticle<M, N>::template vect<K>();
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
