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
  int size = 0;

  size += GenericParticle<M, N>::size();
  size += SpaceDim * sizeof(Real);
  size += sizeof(Real);

  return size;
}

template <size_t M, size_t N>
inline void
TracerParticle<M, N>::linearOut(void* const a_buffer) const
{
  uint8_t* p = static_cast<uint8_t*>(a_buffer);

  detail::pushParticleProperty(p, this->m_particleID);
  detail::pushParticleProperty(p, this->m_rankID);

  detail::pushParticleProperty(p, this->m_position[0]);
  detail::pushParticleProperty(p, this->m_position[1]);
#if CH_SPACEDIM == 3
  detail::pushParticleProperty(p, this->m_position[2]);
#endif

  for (size_t i = 0; i < M; i++) {
    detail::pushParticleProperty(p, this->m_scalars[i]);
  }

  for (size_t i = 0; i < M; i++) {
    detail::pushParticleProperty(p, this->m_vectors[i][0]);
    detail::pushParticleProperty(p, this->m_vectors[i][1]);
#if CH_SPACEDIM == 3
    detail::pushParticleProperty(p, this->m_vectors[i][2]);
#endif
  }

  detail::pushParticleProperty(p, m_weight);

  detail::pushParticleProperty(p, m_velocity[0]);
  detail::pushParticleProperty(p, m_velocity[1]);
#if CH_SPACEDIM == 3
  detail::pushParticleProperty(p, m_velocity[2]);
#endif
}

template <size_t M, size_t N>
inline void
TracerParticle<M, N>::linearIn(const void* const a_buffer)
{
  const uint8_t* p = static_cast<const uint8_t*>(a_buffer);

  detail::pullParticleProperty(p, this->m_particleID);
  detail::pullParticleProperty(p, this->m_rankID);

  detail::pullParticleProperty(p, this->m_position[0]);
  detail::pullParticleProperty(p, this->m_position[1]);
#if CH_SPACEDIM == 3
  detail::pullParticleProperty(p, this->m_position[2]);
#endif

  for (size_t i = 0; i < M; i++) {
    detail::pullParticleProperty(p, this->m_scalars[i]);
  }

  for (size_t i = 0; i < M; i++) {
    detail::pullParticleProperty(p, this->m_vectors[i][0]);
    detail::pullParticleProperty(p, this->m_vectors[i][1]);
#if CH_SPACEDIM == 3
    detail::pullParticleProperty(p, this->m_vectors[i][2]);
#endif
  }

  detail::pullParticleProperty(p, m_weight);

  detail::pullParticleProperty(p, m_velocity[0]);
  detail::pullParticleProperty(p, m_velocity[1]);
#if CH_SPACEDIM == 3
  detail::pullParticleProperty(p, m_velocity[2]);
#endif
}

template <size_t M, size_t N>
inline int
TracerParticle<M, N>::H5size() const
{
  int size = 0;

  size += SpaceDim * sizeof(Real);     // Position.
  size += M * sizeof(Real);            // m_scalars
  size += N * SpaceDim * sizeof(Real); // m_vectors
  size += sizeof(Real);                // m_weight
  size += SpaceDim * sizeof(Real);     // m_velocity

  return size;
}

template <size_t M, size_t N>
inline void
TracerParticle<M, N>::H5linearOut(void* const a_buffer) const
{
  uint8_t* p = static_cast<uint8_t*>(a_buffer);

  detail::pushParticleProperty(p, this->m_position[0]);
  detail::pushParticleProperty(p, this->m_position[1]);
#if CH_SPACEDIM == 3
  detail::pushParticleProperty(p, this->m_position[2]);
#endif

  for (size_t i = 0; i < M; i++) {
    detail::pushParticleProperty(p, this->m_scalars[i]);
  }

  for (size_t i = 0; i < M; i++) {
    detail::pushParticleProperty(p, this->m_vectors[i][0]);
    detail::pushParticleProperty(p, this->m_vectors[i][1]);
#if CH_SPACEDIM == 3
    detail::pushParticleProperty(p, this->m_vectors[i][2]);
#endif
  }

  detail::pushParticleProperty(p, m_weight);

  detail::pushParticleProperty(p, m_velocity[0]);
  detail::pushParticleProperty(p, m_velocity[1]);
#if CH_SPACEDIM == 3
  detail::pushParticleProperty(p, m_velocity[2]);
#endif
}

template <size_t M, size_t N>
inline void
TracerParticle<M, N>::H5linearIn(const void* const a_buffer)
{
  const uint8_t* p = static_cast<const uint8_t*>(a_buffer);

  detail::pullParticleProperty(p, this->m_position[0]);
  detail::pullParticleProperty(p, this->m_position[1]);
#if CH_SPACEDIM == 3
  detail::pullParticleProperty(p, this->m_position[2]);
#endif

  for (size_t i = 0; i < M; i++) {
    detail::pullParticleProperty(p, this->m_scalars[i]);
  }

  for (size_t i = 0; i < M; i++) {
    detail::pullParticleProperty(p, this->m_vectors[i][0]);
    detail::pullParticleProperty(p, this->m_vectors[i][1]);
#if CH_SPACEDIM == 3
    detail::pullParticleProperty(p, this->m_vectors[i][2]);
#endif
  }

  detail::pullParticleProperty(p, m_weight);

  detail::pullParticleProperty(p, m_velocity[0]);
  detail::pullParticleProperty(p, m_velocity[1]);
#if CH_SPACEDIM == 3
  detail::pullParticleProperty(p, m_velocity[2]);
#endif
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
