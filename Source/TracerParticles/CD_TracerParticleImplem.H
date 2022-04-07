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
inline
TracerParticle<M,N>::TracerParticle() {
  for (auto& s : m_scalars) {
    s = 0.0;
  }

  for (auto& v : m_vectors) {
    v = RealVect::Zero;
  }

  m_position = RealVect::Zero;
  m_process  = false;
}

template <size_t M, size_t N>
inline  
TracerParticle<M, N>::TracerParticle(const TracerParticle<M, N>& a_other) {
  this->m_position = a_other.m_position;
  this->m_scalars  = a_other.m_scalars;
  this->m_vectors  = a_other.m_vectors;    
}

template <size_t M, size_t N>
inline
TracerParticle<M,N>::~TracerParticle() {
}

template <size_t M, size_t N>
inline
bool& TracerParticle<M, N>::processParticle() {
  return m_process;
}

template <size_t M, size_t N>
inline
const bool& TracerParticle<M, N>::processParticle() const {
  return m_process;
}

template<size_t M, size_t N>
template<size_t which>
inline
Real& TracerParticle<M, N>::getScalar() {
  return m_scalars[which];
}

template<size_t M, size_t N>
template<size_t which>
inline
const Real& TracerParticle<M, N>::getScalar() const {
  return m_scalars[which];
}

template <size_t M, size_t N>
inline  
Real& TracerParticle<M, N>::getScalar(const size_t a_which) {
  return m_scalars[a_which];
}

template <size_t M, size_t N>
inline  
const Real& TracerParticle<M, N>::getScalar(const size_t a_which) const {
  return m_scalars[a_which];
}

template <size_t M, size_t N>
inline  
RealVect& TracerParticle<M, N>::getVector(const size_t a_which) {
  return m_vectors[a_which];
}

template <size_t M, size_t N>
inline
const RealVect& TracerParticle<M, N>::getVector(const size_t a_which) const {
  return m_vectors[a_which];
}

template <size_t M, size_t N>
inline
bool TracerParticle<M, N>::operator==(const TracerParticle<M, N>& a_p) const{
  return ( this->m_position == a_p.m_position &&
	   this->m_scalars  == a_p.m_scalars  &&
	   this->m_vectors  == a_p.m_vectors);
}

template <size_t M, size_t N>
inline
bool TracerParticle<M, N>::operator!=(const TracerParticle<M, N>& a_p) const{
  return !(*this == a_p);
}

template <size_t M, size_t N>
inline
int TracerParticle<M, N>::size() const{
  return BinItem::size() + sizeof(m_scalars) + sizeof(m_vectors) + sizeof(m_process);
}

template <size_t M, size_t N>
inline
void TracerParticle<M, N>::linearOut(void* buf) const{
  Real* buffer = (Real*)buf;

  // Linearize position onto buffer. 
  *buffer++ = m_position[0];
  *buffer++ = m_position[1];
#if CH_SPACEDIM==3
  *buffer++ = m_position[2];
#endif

  // Linearize scalars onto buffer. 
  for (size_t i = 0; i < M; i++) {
    *buffer++ = m_scalars[i];
  }  

  // Linearize vectors onto buffer. 
  for (size_t i = 0; i < N; i++) {
    const RealVect& v = m_vectors[i];

    *buffer++ = v[0];
    *buffer++ = v[1];
#if CH_SPACEDIM==3
    *buffer++ = v[2];
#endif
  }

  // Linearize the processing flag onto the buffer. 
  bool* buf2 = (bool*) buf;
  *buf2++ = m_process;
}

template <size_t M, size_t N>
inline
void TracerParticle<M, N>::linearIn(void* buf) {
  Real* buffer = (Real*) buf;

  m_position[0] = *buffer++;
  m_position[1] = *buffer++;
#if CH_SPACEDIM==3
  m_position[2] = *buffer++;
#endif

  for (size_t i = 0; i < M; i++) {
    m_scalars[i] = *buffer++;
  }

  for (size_t i = 0; i < N; i++) {
    RealVect& v = m_vectors[i];

    v[0] = *buffer++;
    v[1] = *buffer++;
#if CH_SPACEDIM==3
    v[2] = *buffer++;
#endif
  }

  // Linearize the processing flag onto the buffer. 
  bool* buf2 = (bool*) buf;
  m_process = *buf2++;
}

template <size_t M, size_t N>
inline
std::ostream & operator<<(std::ostream& ostr, const TracerParticle<M, N>& p) {
  ostr << "TracerParticle : \n";

  // Print position. 
  ostr << "\tPosition  = " << p.position() << "\n";

  // Print scalars. 
  for (size_t i = 0; i < M; i++){
    ostr << "\tScalar #" << i << " = " << p.getScalar(i) << "\n";
  }

  // Print vectors. 
  for (size_t i = 0; i < N; i++){
    ostr << "\tVector #" << i << " = " << p.getVector(i) << "\n";
  }  
  
  return ostr;
}

#include <CD_NamespaceFooter.H>

#endif
