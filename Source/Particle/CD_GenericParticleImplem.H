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

  for (auto& s : m_compileTimeScalars) {
    s = 0.0;
  }

  for (auto& v : m_compileTimeVectors) {
    v = RealVect::Zero;
  }

  this->setRuntimeStorage(0, 0);
}

template <size_t M, size_t N>
inline GenericParticle<M, N>::GenericParticle(const GenericParticle<M, N>& a_other)
{
  this->m_position = a_other.m_position;
  this->m_compileTimeScalars  = a_other.m_compileTimeScalars;
  this->m_compileTimeVectors  = a_other.m_compileTimeVectors;

  this->m_runTimeScalars = a_other.m_runTimeScalars;
  this->m_runTimeVectors = a_other.m_runTimeVectors;
}

template <size_t M, size_t N>
inline GenericParticle<M, N>::~GenericParticle()
{}

template <size_t M, size_t N>
inline void
GenericParticle<M, N>::setRuntimeStorage(const unsigned int a_numReal, const unsigned int a_numRealVect) noexcept
{
  m_runTimeScalars.resize(a_numReal);
  m_runTimeVectors.resize(a_numRealVect);
}

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
  return std::get<K>(m_compileTimeScalars);
}

template <size_t M, size_t N>
template <size_t K>
inline const Real&
GenericParticle<M, N>::real() const
{
  return std::get<K>(m_compileTimeScalars);
}

template <size_t M, size_t N>
template <size_t K>
inline RealVect&
GenericParticle<M, N>::vect()
{
  return std::get<K>(m_compileTimeVectors);
}

template <size_t M, size_t N>
template <size_t K>
inline const RealVect&
GenericParticle<M, N>::vect() const
{
  return std::get<K>(m_compileTimeVectors);
}

template <size_t M, size_t N>
inline bool
GenericParticle<M, N>::operator==(const GenericParticle<M, N>& a_p) const
{
  return (this->m_position == a_p.m_position && this->m_compileTimeScalars == a_p.m_compileTimeScalars && this->m_compileTimeVectors == a_p.m_compileTimeVectors);
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
  // Compile size is the position plus M Reals and N RealVects
  const int compileSize =  SpaceDim * sizeof(Real) + M * sizeof(Real) + N * SpaceDim * sizeof(Real);

  // Runtime size if number of runtime reals and realvects:
  const int runtimeSize = m_runTimeScalars.size() * sizeof(Real) + m_runTimeVectors.size() * SpaceDim * sizeof(Real);

  return compileSize + runtimeSize;
}

template <size_t M, size_t N>
inline void
GenericParticle<M, N>::linearOut(void* buf) const
{
  Real* buffer = (Real*)buf;

  // clang-format off

  // Linearize m_position onto buffer.
  D_TERM(*buffer++ = m_position[0];,
	 *buffer++ = m_position[1];,
	 *buffer++ = m_position[2];);

  // Linearize m_compileTimeScalars onto the buffer
  for (size_t i = 0; i < M; i++) {
    *buffer++ = m_compileTimeScalars[i];
  }

  // Linearize vectors onto the buffer
  for (size_t i = 0; i < N; i++) {
    const RealVect& v = m_compileTimeVectors[i];

    D_TERM(*buffer++ = v[0];,
	   *buffer++ = v[1];,
	   *buffer++ = v[2];);
  }

  // Linearize the run-time scalars onto the buffer.
  for (size_t i = 0; i < m_runTimeScalars.size(); i++) {
    *buffer++ = m_runTimeScalars[i];
  }

  // Linearize the run-time vectors onto the buffer.
  for (size_t i = 0; i < m_runTimeVectors.size(); i++) {
    const RealVect& v = m_runTimeVectors[i];
    
    D_TERM(*buffer++ = v[0];,
	   *buffer++ = v[1];,
	   *buffer++ = v[2];);
  }
  // clang-format on
}

template <size_t M, size_t N>
inline void
GenericParticle<M, N>::linearIn(void* buf)
{
  Real* buffer = (Real*)buf;

  // clang-format off  
  D_TERM(m_position[0] = *buffer++;,
	 m_position[1] = *buffer++;,
	 m_position[2] = *buffer++);

  // Linearize buffer onto scalars
  for (size_t i = 0; i < M; i++) {
    m_compileTimeScalars[i] = *buffer++;
  }

  // Linearize buffer onto vectors
  for (size_t i = 0; i < N; i++) {
    RealVect& v = m_compileTimeVectors[i];

    D_TERM(v[0] = *buffer++;,
	   v[1] = *buffer++;,
	   v[2] = *buffer++;);
  }

  // Linearize the run-time scalars onto the buffer.
  for (size_t i = 0; i < m_runTimeScalars.size(); i++) {
    m_runTimeScalars[i] = *buffer++;
  }

  // Linearize the run-time vectors onto the buffer.
  for (size_t i = 0; i < m_runTimeVectors.size(); i++) {
    RealVect& v = m_runTimeVectors[i];

    D_TERM(v[0] = *buffer++;,
	   v[1] = *buffer++;,
	   v[2] = *buffer++;);    
  }
  // clang-format on
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
