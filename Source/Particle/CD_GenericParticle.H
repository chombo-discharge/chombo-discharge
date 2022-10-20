/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_GenericParticle.H
  @brief  Declaration of a generic particle class. 
  @author Robert Marskar
*/

#ifndef CD_GenericParticle_H
#define CD_GenericParticle_H

// Std includes
#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief A generic particle class, holding the position and a specified number of real and vector values.
*/
template <size_t M, size_t N>
class GenericParticle
{
public:
  /*!
    @brief Default constructor -- initializes everything to zero. 
  */
  inline GenericParticle();

  /*!
    @brief Copy constructor. Copies all fields.
    @param[in] a_other Other particle. 
  */
  inline GenericParticle(const GenericParticle<M, N>& a_other);

  /*!
    @brief Destructor (deallocates runtime memory storage)
  */
  inline virtual ~GenericParticle();

  /*!
    @brief Set up particle run-time storage
    @param[in] a_numReal     Number of extra Reals to allocate
    @param[in] a_numRealVect Number of extra RealVects to allocate
    @note Default is zero, but if set to anything other than that, the user MUST ensure that all particles have the same number of run-time components. 
  */
  inline void
  setRuntimeStorage(const unsigned int a_numReal, const unsigned int a_numRealVect) noexcept;

  /*!
    @brief Get the particle position. 
    @return m_position
  */
  inline RealVect&
  position();

  /*!
    @brief Get the particle position.
    @return m_position
  */
  inline const RealVect&
  position() const;

  /*!
    @brief Get one of the scalars. 
    @details Template parameter is the position in the m_compileTimeScalars array. This is templated so that compilers may throw
    compile-time errors if trying to fetch elements out of range. 
    @return m_compileTimeScalars[K]
  */
  template <size_t K>
  inline Real&
  real();

  /*!
    @brief Get one of the scalars. 
    @details Template parameter is the position in the m_compileTimeScalars array. This is templated so that compilers may throw
    compile-time errors if trying to fetch elements out of range. 
    @return m_compileTimeScalars[K]
  */
  template <size_t K>
  inline const Real&
  real() const;

  /*!
    @brief Get one of the RealVects.
    @details Template parameter is the position in the m_compileTimeVectors array. This is templated so that compilers may throw
    compile-time errors if trying to fetch elements out of range. 
    @return m_compileTimeVectors[K]
  */
  template <size_t K>
  inline RealVect&
  vect();

  /*!
    @brief Get one of the RealVects.
    @details Template parameter is the position in the m_compileTimeVectors array. This is templated so that compilers may throw
    compile-time errors if trying to fetch elements out of range. 
    @return m_compileTimeVectors[K]
  */
  template <size_t K>
  inline const RealVect&
  vect() const;

  /*!
    @brief Get one of the run-time configurable reals
    @param[in] index Variable to fetch
  */
  inline Real&
  runtimeReal(const size_t index) noexcept;

  /*!
    @brief Get one of the run-time configurable reals
    @param[in] index Variable to fetch
  */
  inline const Real&
  runtimeReal(const size_t index) const noexcept;

  /*!
    @brief Get one of the run-time configurable RealVects
    @param[in] index Variable to fetch
  */
  inline RealVect&
  runtimeVect(const size_t index) noexcept;

  /*!
    @brief Get one of the run-time configurable RealVects
    @param[in] index Variable to fetch
  */
  inline const RealVect&
  runtimeVect(const size_t index) const noexcept;

  /*!
    @brief Comparison operator with other particle. 
    @param[in] a_other Other particle
    @note Does not compare the M and N Reals and RealVects. 
  */
  inline bool
  operator==(const GenericParticle<M, N>& a_other) const;

  /*!
    @brief Comparison operator with other particle.
    @param[in] a_other Other particle
    @note Does not compare the M and N Reals and RealVects. 
  */
  inline bool
  operator!=(const GenericParticle<M, N>& a_other) const;

  /*!
    @brief Returns the size, in number of bytes, of a flat representation of the data in this object.
  */
  inline virtual int
  size() const;

  /*!
    @brief Write a linear binary representation of the internal data. Assumes that
    sufficient memory for the buffer has already been allocated by the caller. 
    @param[in] a_buffer Pointer to memory block
  */
  inline virtual void
  linearOut(void* a_buffer) const;

  /*!
    @brief Read a linear binary representation of the internal data. Assumes that the buffer has the correct data.
    @param[in] a_buffer Pointer to memory block
  */
  inline virtual void
  linearIn(void* a_buffer);

protected:
  /*!
    @brief Particle position
  */
  RealVect m_position;

  /*!
    @brief Scalar storage array
  */
  std::array<Real, M> m_compileTimeScalars;

  /*!
    @brief vector storage array
  */
  std::array<RealVect, N> m_compileTimeVectors;

  /*!
    @brief Run-time configurable Real storage
  */
  std::vector<Real> m_runTimeScalars;

  /*!
    @brief Run-time configurable RealVect storage
  */
  std::vector<RealVect> m_runTimeVectors;
};

/*!
  @brief Particle printing function
  @param[in] ostr       Output stream
  @param[in] a_particle Particle to output
*/
template <size_t M, size_t N>
inline std::ostream&
operator<<(std::ostream& ostr, const GenericParticle<M, N>& p);

#include <CD_NamespaceFooter.H>

#include <CD_GenericParticleImplem.H>

#endif
