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
#include <array>

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
    @details Template parameter is the position in the m_scalars array. This is templated so that compilers may throw
    compile-time errors if trying to fetch elements out of range. 
    @return m_scalars[K]
  */
  template <size_t K>
  inline Real&
  real();

  /*!
    @brief Get one of the scalars. 
    @details Template parameter is the position in the m_scalars array. This is templated so that compilers may throw
    compile-time errors if trying to fetch elements out of range. 
    @return m_scalars[K]
  */
  template <size_t K>
  inline const Real&
  real() const;

  /*!
    @brief Get one of the RealVects.
    @details Template parameter is the position in the m_vectors array. This is templated so that compilers may throw
    compile-time errors if trying to fetch elements out of range. 
    @return m_vectors[K]
  */
  template <size_t K>
  inline RealVect&
  vect();

  /*!
    @brief Get one of the RealVects.
    @details Template parameter is the position in the m_vectors array. This is templated so that compilers may throw
    compile-time errors if trying to fetch elements out of range. 
    @return m_vectors[K]
  */
  template <size_t K>
  inline const RealVect&
  vect() const;

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
    @brief Particle comparison operator. Returns lexicographical ordering.
  */
  inline bool
  operator<(const GenericParticle<M, N>& a_other) const noexcept;

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
  std::array<Real, M> m_scalars;

  /*!
    @brief vector storage array
  */
  std::array<RealVect, N> m_vectors;
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
