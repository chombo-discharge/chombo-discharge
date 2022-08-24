/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NonCommParticle.H
  @brief  Declaration of a generic particle class. 
  @author Robert Marskar
*/

#ifndef CD_NonCommParticle_H
#define CD_NonCommParticle_H

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
  @brief A generic non-communicated particle class. 
*/
template <size_t M, size_t N>
class NonCommParticle
{
public:
  /*!
    @brief Default constructor -- initializes everything to zero. 
  */
  inline NonCommParticle();

  /*!
    @brief Copy constructor. Copies all fields.
    @param[in] a_other Other particle. 
  */
  inline NonCommParticle(const NonCommParticle<M, N>& a_other);

  /*!
    @brief Destructor (deallocates runtime memory storage)
  */
  inline virtual ~NonCommParticle();

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
  operator==(const NonCommParticle<M, N>& a_other) const;

  /*!
    @brief Comparison operator with other particle.
    @param[in] a_other Other particle
    @note Does not compare the M and N Reals and RealVects. 
  */
  inline bool
  operator!=(const NonCommParticle<M, N>& a_other) const;

protected:

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
operator<<(std::ostream& ostr, const NonCommParticle<M, N>& p);

#include <CD_NamespaceFooter.H>

#include <CD_NonCommParticleImplem.H>

#endif
