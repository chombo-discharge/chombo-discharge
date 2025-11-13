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
#include <cstdint>

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

namespace detail {
  /*!
    @brief Helper function for linearizing a particle onto a buffer.
    @detail This copies T into the buffer p and increments the pointer.
    @param[in,out] p Buffer
    @param[in] a_x Particle property. Must be a trivial type.
  */
  template <class T>
  inline void
  pushParticleProperty(uint8_t*& p, const T& a_x)
  {
    static_assert(std::is_trivially_copyable<T>::value,
                  "detail::pushParticleProperty expects trivially-copyable fields");

    std::memcpy(p, &a_x, sizeof(T));

    p += sizeof(T);
  };

  /*!
    @brief Helper function for linearizing a buffer onto a particle property.
    @detail This copies from the buffer p into the particle property T. 
    @param[in,out] p Buffer. 
    @param[in] a_x Particle property. Must be a trivial type.
  */
  template <class T>
  inline void
  pullParticleProperty(const uint8_t*& p, T& a_x)
  {
    static_assert(std::is_trivially_copyable<T>::value,
                  "detail::pushParticleProperty expects trivially-copyable fields");

    std::memcpy(&a_x, p, sizeof(T));

    p += sizeof(T);
  };
} // namespace detail

/*!
  @brief A generic particle class, holding the position and a specified number of real and vector values.
  @details For the sake of simplicity, I have added separate linearization functions H5size, H5linearOut, H5linearIn that only linearizes
  the Real components of the particle. This is done for convenience only, so I don't have to rewrite all of Chombo's particle HDF5 checkpointing
  routines. Instead, we simply run a bootleg version of those codes that call the H5 linearization versions instead of the regular linearization
  functions. Users can use this functionality for trimming other types of unnecessary properties that don't need to go into the HDF5 files.
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
    @brief Get the particle ID
  */
  inline const int32_t&
  particleID() const noexcept;

  /*!
    @brief Get the particle ID
  */
  inline int32_t&
  particleID() noexcept;

  /*!
    @brief Get the MPI rank ID
  */
  inline const int32_t&
  rankID() const noexcept;

  /*!
    @brief Get the MPI rank ID
  */
  inline int32_t&
  rankID() noexcept;

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
    @brief Get the M scalars
    @return m_scalars
  */
  inline const std::array<Real, M>&
  getReals() const noexcept;

  /*!
    @brief Get the M scalars
    @return m_scalars
  */
  inline std::array<Real, M>&
  getReals() noexcept;

  /*!
    @brief Get the N vectors
    @return m_vects
  */
  inline const std::array<RealVect, N>&
  getVects() const noexcept;

  /*!
    @brief Get the N vectors
    @return m_vects
  */
  inline std::array<RealVect, N>&
  getVects() noexcept;

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
  linearOut(void* const a_buffer) const;

  /*!
    @brief Read a linear binary representation of the internal data. Assumes that the buffer has the correct data.
    @param[in] a_buffer Pointer to memory block.
  */
  inline virtual void
  linearIn(const void* const a_buffer);

  /*!
    @brief Function that is used when writing particles to HDF5.
    @details This ignores the particleID and rankID, so we can directory plug into Chombo's HDF5 APIs.
    @note Please do not try to output anything other than Real components. Chombo will expect only Reals. 
  */
  inline virtual int
  H5size() const;

  /*!
    @brief Linearize the Real components onto a buffer workable by HDF5.
    @details Ignores the particleID and rankID, so we can directly plug into Chombo's HDF5 API.
    @param[in] a_buffer Pointer to memory block.
    @note Please do not try to output anything other than Real components. Chombo will expect only Reals.     
  */
  inline virtual void
  H5linearOut(void* const a_buffer) const;

  /*!
    @brief Delinearize the buffer onto the real components in the particle class.
    @details Ignores the particleID and rankID, so we can directly plug into Chombo's HDF5 API.
    @param[in] a_buffer Pointer to memory block.
    @note Please do not try to output anything other than Real components. Chombo will expect only Reals.     
  */
  inline virtual void
  H5linearIn(const void* const a_buffer);

protected:
  /*!
    @brief Particle ID
  */
  int32_t m_particleID;

  /*!
    @brief MPI rank owning this particle
  */
  int32_t m_rankID;

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
