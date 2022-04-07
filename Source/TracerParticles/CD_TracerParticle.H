/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TracerParticle.H
  @brief  Declaration of a tracer particle class which can hold N components.
  @author Robert Marskar
*/

#ifndef CD_TracerParticle_H
#define CD_TracerParticle_H

// Std includes
#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

// Chombo includes
#include <BinItem.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief A tracer particle class. The particle class can contain M independent scalar and N independent vector components.
*/
template <size_t M, size_t N>
class TracerParticle : public BinItem {
public:

  /*!
    @brief Default constructor -- initializes everything to zero. 
  */
  inline  
  TracerParticle();

  /*!
    @brief Copy constructor. Copies all fields.
    @param[in] a_other Other particle. 
  */
  inline  
  TracerParticle(const TracerParticle& a_other);

  /*!
    @brief Destructor (deallocates runtime memory storage)
  */
  inline  
  virtual ~TracerParticle();

  /*!
    @brief Check if particle should be processed or not
    @return m_process
  */
  inline
  bool& processParticle();

  /*!
    @brief Check if particle should be processed or not
    @return m_process
  */
  inline
  const bool& processParticle() const;

  /*!
    @brief Get a particle scalar. 
    @details Template because we want to have a function getScalar() available but the other option 'getScalar(int)' is 
    not compile-time constant.
  */
  template <size_t which>
  inline
  Real& getScalar();

  /*!
    @brief Get a particle scalar. 
    @details Template because we want to have a function getScalar() available but the other option 'getScalar(int)' is 
    not compile-time constant.
  */
  template <size_t which>
  inline
  const Real& getScalar() const;
  
  /*!
    @brief Get a particle scalar. 
    @param[in] a_which Which particle scalar to fetch. 
    @return m_scalars[a_which]
  */
  inline
  Real& getScalar(const size_t a_which);

  /*!
    @brief Get a particle scalar. 
    @param[in] a_which Which particle scalar to fetch. 
    @return m_scalars[a_which]
  */
  inline
  const Real& getScalar(const size_t a_which) const;  

  /*!
    @brief Get a particle vector
    @param[in] a_which Which particle vector to fetch. 
    @return m_vectors[a_which]
  */
  inline
  RealVect& getVector(const size_t a_which);

  /*!
    @brief Get a particle vector
    @param[in] a_which Which particle vector to fetch. 
    @return m_vectors[a_which]
  */
  inline
  const RealVect& getVector(const size_t a_which) const;

  /*!
    @brief Comparison operator with other particle
    @param[in] a_other Other particle
  */
  inline  
  bool operator==(const TracerParticle& a_other) const;

  /*!
    @brief Comparison operator with other particle
    @param[in] a_other Other particle
  */
  inline  
  bool operator!=(const TracerParticle& a_other) const;

  /*!
    @brief Returns the size, in number of bytes, of a flat representation of the data in this object.
  */
  inline
  virtual int size() const override;

  /*!
    @brief Write a linear binary representation of the internal data. Assumes that
    sufficient memory for the buffer has already been allocated by
    the caller.
    @param[in] a_buf Pointer to memory block
  */
  inline
  virtual void linearOut(void* buf) const override;

  /*!
    @brief Read a linear binary representation of the internal data.
    Assumes that the buffer has the correct data.
    @param[in] a_buf Pointer to memory block
  */
  inline  
  virtual void linearIn(void* a_buf);

  /*!
    @brief Friend function for letting ostream print the particle.
    @param[in] ostr       Output stream
    @param[in] a_particle Particle to output
  */
  //  inline  
  //  friend std::ostream& operator<<(std::ostream& ostr, const TracerParticle<M, N>& a_particle);
  template<size_t i>  
  inline
  void f(){

  }

protected:

  /*!
    @brief Scalar components
  */
  std::array<Real, M> m_scalars;

  /*!
    @brief Vector components
  */
  std::array<RealVect, N> m_vectors;

  /*!
    @brief Process particle or not
  */
  bool m_process;
};

/*!
  @brief Particle printing function
  @param[in] ostr       Output stream
  @param[in] a_particle Particle to output
*/
template <size_t M, size_t N>
inline
std::ostream & operator<<(std::ostream& ostr, const TracerParticle<M, N>& p);

#include <CD_NamespaceFooter.H>

#include <CD_TracerParticleImplem.H>

#endif

