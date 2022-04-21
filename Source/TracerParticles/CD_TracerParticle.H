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
  @brief A tracer particle class. 
*/
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
    @return m_state
  */
  inline
  Real& state();      

  /*!
    @brief Check if particle should be processed or not
    @return m_state
  */
  inline
  const Real& state() const;

  /*!
    @brief Get the particle "mass"
    @return m_mass
  */
  inline
  Real& mass();

  /*!
    @brief Get the particle "mass"
    @return m_mass
  */
  inline
  const Real& mass() const;

  /*!
    @brief Get the particle velocity. 
    @return m_velocity
  */
  inline
  RealVect& velocity();

  /*!
    @brief Get the particle velocity
    @return m_velocity
  */
  inline
  const RealVect& velocity() const;

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

protected:

  /*!
    @brief Particle "mass".
  */
  Real m_mass;

  /*!
    @brief Particle velocity
  */
  RealVect m_velocity;

  /*!
    @brief Process particle or not
  */
  Real m_state;
};

/*!
  @brief Particle printing function
  @param[in] ostr       Output stream
  @param[in] a_particle Particle to output
*/
inline
std::ostream & operator<<(std::ostream& ostr, const TracerParticle& p);

#include <CD_NamespaceFooter.H>

#include <CD_TracerParticleImplem.H>

#endif

