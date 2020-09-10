/*!
  @file   massive_particle.H
  @brief  Declaration of a particle which only holds position and mass. 
  @author Robert Marskar
  @date   Sept. 2020
*/

#ifndef _MASSIVE_PARTICLE_H_
#define _MASSIVE_PARTICLE_H_

#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

#include "BinItem.H"

/*!
  @brief A particle class for use with Ito diffusion, i.e. drift-diffusion brownian walkers
*/
class massive_particle : public BinItem {
public:

  /*!
    @brief Default constructor
  */
  massive_particle();

  /*!
    @brief Constructor. Subsequently calls define
  */
  massive_particle(const Real a_mass, const RealVect a_position);

  /*!
    @brief Default dtor
  */
  ~massive_particle();

  /*!
    @brief For non-default ctor
  */
  void define(const Real a_mass, const RealVect a_position);

  /*!
    @brief Set particle mass
  */
  void setMass(const Real a_mass);

  /*!
    @brief Get mass
  */
  Real& mass();

  /*! 
    @brief Const retrieve mass
  */
  const Real& mass() const;

  /*!
    @brief Comparison operator
  */
  bool operator==(const massive_particle& a_p) const;

  /*!
    @brief Comparison operator
  */
  bool operator==(const massive_particle* a_p) const;

  /*!
    @brief Comparison operator
  */
  bool operator!=(const massive_particle& a_p) const;

  /*!
    @brief Returns the size, in number of bytes, of a flat representation of the data in this object.
  */
  virtual int size() const;

  /*!
    @brief Write a linear binary representation of the internal data. Assumes that
    sufficient memory for the buffer has already been allocated by
    the caller.
  */
  virtual void linearOut(void* buf) const;

  /*!
    @brief Read a linear binary representation of the internal data.
    Assumes that the buffer has the correct data.
  */
  virtual void linearIn(void* buf);

  /*!
    @brief For outputting particle
  */
  friend std::ostream & operator<<(std::ostream& ostr,const massive_particle& p);

protected:

  Real     m_mass;
};

std::ostream & operator<<(std::ostream& ostr, const massive_particle& p);

#endif
