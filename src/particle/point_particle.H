/*!
  @file   point_particle.H
  @brief  Declaration of a particle class for Ito diffusion
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _POINT_PARTICLE_
#define _POINT_PARTICLE_

#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

#include "BinItem.H"


/*!
  @brief A particle class for use with Ito diffusion, i.e. drift-diffusion brownian walkers
*/
class point_particle : public BinItem {
public:

  /*!
    @brief Default constructor
  */
  point_particle();

  /*!
    @brief Constructor. Subsequently calls define
  */
  point_particle(const RealVect a_position, const Real a_mass);

  /*!
    @brief Default dtor
  */
  virtual ~point_particle();

  /*!
    @brief For non-default ctor
  */
  void define(const RealVect a_position, const Real a_mass);

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

protected:
  
  Real m_mass;
};

#endif

