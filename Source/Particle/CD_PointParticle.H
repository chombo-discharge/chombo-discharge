/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_PointParticle.H
  @brief  Declaration of a computational point particle 
  @author Robert Marskar
*/

#ifndef CD_PointParticle_H
#define CD_PointParticle_H

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
  @brief A particle class for use with Ito diffusion, i.e. drift-diffusion brownian walkers
*/
class PointParticle : public BinItem {
public:

  /*!
    @brief Default constructor
  */
  PointParticle();

  /*!
    @brief Constructor. Subsequently calls define
  */
  PointParticle(const RealVect a_position, const Real a_mass);

  /*!
    @brief Default dtor
  */
  virtual ~PointParticle();

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

#include <CD_NamespaceFooter.H>

#endif
