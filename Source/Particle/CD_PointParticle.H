/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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
  @brief A particle class that only has a position and a mass.
*/
class PointParticle : public BinItem {
public:

  /*!
    @brief Default constructor. Sets mass and position to zero.
  */
  PointParticle();

  /*!
    @brief Constructor. Subsequently calls define
    @param[in] a_position Position (in physical space)
    @param[in] a_mass     Particle mass
  */
  PointParticle(const RealVect a_position, const Real a_mass);

  /*!
    @brief Dtor - does nothing.
  */
  virtual ~PointParticle();

  /*!
    @brief Full define function
    @param[in] a_position Position (in physical space)
    @param[in] a_mass     Particle mass
  */
  void define(const RealVect a_position, const Real a_mass);

  /*!
    @brief Set particle mass
    @param[in] a_mass Particle mass
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
    sufficient memory for the buffer has already been allocated by the caller
    @param[out] buf Particle 
  */
  virtual void linearOut(void* buf) const;

  /*!
    @brief Read a linear binary representation of the internal data. Assumes that the buffer has the correct data.
    @param[in] buf Other particle.
  */
  virtual void linearIn(void* buf);

protected:

  /*!
    @brief Computational particle mass
  */
  Real m_mass;
};

#include <CD_NamespaceFooter.H>

#endif
