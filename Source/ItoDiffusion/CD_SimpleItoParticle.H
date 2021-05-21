/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_SimpleItoParticleParticle.H
  @brief  Declaration of a particle which only holds position and mass. 
  @author Robert Marskar
*/

#ifndef CD_SimpleItoParticleParticle_H
#define  CD_SimpleItoParticleParticle_H

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
class SimpleItoParticle : public BinItem {
public:

  /*!
    @brief Default constructor
  */
  SimpleItoParticle();

  /*!
    @brief Constructor. Subsequently calls define
  */
  SimpleItoParticle(const Real a_mass, const RealVect a_position, const Real a_energy);

  /*!
    @brief Default dtor
  */
  ~SimpleItoParticle();

  /*!
    @brief For non-default ctor
  */
  void define(const Real a_mass, const RealVect a_position, const Real a_energy);

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
    @brief Set particle energy
  */
  void setEnergy(const Real a_energy);

  /*!
    @brief Get energy
  */
  Real& energy();

  /*! 
    @brief Const retrieve energy
  */
  const Real& energy() const;

  /*!
    @brief Comparison operator
  */
  bool operator==(const SimpleItoParticle& a_p) const;

  /*!
    @brief Comparison operator
  */
  bool operator==(const SimpleItoParticle* a_p) const;

  /*!
    @brief Comparison operator
  */
  bool operator!=(const SimpleItoParticle& a_p) const;

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
  friend std::ostream & operator<<(std::ostream& ostr,const SimpleItoParticle& p);

protected:

  Real m_mass;
  Real m_energy;
};

std::ostream & operator<<(std::ostream& ostr, const SimpleItoParticle& p);

#include <CD_NamespaceFooter.H>

#endif
