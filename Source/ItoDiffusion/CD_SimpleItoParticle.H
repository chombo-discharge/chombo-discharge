/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SimpleItoParticleParticle.H
  @brief  Declaration of a particle which only holds position, mass, and energy. 
  @author Robert Marskar
*/

#ifndef CD_SimpleItoParticleParticle_H
#define CD_SimpleItoParticleParticle_H

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
  @brief A particle class for reducing HDF5 I/O file sizes when using ItoSolver. 
  @details This is a computational particle class which only stores the mass, position and energy of a particle. It is used
  with ItoSolver in order to reduce the file size of HDF5 files. 
*/
class SimpleItoParticle : public BinItem
{
public:
  /*!
    @brief Default constructor. Should subsequently set mass/position/energy.
  */
  SimpleItoParticle();

  /*!
    @brief Constructor. Sets mass, position, and energy.
    @param[in] a_mass     Computational mass/weight
    @param[in] a_position Position in physical space. 
    @param[in] a_energy   Average particle energy
  */
  SimpleItoParticle(const Real a_mass, const RealVect a_position, const Real a_energy);

  /*!
    @brief Destructor (does nothing)
  */
  ~SimpleItoParticle();

  /*!
    @brief Define function. Sets mass, position, and energy.
    @param[in] a_mass     Computational mass/weight
    @param[in] a_position Position in physical space. 
    @param[in] a_energy   Average particle energy
  */
  void
  define(const Real a_mass, const RealVect a_position, const Real a_energy);

  /*!
    @brief Set particle mass
    @param[in] a_mass Computational weight
    @details Sets m_mass
  */
  void
  setMass(const Real a_mass);

  /*!
    @brief Get mass
    @return Returns m_mass
  */
  Real&
  mass();

  /*! 
    @brief Const retrieve mass
    @return Returns m_mass
  */
  const Real&
  mass() const;

  /*!
    @brief Set particle energy
    @param[in] a_energy Average particle energy
    @details Sets m_energy
  */
  void
  setEnergy(const Real a_energy);

  /*!
    @brief Get energy
    @return Returns m_energy
  */
  Real&
  energy();

  /*! 
    @brief Const retrieve energy
    @return Returns m_energy
  */
  const Real&
  energy() const;

  /*!
    @brief Comparison operator. Returns true if mass, energy, and position are the same. 
    @param[in] a_p Other particle
  */
  bool
  operator==(const SimpleItoParticle& a_p) const;

  /*!
    @brief Comparison operator. Returns true if mass, energy, and position are the same. 
    @param[in] a_p Other particle
  */
  bool
  operator==(const SimpleItoParticle* a_p) const;

  /*!
    @brief Comparison operator. Returns false if mass, energy, or position are not the same.
    @param[in] a_p Other particle
  */
  bool
  operator!=(const SimpleItoParticle& a_p) const;

  /*!
    @brief Returns the size, in number of bytes, of a flat representation of the data in this object.
  */
  virtual int
  size() const;

  /*!
    @brief Write a linear binary representation of the internal data. Assumes that
    sufficient memory for the buffer has already been allocated by
    the caller.
  */
  virtual void
  linearOut(void* buf) const;

  /*!
    @brief Read a linear binary representation of the internal data.
    Assumes that the buffer has the correct data.
  */
  virtual void
  linearIn(void* buf);

  /*!
    @brief For outputting particle.
  */
  friend std::ostream&
  operator<<(std::ostream& ostr, const SimpleItoParticle& p);

protected:
  /*!
    @brief Particle mass (i.e., computational weight)
  */
  Real m_mass;

  /*!
    @brief Particle average energy
  */
  Real m_energy;
};

/*!
  @brief For printing particle to screen
*/
std::ostream&
operator<<(std::ostream& ostr, const SimpleItoParticle& p);

#include <CD_NamespaceFooter.H>

#endif
