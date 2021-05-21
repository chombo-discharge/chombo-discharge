/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ItoParticle.H
  @brief  Declaration of a particle class for Ito diffusion
  @author Robert Marskar
*/

#ifndef CD_ItoParticle_H
#define CD_ItoParticle_H

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
class ItoParticle : public BinItem {
public:

  /*!
    @brief Set the number of run-time defined real variables
  */
  static void setNumRunTimeScalars(const int a_num);

  /*!
    @brief Set the number of run-time defined RealVect variables
  */
  static void setNumRunTimeVectors(const int a_num);
  
  /*!
    @brief Default constructor
  */
  ItoParticle();

  /*!
    @brief Constructor. Subsequently calls define
  */
  ItoParticle(const Real      a_mass,
	       const RealVect& a_position,
	       const RealVect& a_velocity  = RealVect::Zero,
	       const Real      a_diffusion = 0.0,
	       const Real      a_mobility  = 1.0,
	       const Real      a_energy    = 0.0);

  /*!
    @brief Copy constructor
  */
  ItoParticle(const ItoParticle& a_other);

  /*!
    @brief Default dtor
  */
  virtual ~ItoParticle();

  /*!
    @brief For non-default ctor
  */
  void define(const Real      a_mass,
	      const RealVect& a_position,
	      const RealVect& a_velocity  = RealVect::Zero,
	      const Real      a_diffusion = 0.0,
	      const Real      a_mobility  = 1.0,
	      const Real      a_energy    = 0.0);



  /*!
    @brief Set the old position
  */
  void setOldPosition(const RealVect a_oldPosition);

  /*!
    @brief Get the old position
  */
  RealVect& oldPosition();

  /*!
    @brief Old position
  */
  const RealVect& oldPosition() const;

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
    @brief Set particle diffusion
  */
  void setDiffusion(const Real a_diffusion);

  /*!
    @brief Get diffusion
  */
  Real& diffusion();

  /*! 
    @brief Const retrieve diffusion
  */
  const Real& diffusion() const;

  /*!
    @brief Set particle mobility
  */
  void setMobility(const Real a_mobility);

  /*!
    @brief Get mobility
  */
  Real& mobility();

  /*! 
    @brief Const retrieve mobility
  */
  const Real& mobility() const;

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
    @brief Set the velocity
  */
  void setVelocity(const RealVect& a_velocity);

  /*!
    @brief Set the velocity
  */
  void setVelocity(const Real& a_velocity, const int a_dir);

  /*!
    @brief Get reference velocity field
  */
  RealVect& velocity();

  /*!
    @brief Get reference velocity field. Const version. 
  */
  const RealVect& velocity() const;

  /*!
    @brief Get reference velocity field. Const version for one component only. 
  */
  Real velocity(const int a_dir) const;

  /*!
    @brief Get temporary
  */
  Real& tmp();

  /*! 
    @brief Const retrieve temporary
  */
  const Real& tmp() const;

  /*!
    @brief Get a runtime-defined variable
  */
  Real& runtimeScalar(const int a_num);

  /*!
    @brief Get a runtime-defined variable
  */
  const Real& runtimeScalar(const int a_num) const;

  /*!
    @brief Get a runtime-defined vector
  */
  RealVect& runtimeVector(const int a_num);

  /*!
    @brief Get a runtime-defined vector
  */
  const RealVect& runtimeVector(const int a_num) const;

  /*!
    @brief Comparison operator
  */
  bool operator==(const ItoParticle& a_p) const;

  /*!
    @brief Comparison operator
  */
  bool operator==(const ItoParticle* a_p) const;

  /*!
    @brief Comparison operator
  */
  bool operator!=(const ItoParticle& a_p) const;

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
  friend std::ostream & operator<<(std::ostream& ostr,const ItoParticle& p);

protected:

  /*!
    @brief Static for defining the number of runtime Real components
  */
  static int s_num_runtimeScalars;

  /*!
    @brief Static for defining the number of runtime RealVect components
  */
  static int s_num_runtimeVectors;
  
  Real     m_tmp; // Not communicated!
  Real     m_mass;
  Real     m_diffusion;
  Real     m_energy;
  Real     m_mobility;
  RealVect m_velocity;
  RealVect m_oldPosition;

  Real*     m_runtimeScalars;
  RealVect* m_runtimeVectors;

  /*!
    @brief Set the run-time components
  */
  void allocateRuntimeBuffers();
};

std::ostream & operator<<(std::ostream& ostr, const ItoParticle& p);

#include <CD_NamespaceFooter.H>

#endif

