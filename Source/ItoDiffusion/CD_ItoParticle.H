/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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
  @brief A particle class for use with ItoSolvers, i.e. drifting Brownian walkers.
  @details This class is used to encapsulate the requirements for running an ItoSolver. This computational particle contains position, mass, velocity, as well as a
  diffusion coefficient, a mobility, energy, and previous particle position. It also contains a temporary, and users can choose to allocate additional memory at 
  runtime. This memory is added on top of the other particle fields. In total, this class stores 5 + 3*SpaceDim fields. So in 3D this particle class contains
  14 Reals => 112 bytes. 
*/
class ItoParticle : public BinItem {
public:

  /*!
    @brief Set the number of run-time defined real variables. 
    @param[in] a_numRuntimeScalars Number of run-time scalar variables.
  */
  static void setNumRuntimeScalars(const int a_numRuntimeScalars);

  /*!
    @brief Set the number of run-time defined RealVect variables
    @param[in] a_numRuntimeVectors Number of run-time vector variables. 
  */
  static void setNumRuntimeVectors(const int a_numRuntimeVectors);
  
  /*!
    @brief Default constructor -- user should subsequently set the variables or call define.
    @note If the user has called setNumRuntimeScalars/Vectors then  will allocate runtime buffers
  */
  ItoParticle();

  /*!
    @brief Constructor. This calls the define function.
    @param[in] a_mass      Particle mass (i.e., computational weight)
    @param[in] a_position  Particle position
    @param[in] a_velocity  Particle velocity
    @param[in] a_diffusion Particle diffusion coefficient
    @param[in] a_mobility  Particle mobility coefficient
    @param[in] a_energy    Particle average energy
    @note If the user has called setNumRuntimeScalars/Vectors then  will allocate runtime buffers
  */
  ItoParticle(const Real      a_mass,
	      const RealVect& a_position,
	      const RealVect& a_velocity  = RealVect::Zero,
	      const Real      a_diffusion = 0.0,
	      const Real      a_mobility  = 1.0,
	      const Real      a_energy    = 0.0);

  /*!
    @brief Copy constructor. Copies all fields.
    @param[in] a_other Other particle. 
  */
  ItoParticle(const ItoParticle& a_other);

  /*!
    @brief Destructor -- this will deallocate memory so make sure you don't call setRuntimeScalars/setRuntimeVectors between particle creation and destruction. 
  */
  virtual ~ItoParticle();

  /*!
    @brief Full define function. 
    @param[in] a_mass      Particle mass (i.e., computational weight)
    @param[in] a_position  Particle position
    @param[in] a_velocity  Particle velocity
    @param[in] a_diffusion Particle diffusion coefficient
    @param[in] a_mobility  Particle mobility coefficient
    @param[in] a_energy    Particle average energy
  */
  inline
  void define(const Real      a_mass,
	      const RealVect& a_position,
	      const RealVect& a_velocity  = RealVect::Zero,
	      const Real      a_diffusion = 0.0,
	      const Real      a_mobility  = 1.0,
	      const Real      a_energy    = 0.0);

  /*!
    @brief Set the old position
    @details Set m_oldPosition.
    @param[in] a_oldPosition. 
  */
  inline
  void setOldPosition(const RealVect a_oldPosition);

  /*!
    @brief Get the old particle position
    @return m_oldPosition.
  */
  inline
  RealVect& oldPosition();

  /*!
    @brief Get the old particle position
    @return m_oldPosition.
  */
  inline
  const RealVect& oldPosition() const;

  /*!
    @brief Set particle mass
    @param[in] a_mass Particle mass
  */
  inline
  void setMass(const Real& a_mass);

  /*!
    @brief Get particle mass (i.e., computational weight)
    @return m_mass
  */
  inline
  Real& mass();

  /*! 
    @brief Get particle mass (i.e., computational weight)
    @return m_mass
  */
  inline
  const Real& mass() const;

  /*!
    @brief Set particle diffusion coefficient
  */
  inline
  void setDiffusion(const Real a_diffusion);

  /*!
    @brief Get diffusion
  */
  inline
  Real& diffusion();

  /*! 
    @brief Const retrieve diffusion
  */
  inline
  const Real& diffusion() const;

  /*!
    @brief Set particle mobility
  */
  inline
  void setMobility(const Real a_mobility);

  /*!
    @brief Get mobility
  */
  inline
  Real& mobility();

  /*! 
    @brief Const retrieve mobility
  */
  inline
  const Real& mobility() const;

  /*!
    @brief Get the particle conductivity
    @return Returns mobility*mass
  */
  inline
  Real conductivity() const;

  /*!
    @brief Get the particle diffusivity
    @return Returns diffusion*mass
  */
  inline
  Real diffusivity() const;  

  /*!
    @brief Set particle energy
  */
  inline
  void setEnergy(const Real a_energy);

  /*!
    @brief Get energy
  */
  inline
  Real& energy();

  /*! 
    @brief Const retrieve energy
  */
  inline
  const Real& energy() const;

  /*!
    @brief Set the velocity
  */
  inline
  void setVelocity(const RealVect& a_velocity);

  /*!
    @brief Set the velocity
  */
  inline
  void setVelocity(const Real& a_velocity, const int a_dir);

  /*!
    @brief Get reference velocity field
  */
  inline
  RealVect& velocity();

  /*!
    @brief Get reference velocity field. Const version. 
  */
  inline
  const RealVect& velocity() const;

  /*!
    @brief Get reference velocity field. Const version for one component only. 
  */
  inline
  Real velocity(const int a_dir) const;

  /*!
    @brief Get temporary
  */
  inline
  Real& tmp();

  /*! 
    @brief Const retrieve temporary
  */
  inline
  const Real& tmp() const;

  /*!
    @brief Get a runtime-defined variable
  */
  inline
  Real& runtimeScalar(const int a_num);

  /*!
    @brief Get a runtime-defined variable
  */
  inline
  const Real& runtimeScalar(const int a_num) const;

  /*!
    @brief Get a runtime-defined vector
  */
  inline
  RealVect& runtimeVector(const int a_num);

  /*!
    @brief Get a runtime-defined vector
  */
  inline
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
  static int s_numRuntimeScalars;

  /*!
    @brief Static for defining the number of runtime RealVect components
  */
  static int s_numRuntimeVectors;
  
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

#include <CD_ItoParticleImplem.H>

#endif

