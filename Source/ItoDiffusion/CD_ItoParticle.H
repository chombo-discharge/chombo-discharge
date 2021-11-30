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
    @param[in] a_averageEnergy    Particle average energy
    @note If the user has called setNumRuntimeScalars/Vectors then  will allocate runtime buffers
  */
  ItoParticle(const Real      a_mass,
	      const RealVect& a_position,
	      const RealVect& a_velocity  = RealVect::Zero,
	      const Real      a_diffusion = 0.0,
	      const Real      a_mobility  = 1.0,
	      const Real      a_averageEnergy    = 0.0);

  /*!
    @brief Copy constructor. Copies all fields.
    @param[in] a_other Other particle. 
  */
  ItoParticle(const ItoParticle& a_other);

  /*!
    @brief Destructor (deallocates runtime memory storage)
  */
  virtual ~ItoParticle();

  /*!
    @brief Full define function. 
    @param[in] a_mass          Particle mass (i.e., computational weight)
    @param[in] a_position      Particle position
    @param[in] a_velocity      Particle velocity
    @param[in] a_diffusion     Particle diffusion coefficient
    @param[in] a_mobility      Particle mobility coefficient
    @param[in] a_averageEnergy Particle average energy
  */
  inline
  void define(const Real      a_mass,
	      const RealVect& a_position,
	      const RealVect& a_velocity      = RealVect::Zero,
	      const Real      a_diffusion     = 0.0,
	      const Real      a_mobility      = 1.0,
	      const Real      a_averageEnergy = 0.0);

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
    @param[in] a_diffusion Diffusion coefficient
  */
  inline
  void setDiffusion(const Real a_diffusion);

  /*!
    @brief Get particle diffusion coefficient
    @return m_diffusion
  */
  inline
  Real& diffusion();

  /*! 
    @brief Get particle diffusion coefficient
    @return m_diffusion
  */
  inline
  const Real& diffusion() const;

  /*!
    @brief Set particle mobility
    @param[in] a_mobility Mobility coefficient
  */
  inline
  void setMobility(const Real a_mobility);

  /*!
    @brief Get mobility coefficient
    @return m_mobility
  */
  inline
  Real& mobility();

  /*! 
    @brief Get mobility coefficient
    @return m_mobility
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
    @brief Set average particle energy
    @param[in] a_averageEnergy Average particle energy.
  */
  inline
  void setEnergy(const Real a_averageEnergy);

  /*!
    @brief Get average particle energy
    @return m_averageEnergy
  */
  inline
  Real& energy();

  /*! 
    @brief Get average particle energy
    @return m_averageEnergy
  */
  inline
  const Real& energy() const;

  /*!
    @brief Get the total energy
    @return Returns m_energy*m_mass
  */
  inline
  Real totalEnergy() const;

  /*!
    @brief Set the particle velocity
    @param[in] a_velocity Particle velocity
  */
  inline
  void setVelocity(const RealVect& a_velocity);

  /*!
    @brief Set the velocity in particular direction
    @param[in] a_velocity Velocity component
    @param[in] a_dir      Coordinate direction
  */
  inline
  void setVelocity(const Real& a_velocity, const int a_dir);

  /*!
    @brief Get the particle velocity
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
    @brief Get the particle velocity in a particular coordinate direction
    @param[in] a_dir Coordinate direction
    @return m_velocity[a_dir]
  */
  inline
  Real velocity(const int a_dir) const;

  /*!
    @brief Get temporary storage
    @return m_tmp
  */
  inline
  Real& tmp();

  /*! 
    @brief Get temporary storage
    @return m_tmp
  */
  inline
  const Real& tmp() const;

  /*!
    @brief Get one of the runtime-defined scalar variable
    @param[in] a_which Which run-time scalar to obtain
  */
  inline
  Real& runtimeScalar(const int a_which);

  /*!
    @brief Get one of the runtime-defined scalar variable
    @param[in] a_which Which run-time scalar to obtain
  */
  inline
  const Real& runtimeScalar(const int a_which) const;

  /*!
    @brief Get one of the runtime-defined vector variable
    @param[in] a_which Which run-time vector to obtain
  */
  inline
  RealVect& runtimeVector(const int a_which);

  /*!
    @brief Get one of the runtime-defined vector variable
    @param[in] a_which Which run-time vector to obtain
  */
  inline
  const RealVect& runtimeVector(const int a_which) const;

  /*!
    @brief Comparison operator with other particle
    @param[in] a_other Other particle
  */
  bool operator==(const ItoParticle& a_other) const;

  /*!
    @brief Comparison operator with other particle
    @param[in] a_other Other particle
  */
  bool operator==(const ItoParticle* a_other) const;

  /*!
    @brief Comparison operator with other particle
    @param[in] a_other Other particle
  */
  bool operator!=(const ItoParticle& a_other) const;

  /*!
    @brief Returns the size, in number of bytes, of a flat representation of the data in this object.
  */
  virtual int size() const override;

  /*!
    @brief Write a linear binary representation of the internal data. Assumes that
    sufficient memory for the buffer has already been allocated by
    the caller.
    @param[in] a_buf Pointer to memory block
  */
  virtual void linearOut(void* buf) const override;

  /*!
    @brief Read a linear binary representation of the internal data.
    Assumes that the buffer has the correct data.
    @param[in] a_buf Pointer to memory block
  */
  virtual void linearIn(void* a_buf);

  /*!
    @brief Friend function for letting ostream print the particle.
    @param[in] ostr       Output stream
    @param[in] a_particle Particle to output
  */
  friend std::ostream & operator<<(std::ostream& ostr, const ItoParticle& a_particle);

protected:

  /*!
    @brief Static for defining the number of runtime Real components
  */
  static int s_numRuntimeScalars;

  /*!
    @brief Static for defining the number of runtime RealVect components
  */
  static int s_numRuntimeVectors;

  /*!
    @brief Temporary memory -- this is not communicated
  */
  Real m_tmp;

  /*!
    @brief Particle mass (i.e., computational weight)
  */
  Real m_mass;

  /*!
    @brief Particle diffusion coefficient
  */
  Real m_diffusion;

  /*!
    @brief Particle average energy
  */  
  Real m_averageEnergy;

  /*!
    @brief Particle mobility coefficient
  */
  Real m_mobility;

  /*!
    @brief Particle velocity
  */
  RealVect m_velocity;

  /*!
    @brief Previous particle position
  */  
  RealVect m_oldPosition;

  /*!
    @brief Pointer to run-time Real memory allocations
  */
  Real* m_runtimeScalars;

  /*!
    @brief Pointer to run-time RealVect memory allocations
  */  
  RealVect* m_runtimeVectors;

  /*!
    @brief Allocate run-time memory blocks. 
  */
  void allocateRuntimeBuffers();
};

/*!
  @brief Particle printing function
  @param[in] ostr       Output stream
  @param[in] a_particle Particle to output
*/
std::ostream & operator<<(std::ostream& ostr, const ItoParticle& p);

#include <CD_NamespaceFooter.H>

#include <CD_ItoParticleImplem.H>

#endif

