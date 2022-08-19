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

// Our includes
#include <CD_GenericParticle.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief A particle class for use with ItoSolvers, i.e. drifting Brownian walkers.
  @details This class is used to encapsulate the requirements for running an ItoSolver. This computational particle contains position, weight, velocity, as well as a
  diffusion coefficient, a mobility, energy, and previous particle position. These are stored as follows:

  m_scalars[0] => weight
  m_scalars[1] => mobility
  m_scalars[2] => diffusion
  m_scalars[3] => energy
  m_vectors[0] => oldPosition 
  m_vectors[1] => velocity
*/
class ItoParticle : public GenericParticle<4, 2>
{
public:
  /*!
    @brief Default constructor -- user should subsequently set the variables or call define.
    @note If the user has called setNumRuntimeScalars/Vectors then  will allocate runtime buffers
  */
  inline ItoParticle();

  /*!
    @brief Constructor. This calls the define function.
    @param[in] a_weight    Particle weight
    @param[in] a_position  Particle position
    @param[in] a_velocity  Particle velocity
    @param[in] a_diffusion Particle diffusion coefficient
    @param[in] a_mobility  Particle mobility coefficient
    @param[in] a_energy    Particle average energy
    @note If the user has called setNumRuntimeScalars/Vectors then  will allocate runtime buffers
  */
  inline ItoParticle(const Real      a_weight,
                     const RealVect& a_position,
                     const RealVect& a_velocity  = RealVect::Zero,
                     const Real      a_diffusion = 0.0,
                     const Real      a_mobility  = 1.0,
                     const Real      a_energy    = 0.0);

  /*!
    @brief Copy constructor. Copies all fields.
    @param[in] a_other Other particle. 
  */
  inline ItoParticle(const ItoParticle& a_other);

  /*!
    @brief Destructor (deallocates runtime memory storage)
  */
  inline virtual ~ItoParticle();

  /*!
    @brief Full define function. 
    @param[in] a_weight        Particle weight
    @param[in] a_position      Particle position
    @param[in] a_velocity      Particle velocity
    @param[in] a_diffusion     Particle diffusion coefficient
    @param[in] a_mobility      Particle mobility coefficient
    @param[in] a_energy Particle average energy
  */
  inline void
  define(const Real      a_weight,
         const RealVect& a_position,
         const RealVect& a_velocity  = RealVect::Zero,
         const Real      a_diffusion = 0.0,
         const Real      a_mobility  = 1.0,
         const Real      a_energy    = 0.0);

  /*!
    @brief Get particle weight
    @return m_weight
  */
  inline Real&
  weight();

  /*! 
    @brief Get particle weight
    @return m_weight
  */
  inline const Real&
  weight() const;

  /*!
    @brief Get particle diffusion coefficient
    @return m_diffusion
  */
  inline Real&
  diffusion();

  /*! 
    @brief Get particle diffusion coefficient
    @return m_diffusion
  */
  inline const Real&
  diffusion() const;

  /*!
    @brief Get mobility coefficient
    @return m_mobility
  */
  inline Real&
  mobility();

  /*! 
    @brief Get mobility coefficient
    @return m_mobility
  */
  inline const Real&
  mobility() const;

  /*!
    @brief Get average particle energy
    @return m_energy
  */
  inline Real&
  energy();

  /*! 
    @brief Get average particle energy
    @return m_energy
  */
  inline const Real&
  energy() const;

  /*!
    @brief Get the old particle position
    @return m_oldPosition.
  */
  inline RealVect&
  oldPosition();

  /*!
    @brief Get the old particle position
    @return m_oldPosition.
  */
  inline const RealVect&
  oldPosition() const;

  /*!
    @brief Get the particle velocity
    @return m_velocity
  */
  inline RealVect&
  velocity();

  /*!
    @brief Get the particle velocity
    @return m_velocity
  */
  inline const RealVect&
  velocity() const;

  /*!
    @brief Get the total energy
    @return Returns m_weight * m_energy
  */
  inline Real
  totalEnergy() const;

  /*!
    @brief Get the particle conductivity
    @return Returns m_weight * m_mobility
  */
  inline Real
  conductivity() const;

  /*!
    @brief Get the particle diffusivity
    @return Returns m_weight * m_diffusion
  */
  inline Real
  diffusivity() const;

  /*!
    @brief Get temporary storage
    @return m_tmp
  */
  inline Real&
  tmp();

  /*! 
    @brief Get temporary storage
    @return m_tmp
  */
  inline const Real&
  tmp() const;

protected:
  /*!
    @brief Temporary memory -- this is not communicated
  */
  Real m_tmp;
};

#include <CD_NamespaceFooter.H>

#include <CD_ItoParticleImplem.H>

#endif
