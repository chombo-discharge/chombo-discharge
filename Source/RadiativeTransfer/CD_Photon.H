/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Photon.H
  @brief  Declaration of a photon class for particle methods
  @author Robert Marskar
*/

#ifndef CD_Photon_H
#define CD_Photon_H

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
  @brief Photon particle class. This inherits from Chombo's BinItem so it can be used in parallel with our ParticleContainer classes.
*/
class Photon : public BinItem {
public:

  /*!
    @brief Default constructor. Must subsequently set all parameters
  */
  Photon();

  /*!
    @brief Full constructor. 
    @param[in] a_position Photon position
    @param[in] a_velocity Photon velocity
    @param[in] a_kappa    Mean absorption length
    @param[in] a_mass     Particle weight. a_mass=1 implies a physical photon. 
  */  
  Photon(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_mass = 1.0);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~Photon();

  /*!
    @brief Define function. Sets all parameters.
    @param[in] a_position Photon position
    @param[in] a_velocity Photon velocity
    @param[in] a_kappa    Mean absorption length
    @param[in] a_mass     Particle weight. a_mass=1 implies a physical photon. 
  */
  void define(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_mass = 1.0);

  /*!
    @brief Set mean absorption length
    @param[in] a_kappa Mean absorption length
  */
  void setKappa(const Real a_kappa);

  /*!
    @brief Set particle mass
    @param[in] a_mass Particle weight. a_mass=1 implies a physical photon. 
  */
  void setMass(const Real a_mass);

  /*!
    @brief Set photon velocity
    @param[in] a_velocity Photon velocity
  */
  void setVelocity(const RealVect& a_velocity);

  /*!
    @brief Get photon mass
  */
  Real& mass();

  /*!
    @brief Get photon absorption length
  */  
  Real& kappa();

  /*!
    @brief Get photon velocity
  */    
  RealVect& velocity();

  /*!
    @brief Get photon mass
  */    
  const Real& mass()   const;  

  /*!
    @brief Get mean absorption length
  */  
  const Real& kappa()    const;

  /*!
    @brief Get photon velocity
  */
  const RealVect& velocity() const;

  /*!
    @brief Comparison operator
    @param[in] a_p Other photon
  */  
  virtual bool operator==(const Photon& a_p) const;
  
  /*!
    @brief Comparison operator
    @param[in] a_p Other photon
  */    
  virtual bool operator==(const Photon* a_p) const;

  /*!
    @brief Comparison operator
    @param[in] a_p Other photon
  */      
  virtual bool operator!=(const Photon& a_p) const;

  /*!
    @brief MPI linearization function, returns the size of the object (this is BinItem + m_kappa, m_mass, m_velocity)
  */
  virtual int size() const;

  /*!
    @brief MPI linearization function. Writes this object into the buffer
  */  
  virtual void linearOut(void* buf) const;

  /*!
    @brief MPI linearization function. Takes the input buffer and writes into this object. 
  */    
  virtual void linearIn(void* buf);

protected:

  /*!
    @brief Mean absorption length
  */
  Real m_kappa;
  
  /*!
    @brief Photon mass (computational weight)
  */  
  Real m_mass;

  /*!
    @brief Photon velocity
  */
  RealVect m_velocity;
};

#include <CD_NamespaceFooter.H>

#endif 
