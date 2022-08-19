/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
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

// Our includes
#include <CD_GenericParticle.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Particle class for usage with Monte Carlo radiative transfer.
*/
class Photon : public GenericParticle<2, 1>
{
public:
  /*!
    @brief Default constructor. Must subsequently set all parameters
  */
  inline Photon();

  /*!
    @brief Full constructor. 
    @param[in] a_position Photon position
    @param[in] a_velocity Photon velocity
    @param[in] a_kappa    Mean absorption length
    @param[in] a_weight   Particle weight. a_weight=1 implies a physical photon. 
  */
  inline Photon(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_weight = 1.0);

  /*!
    @brief Destructor (does nothing)
  */
  inline virtual ~Photon();

  /*!
    @brief Define function. Sets all parameters.
    @param[in] a_position Photon position
    @param[in] a_velocity Photon velocity
    @param[in] a_kappa    Mean absorption length
    @param[in] a_weight   Particle weight. a_weight=1 implies a physical photon. 
  */
  inline void
  define(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_weight = 1.0);

  /*!
    @brief Get photon weight
  */
  inline Real&
  weight();

  /*!
    @brief Get particle weight
  */
  inline const Real&
  weight() const;

  /*!
    @brief Get photon absorption length
  */
  inline Real&
  kappa();

  /*!
    @brief Get mean absorption length
  */
  inline const Real&
  kappa() const;

  /*!
    @brief Get photon velocity
  */
  inline RealVect&
  velocity();

  /*!
    @brief Get photon velocity
  */
  inline const RealVect&
  velocity() const;
};

#include <CD_NamespaceFooter.H>

#include <CD_PhotonImplem.H>

#endif
