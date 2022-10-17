/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCPlasmaReactionState.H
  @brief  Declaration of a simple state vector for running Kinetic Monte Carlo for plasma problems. 
  @author Robert Marskar
*/

#ifndef CD_KMCPlasmaReactionState_H
#define CD_KMCPlasmaReactionState_H

// Std includes
#include <vector>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Declaration of a "plasma state" for advancing with the Kinetic Monte Carlo module. 
  @details This state consists of both "plasma" and "photons". 
  @note The template parameter indicates the integer type used for tracking the states. 
*/
template <typename T = long long>
class KMCPlasmaReactionState
{
public:
  using State = std::vector<T>;

  /*!
    @brief Default constructor. 
  */
  inline KMCPlasmaReactionState() = delete;

  /*!
    @brief Copy constructor
  */
  inline KMCPlasmaReactionState(KMCPlasmaReactionState&) = default;

  /*!
    @brief Disallowed move constructor
  */
  inline KMCPlasmaReactionState(KMCPlasmaReactionState&&) = delete;

  /*!
    @brief Full constructor. Initializes the state with zero particles
  */
  inline KMCPlasmaReactionState(const size_t a_numParticleSpecies, const size_t a_numPhotonSpecies) noexcept;

  /*!
    @brief Destructor
  */
  virtual ~KMCPlasmaReactionState();

  /*!
    @brief Copy assignment constructor
  */
  inline KMCPlasmaReactionState&
  operator=(const KMCPlasmaReactionState&) = default;

  /*!
    @brief Disallowed move assignement
  */
  inline KMCPlasmaReactionState&
  operator=(const KMCPlasmaReactionState&&) = delete;

  /*!
    @brief Check if state is a valid state. An invalid state will have a negative number of particles/photons. 
  */
  inline bool
  isValidState() const noexcept;

  /*!
    @brief Get modifiable particle state
    @return m_particleState
  */
  inline State&
  getParticleState() noexcept;

  /*!
    @brief Get particle state. 
    @return m_particleState
  */
  inline const State&
  getParticleState() const noexcept;

  /*!
    @brief Get modifiable photon state
    @return m_photonState
  */
  inline State&
  getPhotonState() noexcept;

  /*!
    @brief Get photon state
    @return m_photonState
  */
  inline const State&
  getPhotonState() const noexcept;

protected:
  /*!
    @brief Particle state
  */
  State m_particleState;

  /*!
    @brief "Photon" state
  */
  State m_photonState;
};

#include <CD_NamespaceFooter.H>

#include <CD_KMCPlasmaReactionStateImplem.H>

#endif
