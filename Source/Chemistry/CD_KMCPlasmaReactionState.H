/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCPlasmaReactionState.H
  @brief  Declaration of a state vector for running Kinetic Monte Carlo for plasma problems. 
  @author Robert Marskar
*/

#ifndef CD_KMCPlasmaReactionState_H
#define CD_KMCPlasmaReactionState_H

#include <CD_NamespaceHeader.H>

/*!
  @brief Declaration of a "plasma state" for advancing with the Kinetic Monte Carlo module. 
  @details This state consists of both "plasma" and "photons". Plasma can be produced in 
*/
template <typename T = long long>
class KMCPlasmaReactionState
{
public:

  using State = std::vector<T>;

  /*!
    @brief Default constructor. 
  */
  KMCPlasmaReactionState() = delete;

  /*!
    @brief Full constructor. Initializes the state with zero particles
  */
  KMCPlasmaReactionState(const size_t a_numParticleSpecies, const size_t a_numPhotonSpecies) noexcept;

  /*!
    @brief Destructor
  */
  virtual ~KMCPlasmaReactionState();

  /*!
    @brief Get modifiable particle state
    @return m_particleState
  */
  State& getParticleState() noexcept;

  /*!
    @brief Get particle state. 
    @return m_particleState
  */
  const State& getParticleState() const noexcept;  

  /*!
    @brief Get modifiable photon state
    @return m_photonState
  */
  State& getPhotonState() noexcept;

  /*!
    @brief Get photon state
    @return m_photonState
  */
  const State& getPhotonState() const noexcept;
  
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

//#include <CD_KMCPlasmaReactionStateImplem.H>

#endif
