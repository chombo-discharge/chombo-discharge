/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PlasmaReaction.H
  @brief  Declaration of a plasma reaction type. 
  @author Robert Marskar
*/

#ifndef CD_PlasmaReaction_H
#define CD_PlasmaReaction_H

// Std includes
#include <map>
#include <vector>
#include <list>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Plasma class for encapsulating stoichiometric reactive processes of the type A + B -> C with rate k. 
  @details This reaction type is fairly raw; it simply associates indices in some C++ random-access container
  with the left hand and right hand side of a reaction. The user can choose if the state represents a density or 
  particles by changing the underlying template parameter T. User can also change the underlying container type. 
*/
template <typename T = long long, typename State = std::vector<T>>
class PlasmaReaction
{
public:
  /*!
    @brief Disallowed constructor. Use the full constructor. 
  */
  PlasmaReaction() = delete;

  /*!
    @brief Copy constructor (uses default). 
  */
  PlasmaReaction(const PlasmaReaction&) = default;

  /*!
    @brief Disallowed move constructor
  */
  PlasmaReaction(const PlasmaReaction&&) = delete;

  /*!
    @brief Full constructor
    @param[in] a_particleReactants Particle reactants (left-hand side of reaction)
    @param[in] a_particleProducts  Particle products (right-hand side of reactions)
    @param[in] a_photonProducts    Photon products (right-hand side of reactions)
  */
  inline PlasmaReaction(const std::list<size_t>& a_particleReactants,
                        const std::list<size_t>& a_particleProducts,
                        const std::list<size_t>& a_photonProducts) noexcept;

  /*!
    @brief Destructor
  */
  inline
  virtual ~PlasmaReaction();

  /*!
    @brief Copy assignment operator
  */
  PlasmaReaction&
  operator=(const PlasmaReaction&) = default;

  /*!
    @brief Disallowed move operator
  */
  PlasmaReaction&
  operator=(const PlasmaReaction&&) = delete;

  /*!
    @brief Get modifiable reaction rate.
    @return m_rate
  */
  inline Real&
  rate() const noexcept;

  /*!
    @brief Compute the propensity function for this reaction type. 
    @details This computes k * a_trackedState[r1] * ... * a_untrackedState[r2] * .... 
    @note User should set k before calling this routine. 
  */
  inline Real
  propensity(const State& a_particleState) const noexcept;

  /*!
    @brief Advance the incoming state with the number of reactions. 
  */
  inline void
  advanceState(State& a_particleState, State& a_photonState, const T& a_numReactions) const noexcept;

  /*!
    @brief Get the particle state change.
    @return m_particleStateChange;
  */
  inline const std::map<size_t, int>&
  getParticleStateChange() const noexcept;

  /*!
    @brief Get the photon state change.
    @return m_photonStateChange;
  */
  inline const std::map<size_t, int>&
  getPhotonStateChange() const noexcept;

  /*!
    @brief Get the particle reactants
    @return Returns m_particleReactan
  */
  inline const std::list<size_t>&
  getParticleReactants() const noexcept;

  /*!
    @brief Get the untracked Reactants
    @return Returns m_untrackedReactants
  */
  inline const std::list<size_t>&
  getParticleProducts() const noexcept;

  /*!
    @brief Get the plasma products
    @return Returns m_plasmaProducts
  */
  inline const std::list<size_t>&
  getPhotonProducts() const noexcept;

protected:
  
  /*!
    @brief Reaction rate
  */
  mutable Real m_rate;

  /*!
    @brief Particle reactants. 
  */
  std::list<size_t> m_particleReactants;

  /*!
    @brief Particle products. 
  */
  std::list<size_t> m_particleProducts;

  /*!
    @brief Photon products.
  */
  std::list<size_t> m_photonProducts;

  /*!
    @brief State change for particles. 
  */
  std::map<size_t, int> m_particleStateChange;

  /*!
    @brief State change for photons. 
  */
  std::map<size_t, int> m_photonStateChange;

  /*!
    @brief Compute state change
  */
  inline void
  computeStateChanges() noexcept;

  /*!
    @brief Debugging function which ensures that the class data holders do not reach out of the incoming state. 
    @details This is necessary because PlasmaReaction does not have compile-time size restrictions on the incoming state. The
    internals of this method uses Chombo assertions so the compiler should be able to take this function out of the executable
    if compiling with the correct flags.
  */
  inline void
  sanityCheckParticles(const State& a_particleState) const noexcept;

  /*!
    @brief Debugging function which ensures that the class data holders do not reach out of the incoming state. 
    @details This is necessary because PlasmaReaction does not have compile-time size restrictions on the incoming state. The
    internals of this method uses Chombo assertions so the compiler should be able to take this function out of the executable
    if compiling with the correct flags.
  */  
  inline void
  sanityCheckPhotons(const State& a_photonState) const noexcept;  
};

#include <CD_NamespaceFooter.H>

#include <CD_PlasmaReactionImplem.H>

#endif
