/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCPlasmaReaction.H
  @brief  Declaration of a simple plasma reaction type for Kinetic Monte Carlo
  @author Robert Marskar
*/

#ifndef CD_KMCPlasmaReaction_H
#define CD_KMCPlasmaReaction_H

// Std includes
#include <map>
#include <vector>
#include <list>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_KMCPlasmaReactionState.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Reaction type for advancing a KMCPlasmaState for Kinetic Monte Carlo
  @note T is the integer type used for the state. 
*/
template <typename T = long long, typename State = KMCPlasmaReactionState<T>>
class KMCPlasmaReaction
{
public:
  /*!
    @brief Disallowed constructor. Use the full constructor. 
  */
  KMCPlasmaReaction() = delete;

  /*!
    @brief Copy constructor (uses default). 
  */
  KMCPlasmaReaction(const KMCPlasmaReaction&) = default;

  /*!
    @brief Disallowed move constructor
  */
  KMCPlasmaReaction(const KMCPlasmaReaction&&) = delete;

  /*!
    @brief Full constructor
    @param[in] a_particleReactants Particle reactants (particles on left-hand side of reaction)
    @param[in] a_particleProducts  Particle products (particles right-hand side of reactions)
    @param[in] a_photonProducts    Photon products (photons right-hand side of reactions)
  */
  inline KMCPlasmaReaction(const std::list<size_t>& a_particleReactants,
                           const std::list<size_t>& a_particleProducts,
                           const std::list<size_t>& a_photonProducts) noexcept;

  /*!
    @brief Destructor
  */
  inline virtual ~KMCPlasmaReaction();

  /*!
    @brief Get modifiable reaction rate.
    @return m_rate
  */
  inline Real&
  rate() const noexcept;

  /*!
    @brief Get the population of the reactant in the input state.
  */
  inline T
  population(const size_t& a_reactant, const State& a_state) const noexcept;

  /*!
    @brief Compute the propensity function for this reaction type. 
    @note User should set the rate before calling this routine. 
  */
  inline Real
  propensity(const State& a_state) const noexcept;

  /*!
    @brief Compute the number of times the reaction can fire before exhausting one of the reactants.
  */
  inline T
  computeCriticalNumberOfReactions(const State& a_state) const noexcept;

  /*!
    @brief Return the a vector that triggers the propensity function. 
    @details The returned state should be populated with non-zero values for all reactant species. 
  */
  inline std::list<size_t>
  getReactants() const noexcept;

  /*!
    @brief Get the state change due to a change in the input particle reactant
  */
  inline T
  getStateChange(const size_t a_particleReactant) const noexcept;

  /*!
    @brief Advance the incoming state with the number of reactions. 
  */
  inline void
  advanceState(State& a_state, const T& a_numReactions) const noexcept;

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
  std::map<size_t, T> m_particleStateChange;

  /*!
    @brief State change for photons. 
  */
  std::map<size_t, T> m_photonStateChange;

  /*!
    @brief Reaction order
  */
  T m_reactionOrder;

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
  sanityCheck(const State& a_state) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
