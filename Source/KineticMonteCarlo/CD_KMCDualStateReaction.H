/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCDualStateReaction.H
  @brief  Declaration of a simple plasma reaction type for Kinetic Monte Carlo
  @author Robert Marskar
*/

#ifndef CD_KMCDualStateReaction_H
#define CD_KMCDualStateReaction_H

// Std includes
#include <map>
#include <vector>
#include <list>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_KMCDualState.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Reaction type for advancing a KMCDualState for Kinetic Monte Carlo. 
  @note T is the integer type used for the state. 
*/
template <typename T = long long, typename State = KMCDualState<T>>
class KMCDualStateReaction
{
public:
  /*!
    @brief Disallowed constructor. Use the full constructor. 
  */
  KMCDualStateReaction() = delete;

  /*!
    @brief Copy constructor (uses default). 
  */
  KMCDualStateReaction(const KMCDualStateReaction&) = default;

  /*!
    @brief Disallowed move constructor
  */
  KMCDualStateReaction(const KMCDualStateReaction&&) = delete;

  /*!
    @brief Full constructor
  */
  inline KMCDualStateReaction(const std::list<size_t>& a_lhsReactives,
                              const std::list<size_t>& a_rhsReactives,
                              const std::list<size_t>& a_rhsNonReactives) noexcept;

  /*!
    @brief Destructor
  */
  inline virtual ~KMCDualStateReaction();

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
    @brief Get the reactants in the reaction.
    @return m_lhsReactives
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
  std::list<size_t> m_lhsReactives;

  /*!
    @brief Particle products. 
  */
  std::list<size_t> m_rhsReactives;

  /*!
    @brief Photon products.
  */
  std::list<size_t> m_rhsNonReactives;

  /*!
    @brief State change for particles. 
  */
  std::map<size_t, T> m_reactiveStateChange;

  /*!
    @brief State change for photons. 
  */
  std::map<size_t, T> m_nonReactiveStateChange;

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
