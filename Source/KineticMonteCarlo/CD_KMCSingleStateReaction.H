/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCSingleStateReaction.H
  @brief  Declaration of a simple reaction type for advancing "single states" in Kinetic Monte Carlo codes. 
  @author Robert Marskar
*/

#ifndef CD_KMCSingleStateReaction_H
#define CD_KMCSingleStateReaction_H

// Std includes
#include <map>
#include <vector>
#include <list>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_KMCSingleState.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Reaction type for advancing a KMCSingleState for Kinetic Monte Carlo. 
  @note T is the integer type used for the state. 
*/
template <typename T = long long, typename State = KMCSingleState<T>>
class KMCSingleStateReaction
{
public:
  /*!
    @brief Disallowed constructor. Use the full constructor. 
  */
  KMCSingleStateReaction() = default;

  /*!
    @brief Copy constructor (uses default). 
  */
  KMCSingleStateReaction(const KMCSingleStateReaction&) = default;

  /*!
    @brief Full constructor, sets reactants and products.
    @details For a reaction X_0 -> X_0 + X_0 we have a_reactants = [0] and a_products [0,0]
    @param[in] a_reactants Species on the left-hand side of the reaction. 
    @param[in] a_products  Species on the right-hand side of the reaction. 
  */
  inline KMCSingleStateReaction(const std::list<size_t>& a_reactants,
				const std::list<size_t>& a_products) noexcept;

  /*!
    @brief Destructor
  */
  inline virtual ~KMCSingleStateReaction();

  /*!
    @brief Get modifiable reaction rate.
    @return m_rate
  */
  inline Real&
  rate() const noexcept;

  /*!
    @brief Get the population of the reactant in the input state.
    @param[in] a_reactant Reactant
    @param[in] a_state    State
  */
  inline T
  population(const size_t& a_reactant, const State& a_state) const noexcept;

  /*!
    @brief Compute the propensity function for this reaction type. 
    @param[in] a_state State vector
    @note User should set the rate before calling this routine. 
  */
  inline Real
  propensity(const State& a_state) const noexcept;

  /*!
    @brief Compute the number of times the reaction can fire before exhausting one of the reactants.
    @param[in] a_state State vector
  */
  inline T
  computeCriticalNumberOfReactions(const State& a_state) const noexcept;

  /*!
    @brief Get the reactants involved in the reaction.
    @return m_reactants
  */
  inline std::list<size_t>
  getReactants() const noexcept;

  /*!
    @brief Get the state change due to a change in the input reactant species.
    @param[in] a_rectant Reactant species. 
  */
  inline T
  getStateChange(const size_t a_reactant) const noexcept;

  /*!
    @brief Advance the incoming state with the number of reactions. 
    @param[in] a_state        State vector
    @param[in] a_numReactions Number of reactions.
  */
  inline void
  advanceState(State& a_state, const T& a_numReactions) const noexcept;

protected:
  /*!
    @brief Reaction rate
  */
  mutable Real m_rate;

  /*!
    @brief Reactants
  */
  std::list<size_t> m_reactants;

  /*!
    @brief Products
  */
  std::list<size_t> m_products;

  /*!
    @brief State change for reactants/products
  */
  std::map<size_t, T> m_stateChange;

  /*!
    @brief Compute state change
  */
  inline void
  computeStateChanges() noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_KMCSingleStateReactionImplem.H>

#endif
