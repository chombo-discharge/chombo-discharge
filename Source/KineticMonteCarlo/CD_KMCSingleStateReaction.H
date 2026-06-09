/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_KMCSingleStateReaction.H
  @brief  Declaration of a simple reaction type for advancing "single states" in Kinetic Monte Carlo codes.
  @author Robert Marskar
*/

#ifndef CD_KMCSINGLESTATEREACTION_H
#define CD_KMCSINGLESTATEREACTION_H

// Std includes
#include <map>
#include <vector>
#include <list>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_KMCSingleState.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Reaction type for advancing a KMCSingleState for Kinetic Monte Carlo.
  @note T is the integer type used for the state.
*/
template <typename State = KMCSingleState<long long>, typename T = long long>
class KMCSingleStateReaction
{
public:
  /**
    @brief Default constructor.
  */
  KMCSingleStateReaction() = default;

  /**
    @brief Copy constructor.
  */
  KMCSingleStateReaction(const KMCSingleStateReaction&) = default;

  /**
    @brief Full constructor, sets reactants and products.
    @details For a reaction X_0 -> X_0 + X_0 we have a_reactants = [0] and a_products = [0,0].
    @param[in] a_reactants Species on the left-hand side of the reaction.
    @param[in] a_products  Species on the right-hand side of the reaction.
  */
  inline KMCSingleStateReaction(const std::list<size_t>& a_reactants, const std::list<size_t>& a_products) noexcept;

  /**
    @brief Destructor.
  */
  inline virtual ~KMCSingleStateReaction();

  /**
    @brief Get modifiable reaction rate.
    @return m_rate
  */
  inline Real&
  rate() const noexcept;

  /**
    @brief Get the population of the reactant in the input state.
    @param[in] a_reactant Reactant.
    @param[in] a_state    State.
    @return Population of the reactant.
  */
  static inline T
  population(const size_t& a_reactant, const State& a_state) noexcept;

  /**
    @brief Compute the propensity function for this reaction type.
    @param[in] a_state State vector.
    @return Reaction propensity.
    @note User should set the rate before calling this routine.
  */
  inline Real
  propensity(const State& a_state) const noexcept;

  /**
    @brief Compute the number of times the reaction can fire before exhausting one of the reactants.
    @param[in] a_state State vector.
    @return Critical number of reactions.
  */
  inline T
  computeCriticalNumberOfReactions(const State& a_state) const noexcept;

  /**
    @brief Get the reactants involved in the reaction.
    @return m_reactants
  */
  inline std::list<size_t>
  getReactants() const noexcept;

  /**
    @brief Get the state change due to a change in the input reactant species.
    @param[in] a_particleReactant Reactant species.
    @return State change for the reactant.
  */
  inline T
  getStateChange(const size_t a_particleReactant) const noexcept;

  /**
    @brief Advance the incoming state with the number of reactions.
    @param[in,out] a_state        State vector.
    @param[in]    a_numReactions Number of reactions.
  */
  inline void
  advanceState(State& a_state, const T& a_numReactions) const noexcept;

protected:
  /**
    @brief Reaction rate.
  */
  mutable Real m_rate;

  /**
    @brief Factor for the case where two or more particles of the same species react.
    @details For reactions X + X -> null the propensity is k * X * (X-1)/2 because there are X * (X-1)/2 unique pairs of particles in the cell.
  */
  Real m_propensityFactor;

  /**
    @brief Reactants.
  */
  std::list<size_t> m_reactants;

  /**
    @brief Products.
  */
  std::list<size_t> m_products;

  /**
    @brief State change for reactants/products.
  */
  std::map<size_t, T> m_stateChange;

  /**
    @brief Compute state change vectors from the reactant/product lists.
  */
  inline void
  computeStateChanges() noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_KMCSingleStateReactionImplem.H>

#endif
