/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KMCSolver.H
  @brief  Class for running Kinetic Monte Carlo functionality. 
  @author Robert Marskar
*/

#ifndef CD_KMCSolver_H
#define CD_KMCSolver_H

// Std includes
#include <vector>
#include <memory>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for running Kinetic Monte-Carlo simulations. 
  @details The template parameter State is the underlying state type that KMC operators on. There are no required member functions on the State parameter,
  but the reaction type (template parameter R) MUST be able to operate on the state through the following functions:

  1. Real R::propensity(State) const -> Computes the reaction propensity for the input state.
  2. T R::computeCriticalNumberOfReactions(State) const -> Computes the minimum number of reactions that exhausts one of the reactants. 
  3. void R::advanceState(State&, const T numReactions) const -> Advance state by numReactions
  4. std::<some_container> getReactants() const -> Get reactants involved in the reactions.
  5. T R::population(const <some_type> reactant, const State& a_state) -> Get the population of the input reactant in the input state. 

  The template parameter T should agree across both both R, State, and KMCSolver. 
*/
template <typename R, typename State, typename T = long long>
class KMCSolver
{
public:
  using ReactionList = std::vector<std::shared_ptr<const R>>;

  /*!
    @brief Disallowed constructor.
  */
  KMCSolver() = delete;

  /*!
    @brief Disallowed copy constructor
  */
  KMCSolver(const KMCSolver&) = default;

  /*!
    @brief Disallowed move constructor
  */
  KMCSolver(const KMCSolver&&) = delete;

  /*!
    @brief Full constructor. 
  */
  inline KMCSolver(const ReactionList& a_reaction) noexcept;

  /*!
    @brief Destructor
  */
  virtual ~KMCSolver() noexcept;

  /*!
    @brief Copy assignment operator. 
  */
  KMCSolver&
  operator=(const KMCSolver&) = default;

  /*!
    @brief Disallowed move assignment operator. 
  */
  KMCSolver&
  operator=(const KMCSolver&&) = delete;

  /*!
    @brief Set solver parameters
    @param[in] a_numCrit Determines critical reactions. This is the number of reactions that need to fire before depleting a reactant.
    @param[in] a_numSSA  Maximum number of SSA steps to run when switching from tau-leaping to SSA (hybrid algorithm only). 
    @param[in] a_eps     Maximum permitted change in propensities when performing tau-leaping for non-critical reactions.
    @param[in] a_SSAlim  Threshold for switching from tau-leaping of non-critical reactions to SSA for all reactions (hybrid algorithm only)
  */
  inline void
  setSolverParameters(const T a_numCrit, const T a_numSSA, const Real a_eps, const Real a_SSAlim) noexcept;

  /*!
    @brief Compute propensities for ALL reactions
  */
  inline std::vector<Real>
  propensities(const State& a_state) const noexcept;

  /*!
    @brief Compute propensities for a subset of reactions
  */
  inline std::vector<Real>
  propensities(const State& a_state, const ReactionList& a_reactions) const noexcept;

  /*!
    @brief Compute the total propensity for ALL reactions
  */
  inline Real
  totalPropensity(const State& a_state) const noexcept;

  /*!
    @brief Compute the total propensity for a subset of reactions
  */
  inline Real
  totalPropensity(const State& a_state, const ReactionList& a_reactions) const noexcept;

  /*!
    @brief Partition reactions into critical and non-critical reactions. 
    @details First member in the pair is the critical reactions, then the non-critical reactions.
    @note Calls the other version with m_reactions
  */
  inline std::pair<ReactionList, ReactionList>
  partitionCriticalReactions(const State& a_state) const noexcept;

  /*!
    @brief Partition reactions into critical and non-critical reactions. 
    @details First member in the pair is the critical reactions, then the non-critical reactions
  */
  inline std::pair<ReactionList, ReactionList>
  partitionCriticalReactions(const State& a_state, const ReactionList& a_reactions) const noexcept;

  /*!
    @brief Get the time to the next critical reaction
    @note Computes the total propensity and calls the other version. 
  */
  inline Real
  getCriticalTimeStep(const State& a_state, const ReactionList& a_criticalReactions) const noexcept;

  /*!
    @brief Get the time to the next critical reaction
    @note Computes the total propensity and calls the other version. 
  */
  inline Real
  getCriticalTimeStep(const State& a_state, const std::vector<Real>& a_propensities) const noexcept;

  /*!
    @brief Get the time to the next critical reaction
    @note Computes the total propensity and calls the other version. 
  */
  inline Real
  getCriticalTimeStep(const Real& a_totalPropensity) const noexcept;

  /*!
    @brief Get the time to the next critical reaction
  */
  inline Real
  getNonCriticalTimeStep(const State&             a_state,
                         const ReactionList&      a_nonCriticalReactions,
                         const std::vector<Real>& a_nonCriticalPropensities) const noexcept;

  /*!
    @brief Perform one tau-leaping step using ALL reactions. 
    @note Calls the other version with m_reactions
  */
  inline void
  stepTau(State& a_state, const Real a_dt) const noexcept;

  /*!
    @brief Perform one tau-leaping step over the input reactions using a time step a_dt
    @note Computes propensities and calls the other version. 
  */
  inline void
  stepTau(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept;

  /*!
    @brief Perform one tau-leaping step over the input reactions using a time step a_dt
  */
  inline void
  stepTau(State&                   a_state,
          const ReactionList&      a_reactions,
          const std::vector<Real>& a_propensities,
          const Real               a_dt) const noexcept;

  /*!
    @brief Advance with tau-leaping step over the input time. This can end up using substepping.
    @note Calls the other version with m_reactions.
  */
  inline void
  advanceTau(State& a_state, const Real a_dt) const noexcept;

  /*!
    @brief Advance with tau-leaping step over the input time. This can end up using substepping
  */
  inline void
  advanceTau(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept;

  /*!
    @brief Perform a single SSA step.
    @note Calls the other version with m_reactions
  */
  inline void
  stepSSA(State& a_state) const noexcept;

  /*!
    @brief Perform a single SSA step.
  */
  inline void
  stepSSA(State& a_state, const ReactionList& a_reactions) const noexcept;

  /*!
    @brief Perform a single SSA step. This version has pre-computed propensities (for optimization reasons)
  */
  inline void
  stepSSA(State& a_state, const ReactionList& a_reactions, const std::vector<Real>& a_propensities) const noexcept;

  /*!
    @brief Advance with the SSA over the input time. This can end up using substepping
    @note Calls the other version with m_reactions
  */
  inline void
  advanceSSA(State& a_state, const Real a_dt) const noexcept;

  /*!
    @brief Advance with the SSA over the input time. This can end up using substepping
  */
  inline void
  advanceSSA(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept;

  /*!
    @brief Advance using hybrid method over the input time. This can end up using substepping.
    @note Calls the other version with m_reactions
  */
  inline void
  advanceHybrid(State& a_state, const Real a_dt) const noexcept;

  /*!
    @brief Advance using hybrid method over the input time. This can end up using substepping.
  */
  inline void
  advanceHybrid(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept;

protected:
  /*!
    @brief List of reactions used when advancing states. 
  */
  ReactionList m_reactions;

  /*!
    @brief Definition of critical reactions. 
    @details A reaction is critical if it is m_Ncrit firings away from depleting a reactant. 
  */
  T m_Ncrit;

  /*!
    @brief Maximum number of SSA steps to run when switching into SSA-based advancement for non-critical reactions. 
  */
  T m_numSSA;

  /*!
    @brief
  */
  Real m_eps;

  /*!
    @brief
  */
  Real m_SSAlim;
};

#include <CD_NamespaceFooter.H>

#include <CD_KMCSolverImplem.H>

#endif
