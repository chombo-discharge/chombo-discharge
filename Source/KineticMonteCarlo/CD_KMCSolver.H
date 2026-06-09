/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_KMCSolver.H
  @brief  Class for running Kinetic Monte Carlo functionality.
  @author Robert Marskar
*/

#ifndef CD_KMCSOLVER_H
#define CD_KMCSOLVER_H

// Std includes
#include <vector>
#include <memory>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Supported propagators for hybrid tau leaping.
  @details These leap propagators are as follows:
  ExplicitEuler = Regular tau leaping.
  Midpoint = Gillespie's midpoint method.
  PRC = Hu and Li's Poisson random correction leap method.
*/
enum class KMCLeapPropagator
{
  ExplicitEuler, ///< Regular tau leaping.
  Midpoint,      ///< Gillespie's midpoint method.
  PRC,           ///< Hu and Li's Poisson random correction method.
  ImplicitEuler  ///< Implicit Euler tau leaping.
};

/**
  @brief Class for running Kinetic Monte-Carlo simulations.
  @details The template parameter State is the underlying state type that KMC operators on. There are rather simple
  required member functions on the State parameter, but the reaction type (template parameter R) MUST be able to
  operate on the state through the following functions:

  1. Real R::propensity(State) const -> Computes the reaction propensity for the input state.
  2. T R::computeCriticalNumberOfReactions(State) const -> Computes the minimum number of reactions that exhausts one of the reactants.
  3. void R::advanceState(State&, const T numReactions) const -> Advance state by numReactions.
  4. std::&lt;some_container&gt; getReactants() const -> Get reactants involved in the reactions.
  5. static T R::population(const &lt;some_type&gt; reactant, const State& a_state) -> Get the population of the input reactant in the input state.

  The template parameter T should agree across both R, State, and KMCSolver. Note that this must be a signed integer type.
*/
template <typename R, typename State, typename T = long long>
class KMCSolver
{
public:
  /**
    @brief Alias for the list of reactions.
  */
  using ReactionList = std::vector<std::shared_ptr<const R>>;

  /**
    @brief Default constructor. Must subsequently call define.
  */
  KMCSolver() noexcept;

  /**
    @brief Copy constructor.
  */
  KMCSolver(const KMCSolver&) = default;

  /**
    @brief Disallowed move constructor.
  */
  KMCSolver(const KMCSolver&&) = delete;

  /**
    @brief Full constructor.
    @param[in] a_reactions List of reactions.
  */
  inline KMCSolver(const ReactionList& a_reactions) noexcept;

  /**
    @brief Destructor.
  */
  virtual ~KMCSolver() noexcept;

  /**
    @brief Copy assignment operator.
    @param[in] a_other Other solver.
    @return Reference to this object.
  */
  KMCSolver&
  operator=(const KMCSolver& a_other) = default;

  /**
    @brief Disallowed move assignment operator.
  */
  KMCSolver&
  operator=(const KMCSolver&&) = delete;

  /**
    @brief Define function. Sets the reactions.
    @param[in] a_reactions List of reactions.
  */
  inline void
  define(const ReactionList& a_reactions) noexcept;

  /**
    @brief Set solver parameters.
    @param[in] a_numCrit Determines critical reactions. This is the number of reactions that need to fire before depleting a reactant.
    @param[in] a_numSSA  Maximum number of SSA steps to run when switching from tau-leaping to SSA (hybrid algorithm only).
    @param[in] a_maxIter Maximum permitted number of iterations for the implicit solver.
    @param[in] a_eps     Maximum permitted change in propensities when performing tau-leaping for non-critical reactions.
    @param[in] a_SSAlim  Threshold for switching from tau-leaping of non-critical reactions to SSA for all reactions (hybrid algorithm only).
    @param[in] a_exitTol Exit tolerance for implicit solvers.
  */
  inline void
  setSolverParameters(T a_numCrit, T a_numSSA, T a_maxIter, Real a_eps, Real a_SSAlim, Real a_exitTol) noexcept;

  /**
    @brief Compute the state vector changes for all reactions.
    @param[in] a_state     Input state.
    @param[in] a_reactions Reaction list.
    @return State change vectors for each reaction.
    @note The input parameter a_state can be a dummy state because the state change vector does not depend on it.
  */
  inline std::vector<std::vector<T>>
  getNu(const State& a_state, const ReactionList& a_reactions) const noexcept;

  /**
    @brief Compute propensities for ALL reactions.
    @param[in] a_state State vector.
    @return Propensities for each reaction.
    @note Calls the other version with m_reactions.
  */
  inline std::vector<Real>
  propensities(const State& a_state) const noexcept;

  /**
    @brief Compute propensities for a subset of reactions.
    @param[in] a_state     State vector.
    @param[in] a_reactions Reaction list.
    @return Propensities for each reaction.
  */
  inline std::vector<Real>
  propensities(const State& a_state, const ReactionList& a_reactions) const noexcept;

  /**
    @brief Compute the total propensity for ALL reactions.
    @param[in] a_state State vector.
    @return Sum of all reaction propensities.
    @note Calls the other version with m_reactions.
  */
  inline Real
  totalPropensity(const State& a_state) const noexcept;

  /**
    @brief Compute the total propensity for a subset of reactions.
    @param[in] a_state     State vector.
    @param[in] a_reactions Reaction list.
    @return Sum of propensities for the given reactions.
  */
  inline Real
  totalPropensity(const State& a_state, const ReactionList& a_reactions) const noexcept;

  /**
    @brief Partition reactions into critical and non-critical reactions.
    @details First member in the pair is the critical reactions, then the non-critical reactions.
    @param[in] a_state State vector.
    @return Pair of (critical reactions, non-critical reactions).
    @note Calls the other version with m_reactions.
  */
  inline std::pair<ReactionList, ReactionList>
  partitionReactions(const State& a_state) const noexcept;

  /**
    @brief Partition reactions into critical and non-critical reactions.
    @details First member in the pair is the critical reactions, then the non-critical reactions.
    @param[in] a_state     State vector.
    @param[in] a_reactions Reaction list to be partitioned.
    @return Pair of (critical reactions, non-critical reactions).
  */
  inline std::pair<ReactionList, ReactionList>
  partitionReactions(const State& a_state, const ReactionList& a_reactions) const noexcept;

  /**
    @brief Get the time to the next critical reaction.
    @param[in] a_state State vector.
    @return Time step for the next critical reaction.
    @note Calls the other version with m_reactions.
  */
  inline Real
  getCriticalTimeStep(const State& a_state) const noexcept;

  /**
    @brief Get the time to the next critical reaction.
    @param[in] a_state             State vector.
    @param[in] a_criticalReactions Reaction list.
    @return Time step for the next critical reaction.
    @note Computes the total propensity and calls the other version.
  */
  inline Real
  getCriticalTimeStep(const State& a_state, const ReactionList& a_criticalReactions) const noexcept;

  /**
    @brief Get the time to the next critical reaction.
    @param[in] a_propensities Reaction propensities.
    @return Time step for the next critical reaction.
  */
  inline Real
  getCriticalTimeStep(const std::vector<Real>& a_propensities) const noexcept;

  /**
    @brief Get the time to the next critical reaction.
    @param[in] a_totalPropensity Total propensity.
    @return Time step for the next critical reaction.
  */
  inline Real
  getCriticalTimeStep(const Real& a_totalPropensity) const noexcept;

  /**
    @brief Get the non-critical time step.
    @param[in] a_state State vector.
    @return Non-critical tau-leaping time step.
    @note Calls the other version with m_reactions.
  */
  inline Real
  getNonCriticalTimeStep(const State& a_state) const noexcept;

  /**
    @brief Get the non-critical time step.
    @param[in] a_state     State vector.
    @param[in] a_reactions Reaction list.
    @return Non-critical tau-leaping time step.
    @note Computes propensities and calls the other version.
  */
  inline Real
  getNonCriticalTimeStep(const State& a_state, const ReactionList& a_reactions) const noexcept;

  /**
    @brief Get the non-critical time step.
    @param[in] a_state                   State vector.
    @param[in] a_nonCriticalReactions    Non-critical reactions.
    @param[in] a_nonCriticalPropensities Non-critical propensities.
    @return Non-critical tau-leaping time step.
  */
  inline Real
  getNonCriticalTimeStep(const State&             a_state,
                         const ReactionList&      a_nonCriticalReactions,
                         const std::vector<Real>& a_nonCriticalPropensities) const noexcept;

  /**
    @brief Compute a time step using the leap condition on the mean value.
    @details This is just like getNonCriticalTimeStep EXCEPT that the leap condition is not applied for the deviation.
    @param[in] a_state        State vector.
    @param[in] a_reactions    Reactions.
    @param[in] a_propensities Propensities.
    @param[in] a_epsilon      Epsilon value for maximum change in propensities.
    @return Tau-leaping time step.
  */
  inline Real
  computeDt(const State&             a_state,
            const ReactionList&      a_reactions,
            const std::vector<Real>& a_propensities,
            Real                     a_epsilon) const noexcept;

  /**
    @brief Perform a single SSA step.
    @param[in,out] a_state State vector to advance.
    @note Calls the other version with m_reactions.
  */
  inline void
  stepSSA(State& a_state) const noexcept;

  /**
    @brief Perform a single SSA step.
    @param[in,out] a_state     State vector to advance.
    @param[in]    a_reactions Reactions to advance with.
    @note Computes propensities and calls the other version.
  */
  inline void
  stepSSA(State& a_state, const ReactionList& a_reactions) const noexcept;

  /**
    @brief Perform a single SSA step with pre-computed propensities.
    @param[in,out] a_state        State vector to advance.
    @param[in]    a_reactions    Reactions to advance with.
    @param[in]    a_propensities Propensities for the reactions.
  */
  inline void
  stepSSA(State& a_state, const ReactionList& a_reactions, const std::vector<Real>& a_propensities) const noexcept;

  /**
    @brief Advance with the SSA over the input time. This can end up using substepping.
    @param[in,out] a_state State vector to advance.
    @param[in]    a_dt    Time increment.
    @note Calls the other version with m_reactions.
  */
  inline void
  advanceSSA(State& a_state, const Real a_dt) const noexcept;

  /**
    @brief Advance with the SSA over the input time. This can end up using substepping.
    @param[in,out] a_state     State vector to advance.
    @param[in]    a_reactions Reactions to advance with.
    @param[in]    a_dt        Time increment.
  */
  inline void
  advanceSSA(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept;

  /**
    @brief Perform one plain tau-leaping step using ALL reactions.
    @param[in,out] a_state State vector to be advanced.
    @param[in]    a_dt    Time increment.
    @note Calls the other version with m_reactions.
  */
  inline void
  stepExplicitEuler(State& a_state, const Real a_dt) const noexcept;

  /**
    @brief Perform one plain tau-leaping step over the input reactions.
    @param[in,out] a_state     State vector to be advanced.
    @param[in]    a_reactions List of reactions to advance with.
    @param[in]    a_dt        Time increment.
  */
  inline void
  stepExplicitEuler(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept;

  /**
    @brief Perform one leaping step using the midpoint method for ALL reactions.
    @param[in,out] a_state State vector to be advanced.
    @param[in]    a_dt    Time increment.
    @note Calls the other version with m_reactions.
  */
  inline void
  stepMidpoint(State& a_state, const Real a_dt) const noexcept;

  /**
    @brief Perform one leaping step using the midpoint method for the input reactions.
    @param[in,out] a_state     State vector to be advanced.
    @param[in]    a_reactions List of reactions to advance with.
    @param[in]    a_dt        Time increment.
  */
  inline void
  stepMidpoint(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept;

  /**
    @brief Perform one leaping step using the PRC method for ALL reactions.
    @param[in,out] a_state State vector to be advanced.
    @param[in]    a_dt    Time increment.
    @note Calls the other version with m_reactions.
  */
  inline void
  stepPRC(State& a_state, const Real a_dt) const noexcept;

  /**
    @brief Perform one leaping step using the PRC method for the input reactions.
    @param[in,out] a_state     State vector to be advanced.
    @param[in]    a_reactions List of reactions to advance with.
    @param[in]    a_dt        Time increment.
  */
  inline void
  stepPRC(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept;

  /**
    @brief Perform one implicit Euler tau-leaping step using ALL reactions.
    @param[in,out] a_state State vector to be advanced.
    @param[in]    a_dt    Time increment.
    @note Calls the other version with m_reactions.
  */
  inline void
  stepImplicitEuler(State& a_state, const Real a_dt) const noexcept;

  /**
    @brief Perform one implicit Euler tau-leaping step over the input reactions.
    @param[in,out] a_state     State vector to be advanced.
    @param[in]    a_reactions List of reactions to advance with.
    @param[in]    a_dt        Time increment.
  */
  inline void
  stepImplicitEuler(State& a_state, const ReactionList& a_reactions, const Real a_dt) const noexcept;

  /**
    @brief Advance using a specified tau-leaping algorithm.
    @param[in,out] a_state          State vector to advance.
    @param[in]    a_dt             Time increment.
    @param[in]    a_leapPropagator Which leap propagator to use.
    @note Calls the other version with m_reactions.
  */
  inline void
  advanceTau(State&                   a_state,
             const Real&              a_dt,
             const KMCLeapPropagator& a_leapPropagator = KMCLeapPropagator::ExplicitEuler) const noexcept;

  /**
    @brief Advance using a specified tau-leaping algorithm.
    @param[in,out] a_state          State vector to advance.
    @param[in]    a_reactions      Reactions to advance.
    @param[in]    a_dt             Time increment.
    @param[in]    a_leapPropagator Which leap propagator to use.
    @note Calls the other version with stepTau as the leap propagator.
  */
  inline void
  advanceTau(State&                   a_state,
             const ReactionList&      a_reactions,
             const Real&              a_dt,
             const KMCLeapPropagator& a_leapPropagator = KMCLeapPropagator::ExplicitEuler) const noexcept;

  /**
    @brief Advance using Cao et. al. hybrid algorithm over the input time. This can end up using substepping.
    @param[in,out] a_state          State vector to advance.
    @param[in]    a_dt             Time increment.
    @param[in]    a_leapPropagator Which leap propagator to use.
    @note Calls the other version with m_reactions.
  */
  inline void
  advanceHybrid(State&                   a_state,
                const Real               a_dt,
                const KMCLeapPropagator& a_leapPropagator = KMCLeapPropagator::ExplicitEuler) const noexcept;

  /**
    @brief Advance using Cao et. al. hybrid algorithm over the input time. This can end up using substepping.
    @param[in,out] a_state          State vector to advance.
    @param[in]    a_reactions      Reactions to advance.
    @param[in]    a_dt             Time increment.
    @param[in]    a_leapPropagator Which leap propagator to use.
    @note Calls the other version with stepTau as the leap propagator.
  */
  inline void
  advanceHybrid(State&                   a_state,
                const ReactionList&      a_reactions,
                const Real               a_dt,
                const KMCLeapPropagator& a_leapPropagator = KMCLeapPropagator::ExplicitEuler) const noexcept;

  /**
    @brief Advance using Cao et. al. hybrid algorithm over the input time. This can end up using substepping.
    @param[in,out] a_state      State vector to advance.
    @param[in]    a_reactions  Reactions to advance with.
    @param[in]    a_dt         Time increment.
    @param[in]    a_propagator Leaping propagator function.
  */
  inline void
  advanceHybrid(
    State&                                                                               a_state,
    const ReactionList&                                                                  a_reactions,
    const Real                                                                           a_dt,
    const std::function<void(State&, const ReactionList& a_reactions, const Real a_dt)>& a_propagator) const noexcept;

protected:
  /**
    @brief List of reactions used when advancing states.
  */
  ReactionList m_reactions;

  /**
    @brief Definition of critical reactions.
    @details A reaction is critical if it is m_Ncrit firings away from depleting a reactant.
  */
  T m_Ncrit;

  /**
    @brief Maximum number of SSA steps to run when switching into SSA-based advancement for non-critical reactions.
  */
  T m_numSSA;

  /**
    @brief Maximum number of iterations for implicit solvers.
  */
  T m_maxIter;

  /**
    @brief Maximum permitted change in propensities for non-critical reactions.
  */
  Real m_eps;

  /**
    @brief Threshold for switching to SSA-based algorithm within the Cao algorithm.
  */
  Real m_SSAlim;

  /**
    @brief Exit tolerance for implicit solvers.
  */
  Real m_exitTol;

  /**
    @brief List of state changes for each reaction.
    @details The outer vector corresponds to the reactive index. The inner vector describes the change in the population of the
    state vector, with an indexing given by the State's linearIn/linearOut functions.
  */
  std::vector<std::vector<int>> m_nu;
};

#include <CD_NamespaceFooter.H>

#include <CD_KMCSolverImplem.H>

#endif
