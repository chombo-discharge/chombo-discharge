/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RtLayout.H
  @brief  Declaration of a class that holds a set of RtSolvers
  @author Robert Marskar
*/

#ifndef CD_RtLayout_H
#define CD_RtLayout_H

// Our includes
#include <CD_RtSolver.H>
#include <CD_NamespaceHeader.H>

template <class T>
class RtIterator;

/*!
  @brief Class for holding a set of RtSolvers. T must derive from RtSolver.
*/
template <class T>
class RtLayout
{
public:
  /*!
    @brief Full constructor
    @param[in] a_RtSpecies List of RTE species.
  */
  RtLayout(const Vector<RefCountedPtr<RtSpecies>>& a_RtSpecies);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~RtLayout();

  /*!
    @brief Get iterator. 
    @returns Returns a new iterator for iterating through the RTE solvers.
  */
  virtual RtIterator<T>
  iterator();

  /*!
    @brief Get realm where the solvers are allocated.
  */
  virtual const std::string
  getRealm() const;

  /*!
    @brief Set realm for solvers
    @param[in] a_realm Realm name
  */
  virtual void
  setRealm(const std::string a_realm);

  /*!
    @brief Add a new solver to the solver layout.
    @param[in] a_solver New RTE solver. 
  */
  virtual void
  addSolver(RefCountedPtr<T> a_solver);

  /*!
    @brief Parse options.
  */
  virtual void
  parseOptions();

  /*!
    @brief Parse runtime options.
  */
  virtual void
  parseRuntimeOptions();

  /*!
    @brief Allocate internal storage for solvers
  */
  virtual void
  allocate();

  /*!
    @brief Pre regrid stuff
    @param[in] a_base           Coarsest level which changes.
    @param[in] a_oldFinestLevel Finest AMR level before regrid.
  */
  virtual void
  preRegrid(const int a_base, const int a_oldFinestLevel);

  /*!
    @brief Deallocate internal storage for solvers
  */
  virtual void
  deallocate();

  /*!
    @brief Regrid all solvers.
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level after the regrid. 
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Register operators
  */
  virtual void
  registerOperators();

  /*!
    @brief Set amr
    @param[in] a_amr AmrMesh object. 
  */
  virtual void
  setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set the computational geometry
    @param[in] a_computationalGeometry Computational geometry.
  */
  virtual void
  setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

  /*!
    @brief Do a sanity check
  */
  virtual void
  sanityCheck();

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng setComputationalGeometry
    @param[in] a_phase Phase
  */
  virtual void
  setPhase(phase::which_phase a_phase = phase::gas);

  /*!
    @brief Set verbosity.
    @param[in] a_verbosity Verbosity level for solvers. 
  */
  virtual void
  setVerbosity(const int a_verbosity);

  /*!
    @brief Set time
    @param[in] a_step Time step number
    @param[in] a_time Time (in seconds)
    @param[in] a_dt   Time step increment
    @note This sets m_step=a_step, m_time=a_time, m_dt=a_dt for all solvers. 
  */
  virtual void
  setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Convenience function. Set source terms for all species. Mostly used for debugging. 
    @param[in] a_source Source term. Copied into every solver's source term. 
  */
  virtual void
  setSource(const EBAMRCellData& a_source);

  /*!
    @brief Convenience function. Set constant source terms for all species. Mostly used for debugging. 
    @param[in] a_source Source term set on every solver. 
  */
  virtual void
  setSource(const Real a_source);

  /*!
    @brief Turn on/off stationary
    @param[in] a_stationary If true, every RTE solver will be set to stationary. 
  */
  virtual void
  setStationary(const bool a_stationary);

  /*!
    @brief Convenience function. All solvers write plot files
  */
  virtual void
  writePlotFile();

  /*!
    @brief Convenience function. Call advance method for each solver
    @param[in] a_dt Advance method. 
  */
  virtual void
  advance(const Real a_dt);

  /*!
    @brief Fill with initial data
    @details Only relevant for transient solvers
  */
  virtual void
  initialData();

  /*!
    @brief Check if solvers are stationary
  */
  virtual bool
  isStationary();

  /*!
    @brief Get phase
  */
  virtual phase::which_phase
  getPhase();

  /*!
    @brief Get solvers
  */
  virtual Vector<RefCountedPtr<T>>&
  getSolvers();

  /*!
    @brief Get species
  */
  virtual Vector<RefCountedPtr<RtSpecies>>&
  getSpecies();

  /*!
    @brief Get all source terms
  */
  virtual Vector<EBAMRCellData*>
  getSources();

  /*!
    @brief Get all states
  */
  virtual Vector<EBAMRCellData*>
  getPhis();

protected:
  /*!
    @brief Realm
  */
  std::string m_realm;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Solvers
  */
  Vector<RefCountedPtr<T>> m_solvers;

  /*!
    @brief Species
  */
  Vector<RefCountedPtr<RtSpecies>> m_species;

  /*!
    @brief Verbosity
  */
  int m_verbosity;
};

/*!
  @brief Factory class for RtLayout. 
  @details Factory class is very simple; since we don't want to template RtLayout we use a factory to instantiate solvers of 
  any RtSolver inherited class, and then return a layout with the casted classes. That's about it. 
*/
template <class T, class S>
class RtFactory
{
public:
  /*!
    @brief Constructor (does nothing)
  */
  RtFactory();

  /*!
    @brief Destructor (does nothing)
  */
  ~RtFactory();

  /*!
    @brief Get a new Layout. This will cast S to a specific class (T) 
    @param[in] a_species RTE species. 
  */
  RefCountedPtr<RtLayout<T>>
  newLayout(const Vector<RefCountedPtr<RtSpecies>>& a_species) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_RtLayoutImplem.H>

#endif
