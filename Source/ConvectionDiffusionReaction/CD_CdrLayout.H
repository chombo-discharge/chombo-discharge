/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrLayout.H
  @brief  Declaration of a class that holds a set of CdrSolvers (to cut down on typing). 
  @author Robert Marskar
*/

#ifndef CD_CdrLayout_H
#define CD_CdrLayout_H

// Our includes
#include <CD_CdrSolver.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Forward declare iterator
  @see CD_CdrIterator.H
*/
template <class T>
class CdrIterator;

/*!
  @brief Class for holding a set of CdrSolvers. 
  @details This is primarily a utility class that just cuts down on typing when there are multiple species are being advected/diffused. This class
  simply stores a set of CdrSolvers, and it contains an "iterator" for iterating through them. This is simply meant to cut down on typing. Note that a 
  TimeStepper can keep an instance this class rather than e.g. a Vector<RefCountedPtr<CdrSolver> >, but fundamentally there is no difference between the two. 
  The CdrLayout can be constructed through a factory method, i.e. the user can do a quick-generation of a set of CdrSolvers by using the factory method for
  instantiating an arbitrary number of CdrSolvers. 
  @note A separate use case for this is class is that it allows segregation between groups of CdrSolvers that live on different realms and phases. 
*/
template <class T>
class CdrLayout
{
public:
  /*!
    @brief Disallowed weak constructor. Use the full constructor. 
  */
  CdrLayout() = delete;

  /*!
    @brief Full constructor, setting the number of solvers to be instantiated. 
    @param[in] a_species List of species that defines this solver. 
  */
  CdrLayout(const Vector<RefCountedPtr<CdrSpecies>>& a_species);

  /*!
    @brief Destructor
  */
  virtual ~CdrLayout();

  /*!
    @brief Get phase
    @return Return the phase on which all CdrSolvers are allocated. 
  */
  virtual phase::which_phase
  getPhase() const;

  /*!
    @brief Get an iterator which can iterate through the various Cdr solvers. 
  */
  virtual CdrIterator<T>
  iterator();

  /*!
    @brief Get the realm where the Cdr solvers live. 
    @return Return realm name. 
  */
  virtual std::string
  getRealm() const;

  /*!
    @brief Set realm for all CdrSolvers
    @param[in] a_realm Realm name
  */
  virtual void
  setRealm(const std::string a_realm);

  /*!
    @brief Parse class options
  */
  virtual void
  parseOptions();

  /*!
    @brief Parse runtime options
  */
  virtual void
  parseRuntimeOptions();

  /*!
    @brief Add solver. This appends a_solver to m_solvers
    @param[in] a_solver Solver to add to the layout. 
  */
  virtual void
  addSolver(RefCountedPtr<T>& a_solver);

  /*!
    @brief Allocate internal storage for solvers
  */
  virtual void
  allocate();

  /*!
    @brief Deallocate solver internals
  */
  virtual void
  deallocate();

  /*!
    @brief Call preRegrid(...) method for all the CdrSolvers
    @param[in] a_lbase          Coarsest level which changes during regrid
    @param[in] a_oldFinestLevel Finest level before the regrid (i.e., the current finest level)
  */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel);

  /*!
    @brief Initial data. Fill all solvers with initial data.
  */
  virtual void
  initialData();

  /*!
    @brief Regrid method. Calls regrid(...) for every CdrSolver
    @param[in] a_lmin           Coarsest level that changed during the regrid. 
    @param[in] a_oldFinestLevel Finest level before the regrid
    @param[in] a_newFinestLevel Finest level after the regrid
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Register operators
  */
  virtual void
  registerOperators();

  /*!
    @brief Set AmrMesh for all CdrSolvers
    @param[in] a_amrMesh AmrMesh
  */
  virtual void
  setAmr(const RefCountedPtr<AmrMesh>& a_amrMesh);

  /*!
    @brief Set the computational geometry for all CdrSolvers
    @param[in] a_computationalGeometry Computational geometry
  */
  virtual void
  setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

  /*!
    @brief Set phase
    @param[i] a_phase Phase
  */
  virtual void
  setPhase(const phase::which_phase a_phase = phase::gas);

  /*!
    @brief Set verbosity for all CdrSolvers
    @param[in] a_verbosity Verbosity
  */
  virtual void
  setVerbosity(const int a_verbosity);

  /*!
    @brief Set time
    @param[in] a_step Time step
    @param[in] a_time Time
    @param[in] a_dt   Time increment
  */
  virtual void
  setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Convenience function. All solvers write plot files (for debugging purposes)
  */
  virtual void
  writePlotFile();

  /*!
    @brief Get CFL time for advection
    @return Returns the smallest explicit advective time step (minimized over solvers)
  */
  virtual Real
  computeAdvectionDt();

  /*!
    @brief Get time step for explicit diffusion
    @return Returns the smallest explicit diffusion time step (minimized over solvers)
  */
  virtual Real
  computeDiffusionDt();

  /*!
    @brief Get the time step for explicit advection-diffusion
    @return Returns the smallest explicit advection-diffusion time step (minimized over solvers)
  */
  virtual Real
  computeAdvectionDiffusionDt();

  /*!
    @brief Get solvers
    @return Returns all CdrSolvers in this layout. 
  */
  virtual Vector<RefCountedPtr<T>>&
  getSolvers();

  /*!
    @brief Get species
  */
  virtual Vector<RefCountedPtr<CdrSpecies>>&
  getSpecies();

  /*!
    @brief Get all source terms
  */
  virtual Vector<EBAMRCellData*>
  getSources();

  /*!
    @brief Get all cell-centered states. 
  */
  virtual Vector<EBAMRCellData*>
  getPhis();

  /*!
    @brief Get all cell-centered velocities. 
  */
  virtual Vector<EBAMRCellData*>
  getVelocities();

  /*!
    @brief Get all cell-centered diffusion coefficients. 
  */
  virtual Vector<EBAMRCellData*>
  getCellCenteredDiffusionCoefficients();

  /*!
    @brief Get all face-centered diffusion coefficients. 
  */
  virtual Vector<EBAMRFluxData*>
  getFaceCenteredDiffusionCoefficient();

  /*!
    @brief Get all EB-centered diffusion coefficients
  */
  virtual Vector<EBAMRIVData*>
  getEbCenteredDiffusionCoefficient();

  /*!
    @brief Get all EB flux data holders. 
  */
  virtual Vector<EBAMRIVData*>
  getEbFlux();

  /*!
    @brief Get domain flux data holders. 
  */
  virtual Vector<EBAMRIFData*>
  getDomainFlux();

protected:
  /*!
    @brief Realm where solvers live. 
  */
  std::string m_realm;

  /*!
    @brief Phase where solvers live. 
  */
  phase::which_phase m_phase;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Solver instantiations. 
  */
  Vector<RefCountedPtr<CdrSolver>> m_solvers;

  /*!
    @brief Cdr species
  */
  Vector<RefCountedPtr<CdrSpecies>> m_species;
};

/*!
  @brief Factory class for CdrLayout. T is (usually) CdrSolver and S is the implementation class (e.g., CdrCTU). 
  @details Factory class is very simple; since we don't want to template CdrLayout further we use a factory to instantiate solvers of 
  any CdrSolver inherited class, and then return a layout with the casted classes. That's about it. This is useful because it allows
  the factory to act as a "casting-point" for switching between various solver implementation. 
*/
template <class T, class S>
class CdrFactory
{
public:
  /*!
    @brief Default constructor (does nothing)
  */
  CdrFactory();

  /*!
    @brief Default constructor (does nothing)
  */
  ~CdrFactory();

  /*!
    @brief Factory method, create a new CdrLayout
    @details This creates a collection of CdrSolvers. Run-time errors will occur if the template S is not a derived class of T. 
  */
  RefCountedPtr<CdrLayout<T>>
  newLayout(const Vector<RefCountedPtr<CdrSpecies>>& a_species) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_CdrLayoutImplem.H>

#endif
