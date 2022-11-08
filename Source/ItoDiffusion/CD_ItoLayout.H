/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoLayout.H
  @brief  Declaration of a class that holds a set of ItoSolvers
  @author Robert Marskar
*/

#ifndef CD_ItoLayout_H
#define CD_ItoLayout_H

// Our includes
#include <CD_ItoSolver.H>
#include <CD_NamespaceHeader.H>

template <class T>
class ItoIterator;

/*!
  @brief Class for holding a set of ItoSolvers.
  @details This is primarily a utility class just cuts down on typing when there are multiple species being advected/diffused. This class just stores
  a set of such solvers and contains an "iterator" for going through them. Fundamentally, there is no difference between this class and a vector of 
  ItoSolvers. This class can be constructed through a factory method. 
*/
template <class T>
class ItoLayout
{
public:
  /*!
    @brief Disallowed weak construction.
  */
  ItoLayout() = delete;

  /*!
    @brief Full constructor. Calls the define function (but does not generate solvers). 
    @param[in] a_species List of ItoSpecies.
  */
  ItoLayout(const Vector<RefCountedPtr<ItoSpecies>>& a_species);

  /*!
    @brief Destructor (does nothing).
  */
  ~ItoLayout();

  /*!
    @brief Return phase where the solvers are defined. 
  */
  virtual phase::which_phase
  getPhase() const;

  /*!
    @brief Return a fresh iterator. The iterator is a simple random access object that can iterate through the solvers. 
  */
  virtual ItoIterator<T>
  iterator();

  /*!
    @brief Get the realm name where the solvers are defined. 
  */
  virtual const std::string
  getRealm() const;

  /*!
    @brief Set the realm where the solvers are defined.
    @param[in] a_realm Realm name. 
  */
  virtual void
  setRealm(const std::string a_realm);

  /*!
    @brief Define function. 
    @param[in] a_species List of species (for which we generate solvers). 
  */
  virtual void
  define(const Vector<RefCountedPtr<ItoSpecies>>& a_species);

  /*!
    @brief Utility function -- parse options for all solvers
  */
  virtual void
  parseOptions();

  /*!
    @brief Utility function -- parse runtime options for all solvers
  */
  virtual void
  parseRuntimeOptions();

  /*!
    @brief Allocate internals for all solvers
    @details This will call allocateInternals for each ItoSolver. 
  */
  virtual void
  allocateInternals();

  /*!
    @brief Add solver to the list of solvers.
    @param[in] a_solver New solver. 
  */
  virtual void
  addSolver(RefCountedPtr<T>& a_solver);

  /*!
    @brief Utility function which caches states before regrid step.
    @details Calls the preRegrid function for each solver. 
    @param[in] a_lbase          Coarsest level that changed during regrid. 
    @param[in] a_oldFinestLevel Finest grid level before the regrid operation. 
  */
  virtual void
  preRegrid(const int a_lbase, const int a_finestLevel);

  /*!
    @brief Fill all solvers with initial data. 
    @details Calls initialData() for each solver. 
  */
  virtual void
  initialData();

  /*!
    @brief All solvers deposit their particles
    @details Calls depositParticles for each solver. 
  */
  virtual void
  depositParticles();

  /*!
    @brief Deposit particles for each solver. Input argument determines which container is deposited. 
    @details Calls depositParticles(ItoSolver::WhichContainer) for each solver. 
    @param[in] a_whichContainer Container type (e.g., ItoSolver::WhichContainer::bulk)
  */
  virtual void
  depositParticles(const ItoSolver::WhichContainer a_whichContainer);

  /*!
    @brief Remap function -- calls remap(ItoSolver::WhichContainer::bulk) for each solver. 
  */
  virtual void
  remap();

  /*!
    @brief Remap function -- all solves remap the input container.
    @param[in] a_whichContainer Container type (e.g., ItoSolver::WhichContainer::bulk)
  */
  virtual void
  remap(const ItoSolver::WhichContainer a_whichContainer);

  /*!
    @brief Interpolate velocities
    @details Calls interpolateVelocities() for each solver. 
  */
  virtual void
  interpolateVelocities();

  /*!
    @brief Regrid method -- calls the regrid method for each solver. 
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level after the regrid. 
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Register operators. Calls registerOperators for each solver. 
  */
  virtual void
  registerOperators();

  /*!
    @brief Set AmrMesh object for each solver. 
    @param[in] a_amrMesh AmrMesh object
  */
  virtual void
  setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set the computational geometry. This sets the computational geometry for each solver.
    @param[in] a_computationalGeometry Computational geometry. 
  */
  virtual void
  setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

  /*!
    @brief Set phase for each solver. 
    @param[in] a_phase Phase
  */
  virtual void
  setPhase(phase::which_phase a_phase);

  /*!
    @brief Set verbosity for each solver.
    @param[in] a_verbosity Verbosity level. 
  */
  virtual void
  setVerbosity(const int a_verbosity);

  /*!
    @brief Set time for each solver. 
    @param[in] a_step Time step number
    @param[in] a_time Time (in seconds)
    @param[in] a_dt   Time step increment
    @note This sets m_step=a_step, m_time=a_time, m_dt=a_dt
  */
  virtual void
  setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Sort ItoSolver particles by cell -- this is done for all solvers. 
    @param[in] a_whichContainer Which ItoSolver container to sort by cell. 
  */
  virtual void
  organizeParticlesByCell(const ItoSolver::WhichContainer a_whichContainer);

  /*!
    @brief Sort ItoSolver particles by patch -- this is done for all solvers. 
    @param[in] a_whichContainer Which ItoSolver container to sort by cell. 
  */
  virtual void
  organizeParticlesByPatch(const ItoSolver::WhichContainer a_whichContainer);

  /*!
    @brief Rearrange the input container particle into new superparticles -- this is done for all solvers. 
    @param[in] a_whichContainer Which ItoSolver container to sort by cell. 
    @param[in] a_ppc            Target particle-per-cell (ppc)
  */
  virtual void
  makeSuperparticles(const ItoSolver::WhichContainer a_whichContainer, const int a_ppc);

  /*!
    @brief Compute smallest possible time step. 
    @details This calls computeDt() for each ItoSolver and minimizes the result over the solvers (i.e. it returns the smallest time step)
  */
  virtual Real
  computeDt();

  /*!
    @brief Compute the classical advection time step for all solvers. This returns dt = dx/max(v) where max(v) takes the largest component. 
    @details This calls computeAdvectionDt() for each ItoSolver and minimizes the result over the solvers (i.e. it returns the smallest time step)
    @note See ItoSolver::computeAdvectionDt() to see how the time step is computed. 
  */
  virtual Real
  computeAdvectiveDt();

  /*!
    @brief Compute the largest dt which restricts all particles to move less than a_maxCellsToMove.
    @details This calls computeHopDt for all solvers -- that function uses the maximum permitted diffusion hop to restrict the time step (the user will have
    restricted the normal distribution to some value). So, the result of this value will depend strongly on the user input for the ItoSolvers. 
    @param[in] a_maxCellsToMove Maximum number of cells to move. 
  */
  virtual Real
  computeHopDt(const Real a_maxCellsToMove);

  /*!
    @brief Compute the classical diffusive time step dt = dx*dx/(2*D) where D is the diffusion coefficient (not dimensional dependence here)
    @details This calls computeDiffusiveDt for all solvers and minimizes the result over all solvers. 
  */
  virtual Real
  computeDiffusiveDt();

  /*!
    @brief Get total number of particles. 
    @details This returns the total number of particles in a specified container. The result is accumulated over all solvers.
    @param[in] a_whichContainer Container type in ItoSolver
    @param[in] a_localOnly      If true, only count particles locally. If false, compute globally. 
  */
  virtual size_t
  getNumParticles(const ItoSolver::WhichContainer a_ptype, const bool a_localOnly);

  /*!
    @brief Get solvers all solvers
    @return Returns m_solvers
  */
  virtual Vector<RefCountedPtr<T>>&
  getSolvers();

  /*!
    @brief Get species
    @return Returns m_species
  */
  virtual Vector<RefCountedPtr<ItoSpecies>>&
  getSpecies();

  /*!
    @brief Get all velocity functions. 
    @details This accumulates a pointer to the velocity mesh data for each solver. 
  */
  virtual Vector<EBAMRCellData*>
  getVelocityFunctions();

  /*!
    @brief Get all densities
    @details This accumulates a pointer to the density mesh data for each solver. 
  */
  virtual Vector<EBAMRCellData*>
  getDensities();

  /*!
    @brief Get all scratch storage
    @details This accumulates a pointer to scratch mesh data for each solver. 
  */
  virtual Vector<EBAMRCellData*>
  getScratch();

  /*!
    @brief Get all diffusion coefficients
    @details This accumulates a pointer to diffusion coefficient mesh data for each solver. 
  */
  virtual Vector<EBAMRCellData*>
  getDiffusionFunctions();

  /*!
    @brief Get all mobility mesh functions
    @details This accumulates a pointer to mobility coefficient mesh data for each solver. 
  */
  virtual Vector<EBAMRCellData*>
  getMobilityFunctions();

  /*!
    @brief Get particle containers of a particular container type
    @param[in] a_whichContainer Which ItoSolver container
  */
  virtual Vector<ParticleContainer<ItoParticle>*>
  getParticles(const ItoSolver::WhichContainer a_whichContainer);

protected:
  /*!
    @brief Realm where all the solvers are registered
  */
  std::string m_realm;

  /*!
    @brief ItoSolver solvers
  */
  Vector<RefCountedPtr<T>> m_solvers;

  /*!
    @brief ItoSolver species
  */
  Vector<RefCountedPtr<ItoSpecies>> m_species;

  /*!
    @brief Verbosity level. 
  */
  int m_verbosity;

  /*!
    @brief Phase where solvers live. 
  */
  phase::which_phase m_phase;

  /*!
    @brief Layout is defined or not. 
  */
  bool m_isDefined;
};

/*!
  @brief Factory class for making ItoLayout. 
  @details This is used for generating a layout -- this class has no fundamental usefulness, it is basically a glorified casting machine. The template parameters T control
  the ItoLayout template parameter while the template parameter S controls the derived class. E.g. if one has a derived class S the factory can automatically create either
  an ItoLayout<S>, or it can automatically cast to an ItoLayout<T>. 
*/
template <class T, class S>
class ItoFactory
{
public:
  /*!
    @brief Constructor.
  */
  ItoFactory();

  /*!
    @brief Destructor (does nothing).
  */
  ~ItoFactory();

  /*!
    @brief Factory method which creates a new layout from a set of species. This can do automated casting between the template parameters (in case there are classes
    that derive from ItoSolver).
    @param[in] a_species List of species which will be solved for.
  */
  RefCountedPtr<ItoLayout<T>>
  newLayout(const Vector<RefCountedPtr<ItoSpecies>>& a_species) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_ItoLayoutImplem.H>

#endif
