/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
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
class ito_iterator;

/*!
  @brief Class for holding a set of ItoSolvers
*/
template <class T>
class ItoLayout{
public:

  /*!
    @brief Default constructor. Must subsequently call define. 
  */
  ItoLayout();

  /*!
    @brief Full constructor. Calls define. 
  */
  ItoLayout(const Vector<RefCountedPtr<ItoSpecies> >& a_species);

  /*!
    @brief Destructor. 
  */
  ~ItoLayout();

  /*!
    @brief Return phase
  */
  virtual phase::which_phase getPhase() const;

  /*!
    @brief Return iterator. The iterator is a simple random access object
  */
  virtual ito_iterator<T> iterator();

  /*!
    @brief Get Realm
  */
  virtual const std::string getRealm() const;

  /*!
    @brief Set Realm
  */
  virtual void setRealm(const std::string a_realm);

  /*!
    @brief Define function
  */
  virtual void define(const Vector<RefCountedPtr<ItoSpecies> >& a_species);

  /*!
    @brief Parse options for all solvers
  */
  virtual void parseOptions();

  /*!
    @brief Parse options for all solvers
  */
  virtual void parseRuntimeOptions();

  /*!
    @brief Set PVR buffer for all solvers
  */
  virtual void setPVRBuffer(const int a_buffer);

  /*!
    @brief Set the halo buffer
  */
  virtual void setHalobuffer(const int a_buffer);

  /*!
    @brief Allocate internals for all solvers
  */
  virtual void allocateInternals();

  /*!
    @brief Add solver
  */
  virtual void addSolver(RefCountedPtr<T> a_solver);

  /*!
    @brief Cache states before regrid step
  */
  virtual void preRegrid(const int a_lbase, const int a_finestLevel);

  /*!
    @brief Initial data
  */
  virtual void initialData();

  /*!
    @brief All solvers deposit their particles
  */
  virtual void depositParticles();

  /*!
    @brief Deposit particles
  */
  virtual void depositParticles(const ItoSolver::WhichContainer a_ptype);

  /*!
    @brief All solvers remap. 
  */
  virtual void remap();

  /*!
    @brief All solvers remap. 
  */
  virtual void remap(const ItoSolver::WhichContainer a_ptype);

  /*!
    @brief Interpolate velocities
  */
  virtual void interpolateVelocities();

  /*!
    @brief Regrid method. 
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Register operators
  */
  virtual void registerOperators();

  /*!
    @brief Set amr
  */
  virtual void setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set the computational geometry
  */
  virtual void setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng setComputationalGeometry
  */
  virtual void setPhase(phase::which_phase a_phase);

  /*!
    @brief Set verbosity.
  */
  virtual void setVerbosity(const int a_verbosity);

  /*!
    @brief Set time
  */
  virtual void setTime(const int a_step, const Real a_time, const Real a_dt);
  
  /*!
    @brief Sort the m_particles container by cell
  */
  virtual void sortParticlesByCell(const ItoSolver::WhichContainer a_ptype);

  /*!
    @brief Sort the m_particles container by cell
  */
  virtual void sortParticlesByPatch(const ItoSolver::WhichContainer a_ptype);

  /*!
    @brief Make superparticles
  */
  virtual void makeSuperparticles(const ItoSolver::WhichContainer a_ptype, const int a_ppc);

  /*!
    @brief Compute dt
  */
  virtual Real computeDt();

  /*!
    @brief Compute the classical advection time step
  */
  virtual Real computeAdvectiveDt();

  /*!
    @brief Compute dt
  */
  virtual Real computeMinDt(const Real a_maxCellsToMove);

  /*!
    @brief Get largest possible time step for advection part
  */
  virtual Real computeMinDriftDt(const Real a_maxCellsToMove);

  /*!
    @brief Get largest possible time step for diffusion part
  */
  virtual Real computeMinDiffusionDt(const Real a_maxCellsToMove);

  /*!
    @brief Compute the classical diffusive time step
  */
  virtual Real computeDiffusiveDt();

  /*!
    @brief Get number of particles
  */
  virtual size_t getNumParticles(const ItoSolver::WhichContainer a_ptype, const bool a_local);

  /*!
    @brief Get solvers
  */
  virtual Vector<RefCountedPtr<T> >& getSolvers();

  /*!
    @brief Get species
  */
  virtual Vector<RefCountedPtr<ItoSpecies> >& getSpecies();

  /*!
    @brief Get all velocity functions
  */
  virtual Vector<EBAMRCellData*> getVelocityFunctions();

  /*!
    @brief Get all densities
  */
  virtual Vector<EBAMRCellData*> getDensities();

  /*!
    @brief Get all scratch storage
  */
  virtual Vector<EBAMRCellData*> getScratch();

  /*!
    @brief Get all diffusion coefficients
  */
  virtual Vector<EBAMRCellData*> getDiffusionFunctions();

  /*!
    @brief Get all mobility mesh functions
  */
  virtual Vector<EBAMRCellData*> getMobilityFunctions();

  /*!
    @brief Get particle containers
  */
  virtual Vector<ParticleContainer<ItoParticle>* > getParticles(const ItoSolver::WhichContainer a_ptype);

protected:

  /*!
    @brief Realm
  */
  std::string m_realm;
  
  /*!
    @brief Solvers
  */
  Vector<RefCountedPtr<T> > m_solvers;

  /*!
    @brief Species
  */
  Vector<RefCountedPtr<ItoSpecies> > m_species;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Defined or not
  */
  bool m_isDefined;
  
};

/*!
  @brief Factory class for ItoLayout
*/
template <class T, class S>
class itoFactory {
public:
  itoFactory();
  ~itoFactory();
  RefCountedPtr<ItoLayout<T> > newLayout(const Vector<RefCountedPtr<ItoSpecies> > a_species) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_ItoLayoutImplem.H>

#endif
