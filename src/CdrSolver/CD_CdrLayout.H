/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_CdrLayout.H
  @brief  Declaration of a class that holds a set of CdrSolvers
  @author Robert Marskar
*/

#ifndef CD_CdrLayout_H
#define CD_CdrLayout_H


// Our includes
#include <CD_CdrSolver.H>
#include <CD_NamespaceHeader.H>

// Needs forward declaration.
template <class T>
class cdr_iterator;

/*!
  @brief Class for holding a set of CdrSolvers. 
*/
template <class T>
class CdrLayout {
public:

  CdrLayout();

  /*!
    @brief Full constructor
  */
  CdrLayout(const Vector<RefCountedPtr<CdrSpecies> >& a_plaskin);

  /*!
    @brief Destructor
  */
  virtual ~CdrLayout();

  /*!
    @brief Get phase
  */
  virtual phase::which_phase getPhase() const;

  /*!
    @brief Get the iterator
  */
  virtual cdr_iterator<T> iterator();

  /*!
    @brief Get the Realm
  */
  virtual const std::string getRealm() const;

  /*!
    @brief Set
  */
  virtual void setRealm(const std::string a_realm);

  /*!
    @brief Parse class options
  */
  virtual void parseOptions();

  /*!
    @brief Parse runtime options
  */
  virtual void parseRuntimeOptions();

  /*!
    @brief Add solver
  */
  virtual void addSolver(RefCountedPtr<T> a_solver);

  /*!
    @brief Allocate internal storage for solvers
  */
  virtual void allocateInternals();

  /*!
    @brief Deallocate solver internals
  */
  virtual void deallocateInternals();

  /*!
    @brief Cache states before regrid step
  */
  virtual void preRegrid(const int a_lbase, const int a_finestLevel);

  /*!
    @brief Initial data
  */
  virtual void initialData();

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
  virtual void setComputationalGeometry(const RefCountedPtr<computational_geometry>& a_computationalGeometry);

  /*!
    @brief Do a sanity check
  */
  virtual void sanityCheck();

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng setComputationalGeometry
  */
  virtual void setPhase(phase::which_phase a_phase = phase::gas);

  /*!
    @brief Set phase
  */
  virtual void set_bc(CdrBc a_bc);

  /*!
    @brief Set verbosity.
  */
  virtual void setVerbosity(const int a_verbosity);

  /*!
    @brief Set time
  */
  virtual void setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Convenience function. Set velocity for all species. Mostly used for debugging.
  */
  virtual void setVelocity(const EBAMRCellData& a_velo);

  /*!
    @brief Convenience function. Set constant velocity for all species. Mostly used for debugging. 
  */
  virtual void setVelocity(const RealVect a_velo);

  /*!
    @brief Convenience function. Set diffusion coefficients for all species. Mostly used for debugging. 
  */
  virtual void setDiffusionCoefficient(const EBAMRFluxData& a_diffusionCoefficient, const EBAMRIVData& a_ebDiffusionCoefficient);

  /*!
    @brief Convenience function. Set constant diffusion coefficients for all species. Mostly used for debugging. 
  */
  virtual void setDiffusionCoefficient(const Real a_diffusionCoefficient);

  /*!
    @brief Convenience function. Set source terms for all species. Mostly used for debugging. 
  */
  virtual void setSource(const EBAMRCellData& a_source);

  /*!
    @brief Convenience function. Set constant source terms for all species. Mostly used for debugging. 
  */
  virtual void setSource(const Real a_source);

  /*!
    @brief Convenience function. Set the EB flux for all species. Mostly used for debugging. 
  */
  virtual void setEbFlux(const EBAMRIVData& a_ebFlux);

  /*!
    @brief Convenience function.Set the EB flux for all species. Mostly used for debugging. 
  */
  virtual void setEbFlux(const Real a_ebFlux);

  /*!
    @brief Convenience function. All solvers write plot files
  */
  virtual void writePlotFile();

  /*!
    @brief Get largest possible time step for advection part
  */
  virtual Real computeAdvectionDt();

  /*!
    @brief Get largest possible time step for diffusion part
  */
  virtual Real computeDiffusionDt();

  /*!
    @brief Get largest possible time step for diffusion part
  */
  virtual Real computeAdvectionDiffusionDt();

  /*!
    @brief Get largest possible time step for source term
  */
  virtual Real computeSourceDt(const Real a_tolerance, const bool a_elec_only);

  /*!
    @brief Compute total charge in computational region
  */
  virtual Real computeTotalCharge();

  /*!
    @brief Compute mass for every species
  */
  virtual Vector<Real> computeMass();

  /*!
    @brief Compute charge for each species
  */
  virtual Vector<Real> computeCharge();

  /*!
    @brief Get all names
  */
  virtual Vector<std::string> getNames();

  /*!
    @brief Get solvers
  */
  virtual Vector<RefCountedPtr<T> >& getSolvers();

  /*!
    @brief Get species
  */
  virtual Vector<RefCountedPtr<CdrSpecies> >& getSpecies();

  /*!
    @brief Get all source terms
  */
  virtual Vector<EBAMRCellData*> getSources();

  /*!
    @brief Get all states
  */
  virtual Vector<EBAMRCellData*> getPhis();

  /*!
    @brief Get all velocities
  */
  virtual Vector<EBAMRCellData*> getVelocities();

  /*!
    @brief Get all diffusion coefficients
  */
  virtual Vector<EBAMRFluxData*> getFaceCenteredDiffusionCoefficient();

  /*!
    @brief Get all diffusion coefficients
  */
  virtual Vector<EBAMRIVData*> getEbCenteredDiffusionCoefficient();

  /*!
    @brief Get the EB flux
  */
  virtual Vector<EBAMRIVData*> getEbFlux();

  /*!
    @brief Get domain fluxes
  */
  virtual Vector<EBAMRIFData*> getDomainFlux();

  /*!
    @brief Get mass redistribution
  */
  virtual bool getMassRedist();

protected:

  /*!
    @brief Realm
  */
  std::string m_realm;

  /*!
    @brief Solvers
  */
  Vector<RefCountedPtr<CdrSolver> > m_solvers;

  /*!
    @brief Species
  */
  Vector<RefCountedPtr<CdrSpecies> > m_species;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Mass redistribution for hyperbolic part
  */
  bool m_mass_redist;

  /*!
    @brief phase
  */
  phase::which_phase m_phase;

  /*!
    @brief
  */
  CdrBc m_bc_type;
};

/*!
  @brief Factory class for CdrLayout. 
  @details Factory class is very simple; since we don't want to template CdrLayout2 we use a factory to instantiate solvers of 
  any CdrSolver inherited class, and then return a layout with the casted classes. That's about it. 
*/
template <class T, class S>
class CdrFactory {
public:
  CdrFactory();
  ~CdrFactory();
  RefCountedPtr<CdrLayout<T> > newLayout(const Vector<RefCountedPtr<CdrSpecies> > a_species) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_CdrLayoutImplem.H>

#endif
