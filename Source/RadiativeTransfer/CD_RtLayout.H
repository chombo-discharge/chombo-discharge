/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
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
  @brief Class for holding a set of RtSolvers. 
*/
template <class T>
class RtLayout {
public:

  /*!
    @brief Full constructor
  */
  RtLayout(const Vector<RefCountedPtr<RtSpecies> >& a_RtSpecies);

  /*!
    @brief Destructor
  */
  virtual ~RtLayout();

  /*!
    @brief Get iterator
  */
  virtual RtIterator<T> iterator();

  /*!
    @brief Get Realm
  */
  virtual const std::string getRealm() const;

  /*!
    @brief Set Realm
  */
  virtual void setRealm(const std::string a_realm);

  /*!
    @brief Add solver
  */
  virtual void addSolver(RefCountedPtr<T> a_solver);

  /*!
    @brief Parse options
  */
  virtual void parseOptions();

  /*!
    @brief Parse options
  */
  virtual void parseRuntimeOptions();

  /*!
    @brief Allocate internal storage for solvers
  */
  virtual void allocateInternals();

  /*!
    @brief Pre regrid stuff
  */
  virtual void preRegrid(const int a_base, const int a_oldFinestLevel);

  /*!
    @brief Deallocate internal storage for solvers
  */
  virtual void deallocateInternals();

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
    @brief Do a sanity check
  */
  virtual void sanityCheck();

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng setComputationalGeometry
  */
  virtual void setPhase(phase::which_phase a_phase = phase::gas);

  /*!
    @brief Set verbosity.
  */
  virtual void setVerbosity(const int a_verbosity);

  /*!
    @brief Set time
  */
  virtual void setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Convenience function. Set source terms for all species. Mostly used for debugging. 
  */
  virtual void setSource(const EBAMRCellData& a_source);

  /*!
    @brief Convenience function. Set constant source terms for all species. Mostly used for debugging. 
  */
  virtual void setSource(const Real a_source);

  /*!
    @brief Turn on/off stationary
  */
  virtual void setStationary(const bool a_stationary);

  /*!
    @brief Convenience function. All solvers write plot files
  */
  virtual void writePlotFile();

  /*!
    @brief Convenience function. Call advance method for each solver
  */
  virtual void advance(const Real a_dt);

  /*!
    @brief Fill with initial data
    @details Only relevant for transient solvers
  */
  virtual void initialData();

  /*!
    @brief Check if solvers are stationary
  */
  virtual bool isStationary();

  /*!
    @brief Get phase
  */
  virtual phase::which_phase getPhase();

  /*!
    @brief Get solvers
  */
  virtual Vector<RefCountedPtr<T> >& getSolvers();

  /*!
    @brief Get species
  */
  virtual Vector<RefCountedPtr<RtSpecies> >& getSpecies();

  /*!
    @brief Get all source terms
  */
  virtual Vector<EBAMRCellData*> getSources();

  /*!
    @brief Get all states
  */
  virtual Vector<EBAMRCellData*> getPhis();

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
  Vector<RefCountedPtr<RtSpecies> > m_species;

  /*!
    @brief Verbosity
  */
  int m_verbosity;
};

/*!
  @brief Factory class for RtLayout. 
  @details Factory class is very simple; since we don't want to template RtLayout2 we use a factory to instantiate solvers of 
  any RtSolver inherited class, and then return a layout with the casted classes. That's about it. 
*/
template <class T, class S>
class RtFactory {
public:
  RtFactory();
  ~RtFactory();
  RefCountedPtr<RtLayout<T> > newLayout(const Vector<RefCountedPtr<RtSpecies> >& a_species) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_RtLayoutImplem.H>

#endif
