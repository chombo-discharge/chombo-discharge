/*!
  @file   rte_layout.H
  @brief  Declaration of a class that holds a set of rte_solvers
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _RTE_LAYOUT_
#define _RTE_LAYOUT_

#include "rte_solver.H"

#include "CD_NamespaceHeader.H"
  
template <class T>
class rte_iterator;

/*!
  @brief Class for holding a set of rte_solvers. 
*/
template <class T>
class rte_layout {
public:

  /*!
    @brief Full constructor
  */
  rte_layout(const Vector<RefCountedPtr<rte_species> >& a_rte_species);

  /*!
    @brief Destructor
  */
  virtual ~rte_layout();

  /*!
    @brief Get iterator
  */
  virtual rte_iterator<T> iterator();

  /*!
    @brief Get Realm
  */
  virtual const std::string get_Realm() const;

  /*!
    @brief Set Realm
  */
  virtual void set_Realm(const std::string a_Realm);

  /*!
    @brief Add solver
  */
  virtual void add_solver(RefCountedPtr<T> a_solver);

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
  virtual void pre_regrid(const int a_base, const int a_oldFinestLevel);

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
  virtual void setComputationalGeometry(const RefCountedPtr<computational_geometry>& a_computationalGeometry);

  /*!
    @brief Do a sanity check
  */
  virtual void sanityCheck();

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng setComputationalGeometry
  */
  virtual void set_phase(phase::which_phase a_phase = phase::gas);

  /*!
    @brief Set verbosity.
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Set time
  */
  virtual void set_time(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Convenience function. Set source terms for all species. Mostly used for debugging. 
  */
  virtual void set_source(const EBAMRCellData& a_source);

  /*!
    @brief Convenience function. Set constant source terms for all species. Mostly used for debugging. 
  */
  virtual void set_source(const Real a_source);

  /*!
    @brief Turn on/off stationary
  */
  virtual void set_stationary(const bool a_stationary);

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
  virtual void initial_data();

  /*!
    @brief Check if solvers are stationary
  */
  virtual bool is_stationary();

  /*!
    @brief Get phase
  */
  virtual phase::which_phase get_phase();

  /*!
    @brief Get solvers
  */
  virtual Vector<RefCountedPtr<T> >& get_solvers();

  /*!
    @brief Get species
  */
  virtual Vector<RefCountedPtr<rte_species> >& get_species();

  /*!
    @brief Get all source terms
  */
  virtual Vector<EBAMRCellData*> get_sources();

  /*!
    @brief Get all states
  */
  virtual Vector<EBAMRCellData*> get_states();

protected:

  /*!
    @brief Realm
  */
  std::string m_Realm; 

  /*!
    @brief Solvers
  */
  Vector<RefCountedPtr<T> > m_solvers;

  /*!
    @brief Species
  */
  Vector<RefCountedPtr<rte_species> > m_species;

  /*!
    @brief Verbosity
  */
  int m_verbosity;
};

/*!
  @brief Factory class for rte_layout. 
  @details Factory class is very simple; since we don't want to template rte_layout2 we use a factory to instantiate solvers of 
  any rte_solver inherited class, and then return a layout with the casted classes. That's about it. 
*/
template <class T, class S>
class rte_factory {
public:
  rte_factory();
  ~rte_factory();
  RefCountedPtr<rte_layout<T> > new_layout(const Vector<RefCountedPtr<rte_species> >& a_species) const;
};
#include "CD_NamespaceFooter.H"

#include "rte_layoutI.H"

#endif
