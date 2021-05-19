/*!
  @file   ito_layout.H
  @brief  Declaration of a class that holds a set of ito_solvers
  @author Robert Marskar
  @date   May 2020
*/

#ifndef _ITO_LAYOUT_
#define _ITO_LAYOUT_

#include "ito_solver.H"

#include "CD_NamespaceHeader.H"
  
template <class T>
class ito_iterator;

/*!
  @brief Class for holding a set of ito_solvers
*/
template <class T>
class ito_layout{
public:

  /*!
    @brief Default constructor. Must subsequently call define. 
  */
  ito_layout();

  /*!
    @brief Full constructor. Calls define. 
  */
  ito_layout(const Vector<RefCountedPtr<ito_species> >& a_species);

  /*!
    @brief Destructor. 
  */
  ~ito_layout();

  /*!
    @brief Return phase
  */
  virtual phase::which_phase get_phase() const;

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
  virtual void define(const Vector<RefCountedPtr<ito_species> >& a_species);

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
  virtual void add_solver(RefCountedPtr<T> a_solver);

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
  virtual void deposit_particles();

  /*!
    @brief Deposit particles
  */
  virtual void deposit_particles(const ito_solver::which_container a_ptype);

  /*!
    @brief All solvers remap. 
  */
  virtual void remap();

  /*!
    @brief All solvers remap. 
  */
  virtual void remap(const ito_solver::which_container a_ptype);

  /*!
    @brief Interpolate velocities
  */
  virtual void interpolate_velocities();

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
  virtual void sortParticlesByCell(const ito_solver::which_container a_ptype);

  /*!
    @brief Sort the m_particles container by cell
  */
  virtual void sortParticlesByPatch(const ito_solver::which_container a_ptype);

  /*!
    @brief Make superparticles
  */
  virtual void make_superparticles(const ito_solver::which_container a_ptype, const int a_ppc);

  /*!
    @brief Compute dt
  */
  virtual Real computeDt();

  /*!
    @brief Compute the classical advection time step
  */
  virtual Real compute_advective_dt();

  /*!
    @brief Compute dt
  */
  virtual Real compute_min_dt(const Real a_maxCellsToMove);

  /*!
    @brief Get largest possible time step for advection part
  */
  virtual Real compute_min_drift_dt(const Real a_maxCellsToMove);

  /*!
    @brief Get largest possible time step for diffusion part
  */
  virtual Real compute_min_diffusion_dt(const Real a_maxCellsToMove);

  /*!
    @brief Compute the classical diffusive time step
  */
  virtual Real compute_diffusive_dt();

  /*!
    @brief Get number of particles
  */
  virtual size_t get_num_particles(const ito_solver::which_container a_ptype, const bool a_local);

  /*!
    @brief Get solvers
  */
  virtual Vector<RefCountedPtr<T> >& get_solvers();

  /*!
    @brief Get species
  */
  virtual Vector<RefCountedPtr<ito_species> >& getSpecies();

  /*!
    @brief Get all velocity functions
  */
  virtual Vector<EBAMRCellData*> get_velo_funcs();

  /*!
    @brief Get all densities
  */
  virtual Vector<EBAMRCellData*> get_densities();

  /*!
    @brief Get all scratch storage
  */
  virtual Vector<EBAMRCellData*> get_scratch();

  /*!
    @brief Get all diffusion coefficients
  */
  virtual Vector<EBAMRCellData*> get_diffco_funcs();

  /*!
    @brief Get all mobility mesh functions
  */
  virtual Vector<EBAMRCellData*> get_mobility_funcs();

  /*!
    @brief Get particle containers
  */
  virtual Vector<ParticleContainer<ito_particle>* > getParticles(const ito_solver::which_container a_ptype);

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
  Vector<RefCountedPtr<ito_species> > m_species;

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
  @brief Factory class for ito_layout
*/
template <class T, class S>
class ito_factory {
public:
  ito_factory();
  ~ito_factory();
  RefCountedPtr<ito_layout<T> > new_layout(const Vector<RefCountedPtr<ito_species> > a_species) const;
};
#include "CD_NamespaceFooter.H"

#include "ito_layoutI.H"

#endif
