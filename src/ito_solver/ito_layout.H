/*!
  @file   ito_layout.H
  @brief  Declaration of a class that holds a set of ito_solvers
  @author Robert Marskar
  @date   May 2020
*/

#ifndef _ITO_LAYOUT_
#define _ITO_LAYOUT_

#include "ito_solver.H"

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
    @brief Define function
  */
  virtual void define(const Vector<RefCountedPtr<ito_species> >& a_species);

  /*!
    @brief Parse options for all solvers
  */
  virtual void parse_options();

  /*!
    @brief Allocate internals for all solvers
  */
  virtual void allocate_internals();

  /*!
    @brief Add solver
  */
  virtual void add_solver(RefCountedPtr<T> a_solver);

  /*!
    @brief Cache states before regrid step
  */
  virtual void pre_regrid(const int a_lbase, const int a_finest_level);

  /*!
    @brief Initial data
  */
  virtual void initial_data();

  /*!
    @brief All solvers deposit their particles
  */
  virtual void deposit_particles();

  /*!
    @brief All solvers remap. 
  */
  virtual void remap();

  /*!
    @brief Interpolate velocities
  */
  virtual void interpolate_velocities();

  /*!
    @brief Regrid method. 
  */
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Register operators
  */
  virtual void register_operators();

  /*!
    @brief Set amr
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Set the computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom);

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng set_computational_geometry
  */
  virtual void set_phase(phase::which_phase a_phase);

  /*!
    @brief Set verbosity.
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Set time
  */
  virtual void set_time(const int a_step, const Real a_time, const Real a_dt);
  
  /*!
    @brief Sort the m_particles container by cell
  */
  virtual void sort_particles_by_cell();

  /*!
    @brief Sort the m_particles container by cell
  */
  virtual void sort_particles_by_patch();

  /*!
    @brief Make superparticles
  */
  virtual void make_superparticles(const int a_ppc);

  /*!
    @brief Compute dt
  */
  virtual Real compute_dt();

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
    @brief Get number of particles
  */
  virtual size_t get_num_particles(const bool a_local);

  /*!
    @brief Get number of particles in m_particles
  */
  virtual size_t get_num_eb_particles(const bool a_local);

  /*!
    @brief Get number of particles in m_particles
  */
  virtual size_t get_num_domain_particles(const bool a_local);

  /*!
    @brief Get number of particles in m_particles
  */
  virtual size_t get_num_source_particles(const bool a_local);

  /*!
    @brief Get solvers
  */
  virtual Vector<RefCountedPtr<T> >& get_solvers();

  /*!
    @brief Get species
  */
  virtual Vector<RefCountedPtr<ito_species> >& get_species();

  /*!
    @brief Get all velocities
  */
  virtual Vector<EBAMRCellData*> get_velocities();

  /*!
    @brief Get all densities
  */
  virtual Vector<EBAMRCellData*> get_densities();

  /*!
    @brief Get all diffusion coefficients
  */
  virtual Vector<EBAMRCellData*> get_diffusion();

  /*!
    @brief Get all particles
  */
  virtual Vector<particle_container<ito_particle>* > get_particles();

  /*!
    @brief Get EB particles
  */
  virtual Vector<particle_container<ito_particle>* > get_eb_particles();

  /*!
    @brief Get domain particles
  */
  virtual Vector<particle_container<ito_particle>* > get_domain_particles();

  /*!
    @brief Get source particles
  */
  virtual Vector<particle_container<ito_particle>* > get_source_particles();

protected:

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

#include "ito_layoutI.H"

#endif
