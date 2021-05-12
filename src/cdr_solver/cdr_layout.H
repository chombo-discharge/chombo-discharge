/*!
  @file   cdr_layout.H
  @brief  Declaration of a class that holds a set of cdr_solvers
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _CDR_LAYOUT_
#define _CDR_LAYOUT_

#include "cdr_solver.H"

#include "CD_NamespaceHeader.H"

template <class T>
class cdr_iterator;

/*!
  @brief Class for holding a set of cdr_solvers. 
*/
template <class T>
class cdr_layout {
public:

  cdr_layout();

  /*!
    @brief Full constructor
  */
  cdr_layout(const Vector<RefCountedPtr<cdr_species> >& a_plaskin);

  /*!
    @brief Destructor
  */
  virtual ~cdr_layout();

  /*!
    @brief Get phase
  */
  virtual phase::which_phase get_phase() const;

  /*!
    @brief Get the iterator
  */
  virtual cdr_iterator<T> iterator();

  /*!
    @brief Get the Realm
  */
  virtual const std::string get_Realm() const;

  /*!
    @brief Set
  */
  virtual void set_Realm(const std::string a_Realm);

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
  virtual void add_solver(RefCountedPtr<T> a_solver);

  /*!
    @brief Allocate internal storage for solvers
  */
  virtual void allocate_internals();

  /*!
    @brief Deallocate solver internals
  */
  virtual void deallocate_internals();

  /*!
    @brief Cache states before regrid step
  */
  virtual void pre_regrid(const int a_lbase, const int a_finest_level);

  /*!
    @brief Initial data
  */
  virtual void initial_data();

  /*!
    @brief Regrid method. 
  */
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Register operators
  */
  virtual void registerOperators();

  /*!
    @brief Set amr
  */
  virtual void set_amr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set the computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom);

  /*!
    @brief Do a sanity check
  */
  virtual void sanityCheck();

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng set_computational_geometry
  */
  virtual void set_phase(phase::which_phase a_phase = phase::gas);

  /*!
    @brief Set phase
  */
  virtual void set_bc(cdr_bc::which_bc a_bc);

  /*!
    @brief Set verbosity.
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Set time
  */
  virtual void set_time(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Convenience function. Set velocity for all species. Mostly used for debugging.
  */
  virtual void set_velocity(const EBAMRCellData& a_velo);

  /*!
    @brief Convenience function. Set constant velocity for all species. Mostly used for debugging. 
  */
  virtual void set_velocity(const RealVect a_velo);

  /*!
    @brief Convenience function. Set diffusion coefficients for all species. Mostly used for debugging. 
  */
  virtual void set_diffco(const EBAMRFluxData& a_diffco, const EBAMRIVData& a_diffco_eb);

  /*!
    @brief Convenience function. Set constant diffusion coefficients for all species. Mostly used for debugging. 
  */
  virtual void set_diffco(const Real a_diffco);

  /*!
    @brief Convenience function. Set source terms for all species. Mostly used for debugging. 
  */
  virtual void set_source(const EBAMRCellData& a_source);

  /*!
    @brief Convenience function. Set constant source terms for all species. Mostly used for debugging. 
  */
  virtual void set_source(const Real a_source);

  /*!
    @brief Convenience function. Set the EB flux for all species. Mostly used for debugging. 
  */
  virtual void set_ebflux(const EBAMRIVData& a_ebflux);

  /*!
    @brief Convenience function.Set the EB flux for all species. Mostly used for debugging. 
  */
  virtual void set_ebflux(const Real a_ebflux);

  /*!
    @brief Convenience function. All solvers write plot files
  */
  virtual void write_plot_file();

  /*!
    @brief Get largest possible time step for advection part
  */
  virtual Real compute_advection_dt();

  /*!
    @brief Get largest possible time step for diffusion part
  */
  virtual Real compute_diffusion_dt();

  /*!
    @brief Get largest possible time step for diffusion part
  */
  virtual Real compute_advection_diffusion_dt();

  /*!
    @brief Get largest possible time step for source term
  */
  virtual Real compute_source_dt(const Real a_tolerance, const bool a_elec_only);

  /*!
    @brief Compute total charge in computational region
  */
  virtual Real compute_total_charge();

  /*!
    @brief Compute mass for every species
  */
  virtual Vector<Real> compute_mass();

  /*!
    @brief Compute charge for each species
  */
  virtual Vector<Real> compute_charge();

  /*!
    @brief Get all names
  */
  virtual Vector<std::string> get_names();

  /*!
    @brief Get solvers
  */
  virtual Vector<RefCountedPtr<T> >& get_solvers();

  /*!
    @brief Get species
  */
  virtual Vector<RefCountedPtr<cdr_species> >& get_species();

  /*!
    @brief Get all source terms
  */
  virtual Vector<EBAMRCellData*> get_sources();

  /*!
    @brief Get all states
  */
  virtual Vector<EBAMRCellData*> get_states();

  /*!
    @brief Get all velocities
  */
  virtual Vector<EBAMRCellData*> get_velocities();

  /*!
    @brief Get all diffusion coefficients
  */
  virtual Vector<EBAMRFluxData*> get_diffco_face();

  /*!
    @brief Get all diffusion coefficients
  */
  virtual Vector<EBAMRIVData*> get_diffco_eb();

  /*!
    @brief Get the EB flux
  */
  virtual Vector<EBAMRIVData*> get_ebflux();

  /*!
    @brief Get domain fluxes
  */
  virtual Vector<EBAMRIFData*> get_domainflux();

  /*!
    @brief Get mass redistribution
  */
  virtual bool get_mass_redist();

protected:

  /*!
    @brief Realm
  */
  std::string m_Realm;

  /*!
    @brief Solvers
  */
  Vector<RefCountedPtr<cdr_solver> > m_solvers;

  /*!
    @brief Species
  */
  Vector<RefCountedPtr<cdr_species> > m_species;

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
  cdr_bc::which_bc m_bc_type;
};

/*!
  @brief Factory class for cdr_layout. 
  @details Factory class is very simple; since we don't want to template cdr_layout2 we use a factory to instantiate solvers of 
  any cdr_solver inherited class, and then return a layout with the casted classes. That's about it. 
*/
template <class T, class S>
class cdr_factory {
public:
  cdr_factory();
  ~cdr_factory();
  RefCountedPtr<cdr_layout<T> > new_layout(const Vector<RefCountedPtr<cdr_species> > a_species) const;
};
#include "CD_NamespaceFooter.H"

#endif
