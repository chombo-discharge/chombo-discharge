/*!
  @file   rte_layout.H
  @brief  Declaration of a class that holds a set of rte_solvers
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _RTE_LAYOUT_
#define _RTE_LAYOUT_

#include "plasma_kinetics.H"
#include "rte_solver.H"

/*!
  @brief Class for holding a set of rte_solvers. 
*/
class rte_layout {
public:

  /*!
    @brief Full constructor
  */
  rte_layout(const RefCountedPtr<plasma_kinetics> a_plaskin);

  /*!
    @brief Destructor
  */
  virtual ~rte_layout();

  /*!
    @brief Allocate internal storage for solvers
  */
  virtual void allocate_internals();



  /*!
    @brief Set amr
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Set the computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom);

  /*!
    @brief Set the physical domain
  */
  virtual void set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom);

  /*!
    @brief Do a sanity check
  */
  virtual void sanity_check();

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng set_computational_geometry
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
    @brief Regrid method. 
  */
  virtual void regrid();

  /*!
    @brief Convenience function. Set source terms for all species. Mostly used for debugging. 
  */
  virtual void set_source(const EBAMRCellData& a_source);

  /*!
    @brief Convenience function. Set constant source terms for all species. Mostly used for debugging. 
  */
  virtual void set_source(const Real a_source);

  /*!
    @brief Convenience function. All solvers write plot files
  */
  virtual void write_plot_file();

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
  virtual Vector<RefCountedPtr<rte_solver> >& get_solvers();

  /*!
    @brief Get species
  */
  virtual Vector<RefCountedPtr<photon_group> >& get_photons();

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
    @brief Solvers
  */
  Vector<RefCountedPtr<rte_solver> > m_solvers;

  /*!
    @brief Species
  */
  Vector<RefCountedPtr<photon_group> > m_photons;

  /*!
    @brief Verbosity
  */
  int m_verbosity;
};

#endif
