/*!
  @file   cdr_layout.H
  @brief  Declaration of a class that holds a set of cdr_solvers
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _CDR_LAYOUT_
#define _CDR_LAYOUT_

#include "plasma_kinetics.H"
#include "cdr_solver.H"

/*!
  @brief Class for holding a set of cdr_solvers. 
*/
class cdr_layout {
public:

  /*!
    @brief Full constructor
  */
  cdr_layout(const RefCountedPtr<plasma_kinetics> a_plaskin);

  /*!
    @brief Destructor
  */
  virtual ~cdr_layout();

  /*!
    @brief Get phase
  */
  virtual phase::which_phase get_phase() const;

  /*!
    @brief Allocate internal storage for solvers
  */
  virtual void allocate_internals();

  /*!
    @brief Initial data
  */
  virtual void initial_data();

  /*!
    @brief Regrid method. 
  */
  virtual void regrid(const int a_old_finest_level, const int a_new_finest_level);

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
    @brief Use mass-weighted redistribution
  */
  virtual void set_mass_redist(const bool a_mass_redist);

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
    @brief Convenience function. Call advance method for each solver
  */
  virtual void advance(const Real a_dt);

  /*!
    @brief Get largest possible time step for advection part
  */
  virtual Real compute_cfl_dt();

  /*!
    @brief Get largest possible time step for diffusion part
  */
  virtual Real compute_diffusive_dt();

  /*!
    @brief Compute total charge in computational region
  */
  virtual Real compute_Q();

  /*!
    @brief Get solvers
  */
  virtual Vector<RefCountedPtr<cdr_solver> >& get_solvers();

  /*!
    @brief Get species
  */
  virtual Vector<RefCountedPtr<species> >& get_species();

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

protected:

  /*!
    @brief Solvers
  */
  Vector<RefCountedPtr<cdr_solver> > m_solvers;

  /*!
    @brief Species
  */
  Vector<RefCountedPtr<species> > m_species;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief phase
  */
  phase::which_phase m_phase;
};

#endif
