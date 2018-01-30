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
  @brief Class for holding a set of cdr_solvers. Comes with an iterator
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
    @brief Regrid method. 
  */
  virtual void regrid();

  /*!
    @brief Get solvers
  */
  virtual Vector<RefCountedPtr<cdr_solver> >& get_solvers();

  /*!
    @brief Get species
  */
  virtual Vector<RefCountedPtr<species> >& get_species();

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
};

#endif
