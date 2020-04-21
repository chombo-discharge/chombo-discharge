/*!
  @file   godunov.H
  @brief  Declaration of the Godunov time stepper method
  @author Robert Marskar
  @date   Aug. 2019
*/

#ifndef _GODUNOV_
#define _GODUNOV_

#include "type_definitions.H"
#include "amr_mesh.H"
#include "cdr_plasma_stepper.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"

namespace physics::cdr_plasma {
  
  /*!
    @brief Class for evolving plasma equations using the Euler-Maruyama scheme
  */
  class godunov : public cdr_plasma_stepper {
  public:
  
    // Forward declare nested classes. I don't want to clutter this file.
    class cdr_storage;
    class poisson_storage;
    class rte_storage;
    class sigma_storage;

    godunov();
    godunov(RefCountedPtr<cdr_plasma_physics>& a_physics);
    ~godunov();

    // Pure functions that must be overwritten
    Real advance(const Real a_dt);
    void parse_options();
    void init();
    void regrid_internals(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);
    void allocate_internals();
    void deallocate_internals();
    void compute_dt(Real& a_dt, time_code::which_code& a_timecode);
    bool need_to_regrid();

  protected:

    enum whichTransport{
      euler = 0,
      rk2   = 1,
    };

    enum whichDiffusion{
      Explicit  = 0,
      Implicit  = 1,
      Automatic = 2
    };

    // Scratch storage
    Vector<RefCountedPtr<cdr_storage> > m_cdr_scratch;
    Vector<RefCountedPtr<rte_storage> > m_rte_scratch;
    RefCountedPtr<poisson_storage>      m_poisson_scratch;
    RefCountedPtr<sigma_storage>        m_sigma_scratch;

    bool m_extrap_advect;
    bool m_debug;
    bool m_floor;
    bool m_fhd;
    bool m_implicit_diffusion;

    whichDiffusion m_whichDiffusion;
    whichTransport m_whichTransport;

    // Get functions for storage
    RefCountedPtr<cdr_storage>& get_cdr_storage(const cdr_iterator& a_solverit);
    RefCountedPtr<rte_storage>& get_rte_storage(const rte_iterator& a_solverit);

    Real restrict_dt();

    // Aux functions
    void compute_E_into_scratch();
    void compute_cdr_gradients();
    void compute_cdr_eb_states();
    void compute_cdr_eb_fluxes();
    void compute_cdr_domain_states();
    void compute_cdr_domain_fluxes();
    void compute_sigma_flux();
    void compute_reaction_network(const Real a_dt);

    void advance_transport(const Real a_dt);
    void advance_transport_euler(const Real a_dt);
    void advance_transport_rk2(const Real a_dt);
    void advance_rte(const Real a_dt);

    void compute_cdr_velo(const Real a_time);
    void compute_cdr_diffco(const Real a_time);

    // Various option parsing
    void parse_diffusion();
    void parse_advection();
    void parse_transport();
    void parse_floor();
    void parse_debug();
    void parse_fhd();
  };

}
#endif
