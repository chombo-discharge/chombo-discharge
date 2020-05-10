/*!
  @file   poisson_stepperI.H
  @brief  Implementation of the poisson_stepper
  @author Robert Marskar
  @data   March 2020
*/

#include "poisson_stepper.H"
#include <ParmParse.H>

namespace physics {
  namespace poisson {
  
    template <typename T>
    Real poisson_stepper<T>::s_constant_one(const RealVect a_pos){
      return 1.0;
    }

    template <typename T>
    Real poisson_stepper<T>::s_potential_one(const Real a_time){
      return 1.0;
    }

    template <class T>
    void poisson_stepper<T>::setup_solvers(){
      int verb;
      ParmParse pp("poisson_stepper");
      pp.get("verbosity", verb);

      // Poisson setup
      m_poisson = RefCountedPtr<poisson_solver> (new T());
      m_poisson->set_verbosity(verb);
      m_poisson->parse_options();
      m_poisson->set_amr(m_amr);
      m_poisson->set_computational_geometry(m_compgeom);
      m_poisson->allocate_internals();
      m_poisson->set_potential(s_potential_one);

      // Setup sigma
      m_sigma = RefCountedPtr<sigma_solver> (new sigma_solver());
      m_sigma->set_verbosity(verb);
      m_sigma->set_amr(m_amr);
      m_sigma->set_computational_geometry(m_compgeom);
      m_sigma->allocate_internals();
    }

    template <class T>
    void poisson_stepper<T>::register_operators(){
      m_poisson->register_operators();
    }

    template <class T>
    void poisson_stepper<T>::initial_data(){
      ParmParse pp("poisson_stepper");
      pp.get("init_rho",   m_init_rho);
      pp.get("init_sigma", m_init_sigma);

      MFAMRCellData& state = m_poisson->get_state();
      MFAMRCellData& rho = m_poisson->get_source();
      EBAMRIVData& sigma = m_sigma->get_state();
  
      data_ops::set_value(rho,   m_init_rho);
      data_ops::set_value(sigma, m_init_sigma);

      const bool converged = m_poisson->solve(state, rho, sigma);

      if(!converged){
	MayDay::Warning("poisson_stepper<T>::advance - did not converge");
      }
    }

    template<class T>
    Real poisson_stepper<T>::advance(const Real a_dt){
      MayDay::Abort("poisson_stepper<T>::advance - callling this is an error. Please set driver.max_steps = 0");
      return 1.E99;
    }

    template <class T>
    void poisson_stepper<T>::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl){
      MayDay::Abort("poisson_stepper<T>::read_checkpoint_data - checkpointing not supported for this class");
    }

    template <class T>
    int poisson_stepper<T>::get_num_plot_vars() const{
      int ncomp = m_poisson->get_num_plotvars();
      ncomp += m_sigma->get_num_plotvars();
      return ncomp;
    }

    template <class T>
    void poisson_stepper<T>::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const{
      a_plotvar_names.append(m_poisson->get_plotvar_names());
      m_poisson->write_plot_data(a_output, a_icomp);

      a_plotvar_names.append(m_sigma->get_plotvar_names());
      m_sigma->write_plot_data(a_output, a_icomp);
    }

    template <class T>
    void poisson_stepper<T>::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt){
      m_step = a_step;
      m_time = a_time;
      m_dt   = a_dt;
    }

    template <class T>
    void poisson_stepper<T>::pre_regrid(const int a_lbase, const int a_old_finest_level){
      m_poisson->pre_regrid(a_lbase, a_old_finest_level);
      m_sigma->pre_regrid(a_lbase, a_old_finest_level);
    }

    template <class T>
    void poisson_stepper<T>::deallocate(){
      m_poisson->deallocate_internals();
      m_sigma->deallocate_internals();
    }

    template <class T>
    void poisson_stepper<T>::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
      m_poisson->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
      m_sigma->regrid(a_lmin,   a_old_finest_level, a_new_finest_level);

      m_poisson->solve();
    }
  }
}
