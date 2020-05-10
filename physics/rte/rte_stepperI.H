/*!
  @file   rte_stepperI.H
  @brief  Implementation of rte_stepper.H
  @author Robert Marskar
  @data   May 2020
*/

#include "rte_stepper.H"
#include <ParmParse.H>

namespace physics {
  namespace rte {

    template <typename T>
    rte_stepper<T>::rte_stepper(){}
    
    template <typename T>
    rte_stepper<T>::~rte_stepper(){}

    template <typename T>
    void rte_stepper<T>::setup_solvers(){}

    template <typename T>
    void rte_stepper<T>::initial_data(){}

    template <typename T>
    void rte_stepper<T>::write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const {}

    template <typename T>
    void rte_stepper<T>::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl){}

    template <typename T>
    void rte_stepper<T>::post_checkpoint_setup() {}

    template <typename T>
    int rte_stepper<T>::get_num_plot_vars() const {return 0;}

    template <typename T>
    void rte_stepper<T>::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const {}

    template <typename T>
    void rte_stepper<T>::compute_dt(Real& a_dt, time_code::which_code& a_timecode) {}

    template <typename T>
    Real rte_stepper<T>::advance(const Real a_dt) {}

    template <typename T>
    void rte_stepper<T>::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) {}

    template <typename T>
    void rte_stepper<T>::print_step_report() {}

    // Regrid routines
    template <typename T>
    bool rte_stepper<T>::need_to_regrid() {}

    template <typename T>
    void rte_stepper<T>::register_operators() {}

    template <typename T>
    void rte_stepper<T>::cache() {}

    template <typename T>
    void rte_stepper<T>::deallocate() {}

    template <typename T>
    void rte_stepper<T>::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) {}
  }
}
