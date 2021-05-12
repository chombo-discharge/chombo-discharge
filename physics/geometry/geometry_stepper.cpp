/*!
  @file   geometry_stepper.cpp
  @brief  Implementation of geometry_stepper.H
  @author Robert Marskar
  @date   March 2020
*/

#include "geometry_stepper.H"

#include "CD_NamespaceHeader.H"
using namespace physics::geometry;

// ctor/dtor
geometry_stepper::geometry_stepper(){}
geometry_stepper::~geometry_stepper(){};

// Setup routines
void geometry_stepper::setup_solvers() {};
void geometry_stepper::allocate() {};
void geometry_stepper::initial_data() {};
void geometry_stepper::post_initialize() {};
void geometry_stepper::post_checkpoint_setup() {};

// Registration routines
void geometry_stepper::registerRealms(){}
void geometry_stepper::registerOperators(){ 
  m_amr->registerOperator(s_eb_fill_patch, realm::primal, phase::gas); // Need for averaging down ghost cells in driver output routine
}

// IO routines
void geometry_stepper::write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const {}
void geometry_stepper::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) {}
void geometry_stepper::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const {}
int  geometry_stepper::get_num_plot_vars() const {return 0;}

// Advance routines
void geometry_stepper::compute_dt(Real& a_dt, time_code& a_timecode) {a_dt = 0.0;}
Real geometry_stepper::advance(const Real a_dt) {return 1.0;}
void geometry_stepper::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) {}
void geometry_stepper::print_step_report() {}

// Regrid routines
void geometry_stepper::deallocate() {}
void geometry_stepper::pre_regrid(const int a_lmin, const int a_old_finest_level) {}
void geometry_stepper::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) {}
void geometry_stepper::post_regrid() {}
#include "CD_NamespaceFooter.H"
