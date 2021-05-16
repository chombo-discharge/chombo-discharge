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
void geometry_stepper::initialData() {};
void geometry_stepper::postInitialize() {};
void geometry_stepper::postCheckpointSetup() {};

// Registration routines
void geometry_stepper::registerRealms(){}
void geometry_stepper::registerOperators(){ 
  m_amr->registerOperator(s_eb_fill_patch, Realm::Primal, phase::gas); // Need for averaging down ghost cells in Driver output routine
}

// IO routines
void geometry_stepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const {}
void geometry_stepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl) {}
void geometry_stepper::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const {}
int  geometry_stepper::getNumberOfPlotVariables() const {return 0;}

// Advance routines
void geometry_stepper::computeDt(Real& a_dt, TimeCode& a_timeCode) {a_dt = 0.0;}
Real geometry_stepper::advance(const Real a_dt) {return 1.0;}
void geometry_stepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) {}
void geometry_stepper::printStepReport() {}

// Regrid routines
void geometry_stepper::deallocate() {}
void geometry_stepper::preRegrid(const int a_lmin, const int a_oldFinestLevel) {}
void geometry_stepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) {}
void geometry_stepper::postRegrid() {}
#include "CD_NamespaceFooter.H"
