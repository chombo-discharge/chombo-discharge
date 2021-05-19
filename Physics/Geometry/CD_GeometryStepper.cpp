/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_GeometryStepper.cpp
  @brief  Implementation of CD_GeometryStepper.H
  @author Robert Marskar
*/

// Our includes
#include <CD_GeometryStepper.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::Geometry;

// ctor/dtor
GeometryStepper::GeometryStepper(){}
GeometryStepper::~GeometryStepper(){};

// Setup routines
void GeometryStepper::setup_solvers() {};
void GeometryStepper::allocate() {};
void GeometryStepper::initialData() {};
void GeometryStepper::postInitialize() {};
void GeometryStepper::postCheckpointSetup() {};

// Registration routines
void GeometryStepper::registerRealms(){}
void GeometryStepper::registerOperators(){ 
  m_amr->registerOperator(s_eb_fill_patch, Realm::Primal, phase::gas); // Need for averaging down ghost cells in Driver output routine
}

// IO routines
void GeometryStepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const {}
void GeometryStepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl) {}
void GeometryStepper::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const {}
int  GeometryStepper::getNumberOfPlotVariables() const {return 0;}

// Advance routines
void GeometryStepper::computeDt(Real& a_dt, TimeCode& a_timeCode) {a_dt = 0.0;}
Real GeometryStepper::advance(const Real a_dt) {return 1.0;}
void GeometryStepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) {}
void GeometryStepper::printStepReport() {}

// Regrid routines
void GeometryStepper::preRegrid(const int a_lmin, const int a_oldFinestLevel) {}
void GeometryStepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) {}
void GeometryStepper::postRegrid() {}

#include <CD_NamespaceFooter.H>
