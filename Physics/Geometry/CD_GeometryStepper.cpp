/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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

GeometryStepper::GeometryStepper() {}
GeometryStepper::~GeometryStepper() {}

void
GeometryStepper::setupSolvers()
{}

void
GeometryStepper::allocate()
{}

void
GeometryStepper::initialData()
{}

void
GeometryStepper::postInitialize()
{}

void
GeometryStepper::postCheckpointSetup()
{}

void
GeometryStepper::registerRealms()
{}

void
GeometryStepper::registerOperators()
{}

#ifdef CH_USE_HDF5
void
GeometryStepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
{}
#endif

#ifdef CH_USE_HDF5
void
GeometryStepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl)
{}
#endif

void
GeometryStepper::writePlotData(LevelData<EBCellFAB>& a_output, int& a_icomp, const int a_level) const
{}

int
GeometryStepper::getNumberOfPlotVariables() const
{
  return 0;
}

Vector<std::string>
GeometryStepper::getPlotVariableNames() const
{
  return Vector<std::string>();
}

Real
GeometryStepper::computeDt()
{
  return 0.0;
}

Real
GeometryStepper::advance(const Real a_dt)
{
  return std::numeric_limits<Real>::max();
}

void
GeometryStepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt)
{}

void
GeometryStepper::printStepReport()
{}

void
GeometryStepper::preRegrid(const int a_lmin, const int a_oldFinestLevel)
{}

void
GeometryStepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{}

void
GeometryStepper::postRegrid()
{}

#include <CD_NamespaceFooter.H>
