/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Example of FieldSolverAMG: the ProfiledSurface problem, a live electrode above a profiled dielectric.
 * @author Robert Marskar
 */

// Std includes
#include <string>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_Driver.H>
#include <CD_DiskProfiledPlane.H>
#include <CD_FieldStepper.H>
#include <CD_FieldSolverAMG.H>
#include <CD_FieldSolverGMG.H>

using namespace ChomboDischarge;
using namespace Physics::Electrostatics;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

#ifdef CH_USE_PETSC
  // Scoped so that the PETSc objects are destroyed before ChomboDischarge::finalize finalizes PETSc and MPI.
  {
    std::string fieldSolver;

    ParmParse ppExample("PETScExample");

    ppExample.get("field_solver", fieldSolver);

    // The same problem can be solved with FieldSolverGMG, for comparing the plot files.
    RefCountedPtr<TimeStepper> timestepper;

    if (fieldSolver == "amg") {
      timestepper = RefCountedPtr<TimeStepper>(new FieldStepper<FieldSolverAMG>());
    }
    else if (fieldSolver == "gmg") {
      timestepper = RefCountedPtr<TimeStepper>(new FieldStepper<FieldSolverGMG>());
    }
    else {
      MayDay::Abort("PETScExample -- 'field_solver' must be 'amg' or 'gmg'");
    }

    // The mesh is refined only around the geometry (Driver.refine_electrodes, refine_dielectrics), so there is no
    // tagger.
    auto compgeom = RefCountedPtr<ComputationalGeometry>(new DiskProfiledPlane());
    auto amr      = RefCountedPtr<AmrMesh>(new AmrMesh());
    auto tagger   = RefCountedPtr<CellTagger>(nullptr);
    auto driver   = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

    driver->setupAndRun();
  }
#else
  MayDay::Abort("PETScExample -- this example requires USE_PETSC=TRUE");
#endif

  ChomboDischarge::finalize();
}
