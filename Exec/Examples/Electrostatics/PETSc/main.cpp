/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Example of FieldSolverAMG: a potential difference and a space charge blob in a box, on three AMR levels.
 * @author Robert Marskar
 */

// Std includes
#include <string>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_Driver.H>
#include <CD_RegularGeometry.H>
#include <CD_RefinementBoxTagger.H>
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
    std::string realm;

    ParmParse ppExample("PETScExample");
    ParmParse ppStepper("FieldStepper");

    ppExample.get("field_solver", fieldSolver);
    ppStepper.get("realm", realm);

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

    // There is no embedded boundary, so the refinement comes from the boxes given to RefinementBoxTagger.
    auto compgeom = RefCountedPtr<ComputationalGeometry>(new RegularGeometry());
    auto amr      = RefCountedPtr<AmrMesh>(new AmrMesh());
    auto tagger   = RefCountedPtr<CellTagger>(new RefinementBoxTagger(amr, realm));
    auto driver   = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

    driver->setupAndRun();
  }
#else
  MayDay::Abort("PETScExample -- this example requires USE_PETSC=TRUE");
#endif

  ChomboDischarge::finalize();
}
