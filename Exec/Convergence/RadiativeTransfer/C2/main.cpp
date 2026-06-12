/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   main.cpp
  @brief  Convergence test for radiative transfer on a uniform grid.
  @author Robert Marskar
*/

#include <CD_Driver.H>
#include <CD_EddingtonSP1.H>
#include <CD_CoaxialCable.H>
#include <CD_RadiativeTransferStepper.H>
#include <ParmParse.H>

using namespace ChomboDischarge;
using namespace Physics::RadiativeTransfer;

// This program runs convergence testing on a uniform grid (with EB) for a transient radiative transfer
// problem using the EddingtonSP1 simplification of the RTE. When there are
// solutions with different temporal resolutions available, we can compute the "error" as
//
//    e = phi(T, dt) - phi(T, dt/2)
//
// where phi(T,dt) indicates a solution advanced to time T with a time step dt.

int
main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse         pp(argc - 2, argv + 2, nullptr, input_file.c_str());

  // How much we refine the time step. numRefine = 1 => refine once => two runs. And so on.
  Real      dt;
  const int numRefine = 5;

  ParmParse pp2("RadiativeTransferStepper");
  pp2.get("dt", dt);

  // Storage for max, L1, and L2 solution error norms.
  std::vector<std::array<Real, 3>> norms;

  // Set geometry and AMR
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry>(new CoaxialCable());
  RefCountedPtr<AmrMesh>               amr      = RefCountedPtr<AmrMesh>(new AmrMesh());
  RefCountedPtr<CellTagger>            tagger   = RefCountedPtr<CellTagger>(nullptr);

  // Set up the time stepper.
  auto timestepper = RefCountedPtr<RadiativeTransferStepper<EddingtonSP1>>(
    new RadiativeTransferStepper<EddingtonSP1>());

  // Storage for the error.
  EBAMRCellData error;

  // Run simulations.
  for (int i = 0; i < numRefine; i++) {
    timestepper->forceDt(dt);

    // Run the simulation.
    RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));
    engine->setupAndRun(input_file);

    // Can compute error if i > 0
    if (i > 0) {
      DataOps::incr(error, timestepper->getPhi(), -1.0);

      // Computes max, L1, and L2 norms.
      const Real Linf = DataOps::norm(*error[0], 0, 0, *amr->getVofIterator("primal", phase::gas)[0]);
      const Real L1   = DataOps::norm(*error[0], 1, 0, *amr->getVofIterator("primal", phase::gas)[0]);
      const Real L2   = DataOps::norm(*error[0], 2, 0, *amr->getVofIterator("primal", phase::gas)[0]);

      norms.emplace_back(std::array<Real, 3>{Linf, L1, L2});
    }

    // Allocate storage for the fine/coarse solutions.
    amr->allocate(error, "primal", phase::gas, 1);
    error.copy(timestepper->getPhi());

    dt *= 0.5;
  }

  // Print the solution error.
#ifdef CH_MPI
  if (procID() == 0) {
#endif
    // clang-format off
    std::cout << "# nref\t"
              << "Linf error\t"
              << "L1 error\t"
              << "L2 error\n";

    for (size_t i = 0; i < norms.size(); i++) {
      std::cout << std::pow(2, i + 1) << "\t" << std::get<0>(norms[i]) << "\t" << std::get<1>(norms[i]) << "\t"
                << std::get<2>(norms[i]) << "\n";
    }
    // clang-format on
#ifdef CH_MPI
  }
  MPI_Finalize();
#endif
}
