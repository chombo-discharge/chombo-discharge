#include <CD_Driver.H>
#include <CD_CdrCTU.H>
#include <CD_RodDielectric.H>
#include <CD_AdvectionDiffusionStepper.H>
#include <CD_AdvectionDiffusionTagger.H>
#include <CD_DischargeIO.H>
#include <ParmParse.H>

using namespace ChomboDischarge;
using namespace Physics::AdvectionDiffusion;

// This program runs convergence testing on a uniform grid for an advection-diffusion problem. When there are
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
  ParmParse         pp(argc - 2, argv + 2, NULL, input_file.c_str());

  // How much we refine the time step. numRefine = 1 => refine once => two runs. And so on.
  Real      cfl       = 0.8;
  const int numRefine = 5;

  // Storage for max, L1, and L2 solution error norms.
  std::vector<std::array<Real, 3>> norms;

  // Set geometry and AMR
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry>(new RodDielectric());
  RefCountedPtr<AmrMesh>               amr      = RefCountedPtr<AmrMesh>(new AmrMesh());

  // Set up basic AdvectionDiffusion
  auto solver      = RefCountedPtr<CdrSolver>(new CdrCTU());
  auto timestepper = RefCountedPtr<AdvectionDiffusionStepper>(new AdvectionDiffusionStepper(solver));
  auto tagger      = RefCountedPtr<CellTagger>(new AdvectionDiffusionTagger(solver, amr));

  // Storage for the error.
  EBAMRCellData error;

  // Run simulations.
  for (int i = 0; i < numRefine; i++) {
    timestepper->setCFL(cfl);

    // Run the simulation.
    RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));
    engine->setupAndRun(input_file);

    // Can compute error if i > 0
    if (i > 0) {
      DataOps::incr(error, solver->getPhi(), -1.0);
      DischargeIO::writeEBHDF5(error, "error.hdf5");

      // Computes max, L1, and L2 norms.
      const Real Linf = DataOps::norm(*error[0], amr->getDomains()[0], 0, true);
      const Real L1   = DataOps::norm(*error[0], amr->getDomains()[0], 1, true);
      const Real L2   = DataOps::norm(*error[0], amr->getDomains()[0], 2, true);

      norms.emplace_back(std::array<Real, 3>{Linf, L1, L2});
    }

    // Allocate storage for the fine/coarse solutions.
    amr->allocate(error, "primal", phase::gas, 1);
    error.copy(solver->getPhi());

    cfl *= 0.5;
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

    for (int i = 0; i < norms.size(); i++) {
      std::cout << std::pow(2, i + 1) << "\t" << std::get<0>(norms[i]) << "\t" << std::get<1>(norms[i]) << "\t"
                << std::get<2>(norms[i]) << "\n";
    }
    // clang-format on
#ifdef CH_MPI
  }
  MPI_Finalize();
#endif
}
