#include "CD_Driver.H"
#include <CD_CdrGodunov.H>
#include <CD_CdrCTU.H>
#include <CD_RodDielectric.H>
#include <CD_AdvectionDiffusionStepper.H>
#include <CD_AdvectionDiffusionTagger.H>
#include "ParmParse.H"

using namespace ChomboDischarge;
using namespace Physics::AdvectionDiffusion;

int
main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse         pp(argc - 2, argv + 2, NULL, input_file.c_str());

  // Set geometry and AMR
  RefCountedPtr<ComputationalGeometry> compgeom   = RefCountedPtr<ComputationalGeometry>(new RodDielectric());
  RefCountedPtr<AmrMesh>               amr        = RefCountedPtr<AmrMesh>(new AmrMesh());

  // Set up basic AdvectionDiffusion
  auto solver      = RefCountedPtr<CdrSolver>(new CdrGodunov());
  auto timestepper = RefCountedPtr<AdvectionDiffusionStepper>(new AdvectionDiffusionStepper(solver));
  auto tagger      = RefCountedPtr<CellTagger>(new AdvectionDiffusionTagger(solver, amr));

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
