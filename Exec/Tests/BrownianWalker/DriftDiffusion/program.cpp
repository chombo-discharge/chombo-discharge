#include "CD_Driver.H"
#include <CD_ItoSolver.H>
#include <CD_RodDielectric.H>
#include <CD_BrownianWalkerStepper.H>
#include <CD_BrownianWalkerTagger.H>
#include "ParmParse.H"

using namespace ChomboDischarge;
using namespace Physics::BrownianWalker;

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
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry>(new RodDielectric());
  RefCountedPtr<AmrMesh>               amr      = RefCountedPtr<AmrMesh>(new AmrMesh());

  // Set up basic BrownianWalker
  RefCountedPtr<ItoSolver>   solver      = RefCountedPtr<ItoSolver>(new ItoSolver());
  RefCountedPtr<TimeStepper> timestepper = RefCountedPtr<TimeStepper>(new BrownianWalkerStepper(solver));
  RefCountedPtr<CellTagger>  tagger      = RefCountedPtr<CellTagger>(new BrownianWalkerTagger(solver, amr));

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
