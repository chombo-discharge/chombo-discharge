#include <CD_Driver.H>
#include <CD_ItoSolver.H>
#include <CD_RodDielectric.H>
#include <CD_BrownianWalkerStepper.H>
#include <CD_BrownianWalkerTagger.H>

using namespace ChomboDischarge;
using namespace Physics::BrownianWalker;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new RodDielectric());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto solver      = RefCountedPtr<ItoSolver>(new ItoSolver());
  auto timestepper = RefCountedPtr<TimeStepper>(new BrownianWalkerStepper(solver));
  auto tagger      = RefCountedPtr<CellTagger>(new BrownianWalkerTagger(solver, amr));
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  engine->setupAndRun();

  ChomboDischarge::finalize();
}
