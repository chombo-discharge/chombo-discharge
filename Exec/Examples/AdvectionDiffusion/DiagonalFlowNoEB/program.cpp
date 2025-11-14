#include <CD_Driver.H>
#include <CD_CdrCTU.H>
#include <CD_RegularGeometry.H>
#include <CD_AdvectionDiffusionStepper.H>
#include <CD_AdvectionDiffusionTagger.H>

using namespace ChomboDischarge;
using namespace Physics::AdvectionDiffusion;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);
  
  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new RegularGeometry());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto solver      = RefCountedPtr<CdrSolver>(new CdrCTU());
  auto timestepper = RefCountedPtr<AdvectionDiffusionStepper>(new AdvectionDiffusionStepper(solver));
  auto tagger      = RefCountedPtr<CellTagger>(new AdvectionDiffusionTagger(solver, amr));
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  timestepper->setVelocity([](const RealVect& a_pos) {
    return RealVect::Unit;
  });

  engine->setupAndRun();

  ChomboDischarge::finalize();
}
