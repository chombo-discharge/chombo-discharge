#include <CD_Driver.H>
#include <CD_FieldSolverGMG.H>
#include <CD_ItoSolver.H>
#include <CD_McPhoto.H>
#include <CD_ItoKMCJSON.H>
#include <CD_CdrCTU.H>
#include <CD_WireWire.H>
#include <CD_ItoKMCGodunovStepper.H>
#include <CD_ItoKMCStreamerTagger.H>

using namespace ChomboDischarge;
using namespace Physics::ItoKMC;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  Random::seed();

  auto geometry    = RefCountedPtr<ComputationalGeometry>(new WireWire());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto physics     = RefCountedPtr<ItoKMCPhysics>(new ItoKMCJSON());
  auto timestepper = RefCountedPtr<ItoKMCStepper<>>(new ItoKMCGodunovStepper<>(physics));
  auto tagger      = RefCountedPtr<CellTagger>(new ItoKMCStreamerTagger<ItoKMCStepper<>>(physics, timestepper, amr));
  auto engine      = RefCountedPtr<Driver>(new Driver(geometry, timestepper, amr, tagger));

  Real U0 = 1.0;

  ParmParse pp("ItoKMC");
  pp.get("potential", U0);
  timestepper->setVoltage([=](const Real t) -> Real {
    return U0;
  });

  engine->setupAndRun();

  ChomboDischarge::finalize();
}
