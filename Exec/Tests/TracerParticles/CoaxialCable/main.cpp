#include <CD_Driver.H>
#include <CD_CoaxialCable.H>
#include <CD_TracerParticle.H>
#include <CD_TracerParticleStepper.H>

using namespace ChomboDischarge;
using namespace Physics::TracerParticle;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new CoaxialCable());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto timestepper = RefCountedPtr<TimeStepper>(new TracerParticleStepper<TracerParticle<0, 4>>());
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr));

  engine->setupAndRun();

  ChomboDischarge::finalize();
}
