#include <CD_Driver.H>
#include <CD_CoaxialCable.H>
#include <CD_MeshODEStepper.H>

using namespace ChomboDischarge;
using namespace Physics::MeshODE;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new CoaxialCable());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto timestepper = RefCountedPtr<TimeStepper>(new MeshODEStepper<1>());
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr));

  engine->setupAndRun();

  ChomboDischarge::finalize();
}
