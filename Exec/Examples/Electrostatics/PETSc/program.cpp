#include <CD_Driver.H>
#include <CD_DiskProfiledPlane.H>
#include <CD_FieldStepper.H>
#include <CD_FieldSolverAMG.H>
#include <CD_FieldSolverGMG.H>

using namespace ChomboDischarge;
using namespace Physics::Electrostatics;

using FSOLVE = FieldSolverGMG;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new DiskProfiledPlane());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto timestepper = RefCountedPtr<FieldStepper<FSOLVE>>(new FieldStepper<FSOLVE>());
  auto driver      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr));

  driver->setupAndRun();

  ChomboDischarge::finalize();
}
