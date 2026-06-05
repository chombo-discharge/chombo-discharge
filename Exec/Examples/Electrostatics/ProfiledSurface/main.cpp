#include <CD_Driver.H>
#include <CD_FieldSolverGMG.H>
#include <CD_FieldSolverAMG.H>
#include <CD_DiskProfiledPlane.H>
#include <CD_FieldStepper.H>

using namespace ChomboDischarge;
using namespace Physics::Electrostatics;

#ifdef CH_USE_PETSC
using F = FieldSolverAMG;
#else
using F = FieldSolverGMG;
#endif

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new DiskProfiledPlane());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto tagger      = RefCountedPtr<CellTagger>(nullptr);
  auto timestepper = RefCountedPtr<FieldStepper<F>>(new FieldStepper<F>());
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  engine->setupAndRun();

  ChomboDischarge::finalize();
}
