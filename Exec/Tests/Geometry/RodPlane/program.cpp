#include <CD_Driver.H>
#include <CD_RodPlane.H>
#include <CD_GeometryStepper.H>

using namespace ChomboDischarge;
using namespace Physics::Geometry;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new RodPlane());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto tagger      = RefCountedPtr<CellTagger>(nullptr);
  auto timestepper = RefCountedPtr<GeometryStepper>(new GeometryStepper());
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  engine->setupAndRun();

  ChomboDischarge::finalize();
}
