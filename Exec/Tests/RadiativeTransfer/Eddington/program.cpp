#include "CD_Driver.H"
#include <CD_EddingtonSP1.H>
#include <CD_RodDielectric.H>
#include <CD_RadiativeTransferStepper.H>

using namespace ChomboDischarge;
using namespace Physics::RadiativeTransfer;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new RodDielectric());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto tagger      = RefCountedPtr<CellTagger>(NULL);
  auto timestepper = RefCountedPtr<RadiativeTransferStepper<EddingtonSP1>>(
    new RadiativeTransferStepper<EddingtonSP1>());
  auto engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  engine->setupAndRun();

  ChomboDischarge::finalize();
}
