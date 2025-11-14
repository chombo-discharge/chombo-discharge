#include "CD_Driver.H"
#include <CD_FieldSolverFactory.H>
#include <CD_FieldSolverGMG.H>
#include <CD_CdrLayoutImplem.H>
#include <CD_CdrCTU.H>
#include <CD_RtLayoutImplem.H>
#include <CD_McPhoto.H>
#include <CD_CdrPlasmaJSON.H>
#include <CD_RegularGeometry.H>
#include <CD_CdrPlasmaGodunovStepper.H>
#include <CD_CdrPlasmaStreamerTagger.H>

using namespace ChomboDischarge;
using namespace Physics::CdrPlasma;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  Random::seed();

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new RegularGeometry());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto physics     = RefCountedPtr<CdrPlasmaPhysics>(new CdrPlasmaJSON());
  auto timestepper = RefCountedPtr<CdrPlasmaStepper>(new CdrPlasmaGodunovStepper(physics));
  auto tagger      = RefCountedPtr<CellTagger>(new CdrPlasmaStreamerTagger(physics, timestepper, amr, compgeom));
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  auto poi_fact = new FieldSolverFactory<FieldSolverGMG>();
  auto cdr_fact = new CdrFactory<CdrSolver, CdrCTU>();
  auto rte_fact = new RtFactory<RtSolver, McPhoto>();

  auto poi = poi_fact->newSolver();
  auto cdr = cdr_fact->newLayout(physics->getCdrSpecies());
  auto rte = rte_fact->newLayout(physics->getRtSpecies());

  timestepper->setFieldSolver(poi);
  timestepper->setCdrSolvers(cdr);
  timestepper->setRadiativeTransferSolvers(rte);
  timestepper->setVoltage([&](const Real a_time) -> Real {
    return 1.0;
  });

  engine->setupAndRun();

  // Clean up memory
  delete poi_fact;
  delete cdr_fact;
  delete rte_fact;

  ChomboDischarge::finalize();
}
