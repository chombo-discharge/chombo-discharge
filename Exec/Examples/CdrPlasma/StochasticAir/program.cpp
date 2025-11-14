#include <CD_Driver.H>
#include <CD_FieldSolverFactory.H>
#include <CD_FieldSolverGMG.H>
#include <CD_CdrLayoutImplem.H>
#include <CD_CdrCTU.H>
#include <CD_RtLayoutImplem.H>
#include <CD_McPhoto.H>
#include <CD_CdrPlasmaJSON.H>
#include <CD_Vessel.H>
#include <CD_CdrPlasmaGodunovStepper.H>
#include <CD_CdrPlasmaStreamerTagger.H>

using namespace ChomboDischarge;
using namespace Physics::CdrPlasma;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new Vessel());
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

  Real U0 = 1.0;

  ParmParse pp("StochasticAir");
  pp.get("voltage", U0);

  timestepper->setFieldSolver(poi);
  timestepper->setCdrSolvers(cdr);
  timestepper->setRadiativeTransferSolvers(rte);
  timestepper->setVoltage([=](const Real t) -> Real {
    return U0;
  });

  engine->setupAndRun();

  delete poi_fact;
  delete cdr_fact;
  delete rte_fact;

  ChomboDischarge::finalize();
}
