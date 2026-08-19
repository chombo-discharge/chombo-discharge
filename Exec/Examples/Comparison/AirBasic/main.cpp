#include <CD_Driver.H>
#include <CD_FieldSolverFactory.H>
#include <CD_FieldSolverGMG.H>
#include <CD_CdrLayoutImplem.H>
#include <CD_CdrCTU.H>
#include <CD_RtLayoutImplem.H>
#include <CD_ItoSolver.H>
#include <CD_McPhoto.H>
#include <CD_RodNeedleDisk.H>
#include <CD_CdrPlasmaJSON.H>
#include <CD_CdrPlasmaGodunovStepper.H>
#include <CD_CdrPlasmaStreamerTagger.H>
#include <CD_ItoKMCJSON.H>
#include <CD_ItoKMCGodunovStepper.H>
#include <CD_ItoKMCStreamerTagger.H>

using namespace ChomboDischarge;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  // Seed the RNG. Both models draw their initial particles through ParticleManagement, so with a
  // deterministic seed the two models start from precisely the same set of particles.
  Random::seed();

  // Which model to run, and the voltage applied to the needle.
  std::string model;
  Real        potential;

  ParmParse pp("Comparison");
  pp.get("model", model);
  pp.get("potential", potential);

  auto geometry = RefCountedPtr<ComputationalGeometry>(new RodNeedleDisk());
  auto amr      = RefCountedPtr<AmrMesh>(new AmrMesh());

  // The two models are assembled below -- everything downstream of this only sees a TimeStepper and a CellTagger.
  RefCountedPtr<TimeStepper> timestepper;
  RefCountedPtr<CellTagger>  tagger;

  if (model == "cdr_plasma") {
    using namespace Physics::CdrPlasma;

    auto physics    = RefCountedPtr<CdrPlasmaPhysics>(new CdrPlasmaJSON());
    auto cdrStepper = RefCountedPtr<CdrPlasmaStepper>(new CdrPlasmaGodunovStepper(physics));

    auto fieldFactory = new FieldSolverFactory<FieldSolverGMG>();
    auto cdrFactory   = new CdrFactory<CdrSolver, CdrCTU>();
    auto rteFactory   = new RtFactory<RtSolver, McPhoto>();

    auto poisson = fieldFactory->newSolver();
    auto cdr     = cdrFactory->newLayout(physics->getCdrSpecies());
    auto rte     = rteFactory->newLayout(physics->getRtSpecies());

    cdrStepper->setFieldSolver(poisson);
    cdrStepper->setCdrSolvers(cdr);
    cdrStepper->setRadiativeTransferSolvers(rte);

    delete fieldFactory;
    delete cdrFactory;
    delete rteFactory;

    cdrStepper->setVoltage([potential](const Real t) -> Real {
      return potential;
    });

    timestepper = cdrStepper;
    tagger      = RefCountedPtr<CellTagger>(new CdrPlasmaStreamerTagger(physics, cdrStepper, amr, geometry));
  }
  else if (model == "ito_kmc") {
    using namespace Physics::ItoKMC;

    auto physics    = RefCountedPtr<ItoKMCPhysics>(new ItoKMCJSON());
    auto itoStepper = RefCountedPtr<ItoKMCStepper<>>(new ItoKMCGodunovStepper<>(physics));

    itoStepper->setVoltage([potential](const Real t) -> Real {
      return potential;
    });

    timestepper = itoStepper;
    tagger      = RefCountedPtr<CellTagger>(new ItoKMCStreamerTagger<ItoKMCStepper<>>(physics, itoStepper, amr));
  }
  else {
    MayDay::Abort("Comparison.model must be either 'cdr_plasma' or 'ito_kmc'");
  }

  auto engine = RefCountedPtr<Driver>(new Driver(geometry, timestepper, amr, tagger));

  engine->setupAndRun();

  ChomboDischarge::finalize();
}
