#include "CD_Driver.H"
#include <CD_GeoCoarsener.H>
#include <CD_FieldSolverFactory.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_CdrLayoutImplem.H>
#include <CD_CdrCTU.H>
#include <CD_RtLayoutImplem.H>
#include <CD_McPhoto.H>
#include <CD_CdrPlasmaJSON.H>
#include <CD_Vessel.H>
#include <CD_CdrPlasmaGodunovStepper.H>
#include <CD_CdrPlasmaStreamerTagger.H>
#include <ParmParse.H>

// This is the voltage curve (constant in this case). Modify it if you want to.
Real g_voltage;
Real
voltageCurve(const Real a_time)
{
  return g_voltage;
}

using namespace ChomboDischarge;
using namespace Physics::CdrPlasma;

int
main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse         pp(argc - 2, argv + 2, NULL, input_file.c_str());

  // Get voltage from input script
  std::string basename;
  {
    ParmParse pp("StochasticAir");
    pp.get("voltage", g_voltage);
    pp.get("basename", basename);
    setPoutBaseName(basename);
  }

  // Set geometry and AMR
  RefCountedPtr<ComputationalGeometry> compgeom   = RefCountedPtr<ComputationalGeometry>(new Vessel());
  RefCountedPtr<AmrMesh>               amr        = RefCountedPtr<AmrMesh>(new AmrMesh());
  RefCountedPtr<GeoCoarsener>          geocoarsen = RefCountedPtr<GeoCoarsener>(new GeoCoarsener());

  // Set up physics
  RefCountedPtr<CdrPlasmaPhysics> physics     = RefCountedPtr<CdrPlasmaPhysics>(new CdrPlasmaJSON());
  RefCountedPtr<CdrPlasmaStepper> timestepper = RefCountedPtr<CdrPlasmaStepper>(new CdrPlasmaGodunovStepper(physics));
  RefCountedPtr<CellTagger>       tagger      = RefCountedPtr<CellTagger>(
    new CdrPlasmaStreamerTagger(physics, timestepper, amr, compgeom));

  // Create solver factories
  auto poi_fact = new FieldSolverFactory<FieldSolverMultigrid>();
  auto cdr_fact = new CdrFactory<CdrSolver, CdrCTU>();
  auto rte_fact = new RtFactory<RtSolver, McPhoto>();

  // Instantiate solvers
  auto poi = poi_fact->newSolver();
  auto cdr = cdr_fact->newLayout(physics->getCdrSpecies());
  auto rte = rte_fact->newLayout(physics->getRtSpecies());

  // Send solvers to TimeStepper
  timestepper->setFieldSolver(poi);
  timestepper->setCdrSolvers(cdr);
  timestepper->setRadiativeTransferSolvers(rte);

  // Set voltage
  timestepper->setVoltage(voltageCurve);

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger, geocoarsen));
  engine->setupAndRun(input_file);

  // Clean up memory
  delete poi_fact;
  delete cdr_fact;
  delete rte_fact;

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
