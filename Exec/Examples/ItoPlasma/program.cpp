#include <CD_Driver.H>
#include <CD_GeoCoarsener.H>
#include <CD_FieldSolverFactory.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_ItoLayout.H>
#include <CD_ItoSolver.H>
#include <CD_RtLayout.H>
#include <CD_McPhoto.H>
#include <CD_ItoPlasmaAir3LFA.H>
#include <CD_RodDielectric.H>
#include <CD_ItoPlasmaGodunovStepper.H>
#include <CD_ItoPlasmaStreamerTagger.H>
#include "ParmParse.H"

// This is the potential curve (constant in this case). Modify it if you want to.
Real g_potential;
Real
potential_curve(const Real a_time)
{
  return g_potential;
}

using namespace ChomboDischarge;
using namespace Physics::ItoPlasma;

int
main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse         pp(argc - 2, argv + 2, NULL, input_file.c_str());

  // Get potential from input script
  std::string basename;
  {
    ParmParse pp("ItoPlasma");
    pp.get("potential", g_potential);
    pp.get("basename", basename);
    setPoutBaseName(basename);
  }

  // Set geometry and AMR
  RefCountedPtr<ComputationalGeometry> compgeom   = RefCountedPtr<ComputationalGeometry>(new RodDielectric());
  RefCountedPtr<AmrMesh>               amr        = RefCountedPtr<AmrMesh>(new AmrMesh());
  RefCountedPtr<GeoCoarsener>          geocoarsen = RefCountedPtr<GeoCoarsener>(new GeoCoarsener());

  // Set up physics
  RefCountedPtr<ItoPlasmaPhysics> physics     = RefCountedPtr<ItoPlasmaPhysics>(new ItoPlasmaAir3LFA());
  RefCountedPtr<ItoPlasmaStepper> timestepper = RefCountedPtr<ItoPlasmaStepper>(new ItoPlasmaGodunovStepper(physics));
  RefCountedPtr<CellTagger>       tagger      = RefCountedPtr<CellTagger>(
    new ItoPlasmaStreamerTagger(physics, timestepper, amr, compgeom));

  // Create solver factories
  auto poi_fact = new FieldSolverFactory<FieldSolverMultigrid>();
  auto ito_fact = new ItoFactory<ItoSolver, ItoSolver>();
  auto rte_fact = new RtFactory<McPhoto, McPhoto>();

  // Instantiate solvers
  auto poi = poi_fact->newSolver();
  auto cdr = ito_fact->newLayout(physics->getItoSpecies());
  auto rte = rte_fact->newLayout(physics->getRtSpecies());

  // Send solvers to TimeStepper
  timestepper->setFieldSolver(poi);
  timestepper->setIto(cdr);
  timestepper->setRadiativeTransferSolvers(rte);

  // Set potential
  timestepper->setVoltage(potential_curve);

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger, geocoarsen));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
