#include "ParmParse.H"

#include <CD_Driver.H>
#include <CD_GeoCoarsener.H>
#include <CD_FieldSolverFactory.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_ItoLayout.H>
#include <CD_ItoSolver.H>
#include <CD_RtLayout.H>
#include <CD_McPhoto.H>
#include <CD_ItoPlasmaAir3LFA.H>
#include <CD_CdrCTU.H>
#include <CD_RodDielectric.H>
#include <CD_ItoPlasmaGodunovStepper.H>
#include <CD_ItoPlasmaStreamerTagger.H>

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

  Random::seed();

  // Get potential from input script
  std::string basename;
  {
    ParmParse pp("ItoPlasma");
    pp.get("potential", g_potential);
    pp.get("basename", basename);
    setPoutBaseName(basename);
  }

  auto geometry    = RefCountedPtr<ComputationalGeometry>(new RodDielectric());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto geocoarsen  = RefCountedPtr<GeoCoarsener>(new GeoCoarsener());
  auto physics     = RefCountedPtr<ItoPlasmaPhysics>(new ItoPlasmaAir3LFA());
  auto timestepper = RefCountedPtr<ItoPlasmaStepper<>>(new ItoPlasmaGodunovStepper<>(physics));
  auto tagger      = RefCountedPtr<CellTagger>(new ItoPlasmaStreamerTagger<ItoPlasmaStepper<>>(physics, timestepper, amr));

  // Set potential
  timestepper->setVoltage(potential_curve);

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(geometry, timestepper, amr, tagger, geocoarsen));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
