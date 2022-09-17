#include <CD_Driver.H>
#include <CD_CoaxialCable.H>
#include <CD_TracerParticle.H>
#include <CD_TracerParticleStepper.H>
#include <ParmParse.H>

using namespace ChomboDischarge;
using namespace Physics::TracerParticle;

int
main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse         pp(argc - 2, argv + 2, NULL, input_file.c_str());

  // Set geometry and AMR
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry>(new CoaxialCable());
  RefCountedPtr<AmrMesh>               amr      = RefCountedPtr<AmrMesh>(new AmrMesh());

  // Set up basic tracer particle physics
  RefCountedPtr<TimeStepper> timestepper = RefCountedPtr<TimeStepper>(
    new TracerParticleStepper<TracerParticle<0, 4>>());

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
