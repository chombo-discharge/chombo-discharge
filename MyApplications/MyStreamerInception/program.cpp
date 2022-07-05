#include <CD_Driver.H>
#include <CD_TracerParticleSolver.H>
#include <CD_Cylinder.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_StreamerInceptionStepper.H>
#include <ParmParse.H>

using namespace ChomboDischarge;
using namespace Physics::StreamerInception;

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, input_file.c_str());

  // Set geometry and AMR 
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry> (new Cylinder());
  RefCountedPtr<AmrMesh> amr                    = RefCountedPtr<AmrMesh> (new AmrMesh());

  // Set up time stepper 
  auto timestepper = RefCountedPtr<StreamerInceptionStepper<>> (new StreamerInceptionStepper<>());

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
