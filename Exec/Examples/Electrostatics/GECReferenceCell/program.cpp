#include <CD_Driver.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_GECReferenceCell.H>
#include <CD_FieldStepper.H>
#include <ParmParse.H>

using namespace ChomboDischarge;
using namespace Physics::Electrostatics;

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
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry>(new GECReferenceCell());
  RefCountedPtr<AmrMesh>               amr      = RefCountedPtr<AmrMesh>(new AmrMesh());
  RefCountedPtr<CellTagger>            tagger   = RefCountedPtr<CellTagger>(NULL);

  // Set up basic Poisson, potential = 1
  auto timestepper = RefCountedPtr<FieldStepper<FieldSolverMultigrid>>(new FieldStepper<FieldSolverMultigrid>());

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
