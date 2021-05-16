#include "CD_Driver.H"
#include "ito_solver.H"
#include "rod_dielectric.H"
#include "brownian_walker_stepper.H"
#include "brownian_walker_tagger.H"
#include "ParmParse.H"

using namespace ChomboDischarge;
using namespace physics::brownian_walker;

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, input_file.c_str());

  // Set geometry and AMR 
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new rod_dielectric());
  RefCountedPtr<AmrMesh> amr                    = RefCountedPtr<AmrMesh> (new AmrMesh());
  RefCountedPtr<geo_coarsener> geocoarsen        = RefCountedPtr<geo_coarsener> (new geo_coarsener());

  // Set up basic brownian_walker 
  RefCountedPtr<ito_solver> solver        = RefCountedPtr<ito_solver>   (new ito_solver());
  RefCountedPtr<TimeStepper> timestepper = RefCountedPtr<TimeStepper> (new brownian_walker_stepper(solver));
  RefCountedPtr<cell_tagger> tagger       = RefCountedPtr<cell_tagger>  (new brownian_walker_tagger(solver, amr));

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr, tagger, geocoarsen));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
#include "CD_NamespaceFooter.H"
