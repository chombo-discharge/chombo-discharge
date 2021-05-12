#include "driver.H"
#include "cdr_gdnv.H"
#include "cdr_muscl.H"
#include "rod_dielectric.H"
#include "advection_diffusion_stepper.H"
#include "advection_diffusion_tagger.H"
#include "ParmParse.H"

using namespace ChomboDischarge;
using namespace physics::advection_diffusion;

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

  // Set up basic advection_diffusion 
  RefCountedPtr<cdr_solver> solver        = RefCountedPtr<cdr_solver>   (new cdr_gdnv());
  RefCountedPtr<time_stepper> timestepper = RefCountedPtr<time_stepper> (new advection_diffusion_stepper(solver));
  RefCountedPtr<cell_tagger> tagger       = RefCountedPtr<cell_tagger>  (new advection_diffusion_tagger(solver, amr));

  // Set up the driver and run it
  RefCountedPtr<driver> engine = RefCountedPtr<driver> (new driver(compgeom, timestepper, amr, tagger, geocoarsen));
  engine->setup_and_run(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
#include "CD_NamespaceFooter.H"
