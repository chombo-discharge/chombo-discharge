#include "CD_Driver.H"
#include <CD_CdrGodunov.H>
#include <CD_CdrMuscl.H>
#include "rod_dielectric.H"
#include <CD_AdvectionDiffusionStepper.H>
#include <CD_AdvectionDiffusionTagger.H>
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
  RefCountedPtr<CdrSolver> solver        = RefCountedPtr<CdrSolver>   (new CdrGodunov());
  RefCountedPtr<TimeStepper> timestepper = RefCountedPtr<TimeStepper> (new AdvectionDiffusionStepper(solver));
  RefCountedPtr<CellTagger> tagger       = RefCountedPtr<CellTagger>  (new AdvectionDiffusionTagger(solver, amr));

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr, tagger, geocoarsen));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
#include "CD_NamespaceFooter.H"
