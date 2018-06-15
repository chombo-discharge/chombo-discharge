#include "plasma_engine.H"
#include "morrow_lowke.H"
#include "double_rod.H"
#include "rk2.H"
#include "rod_coarsen.H"
#include "ml_tagger.H"
#include "ParmParse.H"

// This is the potential curve (constant in this case). Modify it if you want to.
Real g_potential;
Real potential_curve(const Real a_time){
  return g_potential;
}

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  char* input_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, input_file);

  { // Get potential from input script 
    ParmParse pp("double_rod");
    pp.get("potential", g_potential);
  }

  // Set up everything 
  RefCountedPtr<plasma_kinetics> plaskin         = RefCountedPtr<plasma_kinetics> (new morrow_lowke());
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new double_rod());
  RefCountedPtr<time_stepper> timestepper        = RefCountedPtr<time_stepper> (new rk2());
  RefCountedPtr<cell_tagger> tagger              = RefCountedPtr<cell_tagger> (new ml_tagger());
  RefCountedPtr<physical_domain> physdom         = RefCountedPtr<physical_domain> (new physical_domain());
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());
  RefCountedPtr<geo_coarsener> geocoarsen        = RefCountedPtr<geo_coarsener> (new geo_coarsener());
  RefCountedPtr<plasma_engine> engine            = RefCountedPtr<plasma_engine> (new plasma_engine(physdom, compgeom, plaskin, timestepper, amr, tagger, geocoarsen));

  // Run the plasma engine
  engine->set_potential(potential_curve);
  engine->setup_and_run();

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
