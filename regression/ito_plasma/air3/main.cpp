#include "driver.H"
#include "geo_coarsener.H"
#include "poisson_factoryI.H"
#include "poisson_multifluid_gmg.H"
#include "ito_layout.H"
#include "ito_solver.H"
#include "rte_layout.H"
#include "mc_photo.H"
#include "ito_plasma_air3.H"
#include "rod_sphere.H"
#include "ito_plasma_godunov.H"
#include "ito_plasma_streamer_tagger.H"
#include "ParmParse.H"

// This is the potential curve (constant in this case). Modify it if you want to.
Real g_potential;
Real potential_curve(const Real a_time){
  return g_potential;
}

using namespace physics::ito_plasma;

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  char* input_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, input_file);

  // Get potential from input script 
  std::string basename; 
  {
     ParmParse pp("air3");
     pp.get("potential", g_potential);
     pp.get("basename",  basename);
     setPoutBaseName(basename);
  }

  // Set geometry and AMR 
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new rod_sphere());
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());
  RefCountedPtr<geo_coarsener> geocoarsen        = RefCountedPtr<geo_coarsener> (new geo_coarsener());

  // Set up physics 
  RefCountedPtr<ito_plasma_physics> physics      = RefCountedPtr<ito_plasma_physics> (new ito_plasma_air3());
  RefCountedPtr<ito_plasma_stepper> timestepper  = RefCountedPtr<ito_plasma_stepper> (new ito_plasma_godunov(physics));
  RefCountedPtr<cell_tagger> tagger              = RefCountedPtr<cell_tagger> (new ito_plasma_streamer_tagger(physics, timestepper, amr, compgeom));

  // Create solver factories
  auto poi_fact = new poisson_factory<poisson_multifluid_gmg>();
  auto ito_fact = new ito_factory<ito_solver, ito_solver>();
  auto rte_fact = new rte_factory<mc_photo, mc_photo>();

  // Instantiate solvers
  auto poi = poi_fact->new_solver();
  auto cdr = ito_fact->new_layout(physics->get_ito_species());
  auto rte = rte_fact->new_layout(physics->get_rte_species());

  // Send solvers to time_stepper 
  timestepper->set_poisson(poi);
  timestepper->set_ito(cdr);
  timestepper->set_rte(rte);

  // Set potential 
timestepper->set_potential(potential_curve);

  // Set up the driver and run it
  RefCountedPtr<driver> engine = RefCountedPtr<driver> (new driver(compgeom, timestepper, amr, tagger, geocoarsen));
  engine->setup_and_run();

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
