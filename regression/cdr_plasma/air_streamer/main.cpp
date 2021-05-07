#include "driver.H"
#include "geo_coarsener.H"
#include "field_solver_factoryI.H"
#include "field_solver_multigrid.H"
#include "cdr_layoutI.H"
#include "cdr_gdnv.H"
#include "rte_layoutI.H"
#include "eddington_sp1.H"
#include "air3_bourdon.H"
#include "rod_dielectric.H"
#include "imex_sdc.H"
#include "streamer_tagger.H"
#include "ParmParse.H"

// This is the potential curve (constant in this case). Modify it if you want to.
Real g_potential;
Real potential_curve(const Real a_time){
  return g_potential;
}

using namespace ChomboDischarge;
using namespace physics::cdr_plasma;

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, input_file.c_str());

  // Get potential from input script 
  std::string basename; 
  {
     ParmParse pp("air_streamer");
     pp.get("potential", g_potential);
     pp.get("basename",  basename);
     setPoutBaseName(basename);
  }

  // Set geometry and AMR 
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new rod_dielectric());
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());
  RefCountedPtr<geo_coarsener> geocoarsen        = RefCountedPtr<geo_coarsener> (new geo_coarsener());

  // Set up physics 
  RefCountedPtr<cdr_plasma_physics> physics      = RefCountedPtr<cdr_plasma_physics> (new air3_bourdon());
  RefCountedPtr<cdr_plasma_stepper> timestepper  = RefCountedPtr<cdr_plasma_stepper> (new imex_sdc(physics));
  RefCountedPtr<cell_tagger> tagger              = RefCountedPtr<cell_tagger> (new streamer_tagger(physics, timestepper, amr, compgeom));

  // Create solver factories
  auto poi_fact = new field_solver_factory<field_solver_multigrid>();
  auto cdr_fact = new cdr_factory<cdr_solver, cdr_gdnv>();
  auto rte_fact = new rte_factory<rte_solver, eddington_sp1>();

  // Instantiate solvers
  auto poi = poi_fact->new_solver();
  auto cdr = cdr_fact->new_layout(physics->get_cdr_species());
  auto rte = rte_fact->new_layout(physics->get_rte_species());

  // Send solvers to time_stepper 
  timestepper->set_poisson(poi);
  timestepper->set_cdr(cdr);
  timestepper->set_rte(rte);

  // Set potential 
timestepper->set_potential(potential_curve);

  // Set up the driver and run it
  RefCountedPtr<driver> engine = RefCountedPtr<driver> (new driver(compgeom, timestepper, amr, tagger, geocoarsen));
  engine->setup_and_run(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
