#include "CD_Driver.H"
#include "geo_coarsener.H"
#include "CD_FieldSolverFactoryImplem.H"
#include "CD_FieldSolverMultigrid.H"
#include "cdr_layoutI.H"
#include <CD_CdrGodunov.H>
#include "rte_layoutI.H"
#include "eddington_sp1.H"
#include "air9eed_bourdon.H"
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
    ParmParse pp("air9eed_bourdon");
    pp.get("potential", g_potential);
    pp.get("basename",  basename);
    setPoutBaseName(basename);
  }

  // Set geometry and AMR 
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new rod_dielectric());
  RefCountedPtr<AmrMesh> amr                    = RefCountedPtr<AmrMesh> (new AmrMesh());
  RefCountedPtr<geo_coarsener> geocoarsen        = RefCountedPtr<geo_coarsener> (new geo_coarsener());

  // Set up physics 
  RefCountedPtr<cdr_plasma_physics> physics      = RefCountedPtr<cdr_plasma_physics> (new air9eed_bourdon());
  RefCountedPtr<cdr_plasma_stepper> timestepper  = RefCountedPtr<cdr_plasma_stepper> (new imex_sdc(physics));
  RefCountedPtr<CellTagger> tagger              = RefCountedPtr<CellTagger> (new streamer_tagger(physics, timestepper, amr, compgeom));

  // Create solver factories
  auto poi_fact = new FieldSolverFactory<FieldSolverMultigrid>();
  auto cdr_fact = new cdr_factory<CdrSolver, CdrGodunov>();
  auto rte_fact = new rte_factory<rte_solver, eddington_sp1>();

  // Instantiate solvers
  auto poi = poi_fact->newSolver();
  auto cdr = cdr_fact->new_layout(physics->get_cdr_species());
  auto rte = rte_fact->new_layout(physics->get_rte_species());

  // Send solvers to TimeStepper 
  timestepper->set_poisson(poi);
  timestepper->set_cdr(cdr);
  timestepper->set_rte(rte);

  // Set potential 
  timestepper->set_potential(potential_curve);

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr, tagger, geocoarsen));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
#include "CD_NamespaceFooter.H"
