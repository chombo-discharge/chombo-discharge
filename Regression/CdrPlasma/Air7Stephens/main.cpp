#include "CD_Driver.H"
#include <CD_GeoCoarsener.H>
#include "CD_FieldSolverFactoryImplem.H"
#include "CD_FieldSolverMultigrid.H"
#include "cdr_layoutI.H"
#include <CD_CdrGodunov.H>
#include <CD_RtLayoutImplem.H>
#include <CD_McPhoto.H>
#include "air7_stephens.H"
#include "rod_dielectric.H"
#include "godunov.H"
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
    ParmParse pp("cdr_plasma");
    pp.get("potential", g_potential);
    pp.get("basename",  basename);
    setPoutBaseName(basename);
  }

  // Set geometry and AMR 
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry> (new rod_dielectric());
  RefCountedPtr<AmrMesh> amr                    = RefCountedPtr<AmrMesh> (new AmrMesh());
  RefCountedPtr<GeoCoarsener> geocoarsen        = RefCountedPtr<GeoCoarsener> (new GeoCoarsener());

  // Set up physics 
  RefCountedPtr<cdr_plasma_physics> physics      = RefCountedPtr<cdr_plasma_physics> (new air7_stephens());
  RefCountedPtr<cdr_plasma_stepper> timestepper  = RefCountedPtr<cdr_plasma_stepper> (new godunov(physics));
  RefCountedPtr<CellTagger> tagger              = RefCountedPtr<CellTagger> (new streamer_tagger(physics, timestepper, amr, compgeom));

  // Create solver factories
  auto poi_fact = new FieldSolverFactory<FieldSolverMultigrid>();
  auto cdr_fact = new cdr_factory<CdrSolver, CdrGodunov>();
  auto rte_fact = new RtFactory<RtSolver, McPhoto>();

  // Instantiate solvers
  auto poi = poi_fact->newSolver();
  auto cdr = cdr_fact->newLayout(physics->get_CdrSpecies());
  auto rte = rte_fact->newLayout(physics->get_RtSpecies());

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
