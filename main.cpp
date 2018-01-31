/*!
  @file main.cpp
  @brief Main file for running the chombo-streamer code
  @author Robert Marskar
*/

#include <ParmParse.H>

#include "plasma_engine.H"
#include "plasma_kinetics.H"
#include "rte_solver.H"
#include "amr_mesh.H"
#include "poisson_multifluid_gmg.H"
#include "poisson_staircase_gmg.H"
#include "cdr_solver.H"
#include "cdr_gdnv.H"
#include "cdr_sg.H"
#include "eddington_sp1.H"
#include "rk2.H"
#include "sphere_sphere_geometry.H"
#include "mechanical_shaft.H"
#include "cdr_layout.H"
#include "rte_layout.H"
#include "morrow_lowke.H"
#include "sigma_solver.H"

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Build argument list from input file
  char* inputFile = argv[1];
  ParmParse PP(argc-2,argv+2,NULL,inputFile);

  EBIndexSpace::s_useMemoryLoadBalance = false;

  // Physical domain, geometry, time stepper, amr, and plasma kinetics
  const RealVect probLo = -RealVect::Unit;
  const RealVect probHi =  RealVect::Unit;


  RefCountedPtr<physical_domain> physdom         = RefCountedPtr<physical_domain> (new physical_domain(probLo, probHi));
  RefCountedPtr<plasma_kinetics> plaskin         = RefCountedPtr<plasma_kinetics>(new morrow_lowke());
  RefCountedPtr<time_stepper> timestepper        = RefCountedPtr<time_stepper>(new time_stepper());
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());
#if CH_SPACEDIM == 2
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new sphere_sphere_geometry());
#else
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new mechanical_shaft());
#endif


  // Refinement ratios
  Vector<int> refrat(5);
  refrat[0] = 4;
  refrat[1] = 2;
  refrat[2] = 2;
  refrat[3] = 2;
  refrat[4] = 2;
  
  // Set up the amr strategey
  amr->set_verbosity(10);                         // Set verbosity
  amr->set_coarsest_num_cells(128*IntVect::Unit);  // Set number of cells on coarsest level
  amr->set_max_amr_depth(0);                      // Set max amr depth
  amr->set_ebcf(false);                           // Tell amr to forget about EBCF.
  amr->set_refinement_ratios(refrat);             // Set refinement ratios
  amr->set_fill_ratio(1.0);                       // Set grid fill ratio
  amr->set_blocking_factor(16);                   // Set blocking factor
  amr->set_buffer_size(1);                        // Set buffer size
  amr->set_max_box_size(32);                      // Set max box size
  amr->set_redist_rad(1);                         // Set redistribution radius
  amr->set_eb_ghost(4);                           // Set EB ghost vectors
  amr->set_physical_domain(physdom);              // Set physical domain
  amr->set_irreg_sten_type(stencil_type::linear); // Set preferred stencil type
  amr->set_irreg_sten_order(1);                   // Set preferred stencil order
  amr->set_irreg_sten_radius(1);                  // Set extrapolation stencil radius

  // Set up the plasma engine
  RefCountedPtr<plasma_engine> engine = RefCountedPtr<plasma_engine> (new plasma_engine(physdom,
											compgeom,
											plaskin,
											timestepper,
											amr));

  engine->set_verbosity(10);
  engine->set_geom_refinement_depth(-1);
  engine->setup_fresh();


#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
