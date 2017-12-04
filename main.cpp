/*!
  @file main.cpp
  @brief Main file for running the chombo-streamer code
  @author Robert Marskar
*/

#include <ParmParse.H>

#include "sphere_sphere_geometry.H"
#include "plasma_engine.H"
#include "plasma_kinetics.H"
#include "rte_solver.H"
#include "amr_mesh.H"
#include "poisson_multifluid_gmg.H"
#include "poisson_staircase_gmg.H"

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Physical domain, geometry, time stepper, amr, and plasma kinetics
  const RealVect probLo = -RealVect::Unit;
  const RealVect probHi =  RealVect::Unit;

  RefCountedPtr<physical_domain> physdom         = RefCountedPtr<physical_domain> (new physical_domain(probLo, probHi));
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new sphere_sphere_geometry());
  RefCountedPtr<plasma_kinetics> plaskin         = RefCountedPtr<plasma_kinetics>(NULL);
  RefCountedPtr<time_stepper> timestepper        = RefCountedPtr<time_stepper>(NULL);
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());

  // Set up the amr strategey
  amr->set_verbosity(10);                         // Set verbosity
  amr->set_coarsest_num_cells(64*IntVect::Unit);  // Set number of cells on coarsest level
  amr->set_max_amr_depth(2);                      // Set max amr depth
  amr->set_ebcf(false);                           // Tell amr to forget about EBCF. 
  amr->set_refinement_ratio(2);                   // Set refinement ratio
  amr->set_fill_ratio(1.0);                       // Set grid fill ratio
  amr->set_blocking_factor(8);                    // Set blocking factor
  amr->set_buffer_size(2);                        // Set buffer size
  amr->set_max_box_size(64);                      // Set max box size
  amr->set_redist_rad(1);                         // Set redistribution radius
  amr->set_eb_ghost(4);                           // Set EB ghost vectors
  amr->set_physical_domain(physdom);              // Set physical domain
  amr->set_irreg_sten_order(1);
  amr->set_irreg_sten_radius(1);

  // Set up plasma engine
  RefCountedPtr<plasma_engine> engine = RefCountedPtr<plasma_engine> (new plasma_engine(physdom,
											compgeom,
											plaskin,
											timestepper,
											amr));

  // Set up a Poisson solver
  //RefCountedPtr<poisson_solver> poisson = RefCountedPtr<poisson_solver> (new poisson_staircase_gmg());
  RefCountedPtr<poisson_solver> poisson = RefCountedPtr<poisson_solver> (new poisson_multifluid_gmg());
  poisson->set_verbosity(10);
  poisson->set_amr(amr);
  poisson->set_computational_geometry(compgeom);
  poisson->set_physical_domain(physdom);

  if(SpaceDim == 2){
    poisson->set_neumann_wall_bc(0,   Side::Lo, 0.0);                  
    poisson->set_neumann_wall_bc(0,   Side::Hi, 0.0);
    poisson->set_dirichlet_wall_bc(1, Side::Lo, Potential::Ground);
    poisson->set_dirichlet_wall_bc(1, Side::Hi, Potential::Live);
  }
  else if(SpaceDim == 3){
    poisson->set_neumann_wall_bc(0,   Side::Lo, 0.0);                  
    poisson->set_neumann_wall_bc(0,   Side::Hi, 0.0);
    poisson->set_neumann_wall_bc(1,   Side::Lo, 0.0);                  
    poisson->set_neumann_wall_bc(1,   Side::Hi, 0.0);
    poisson->set_dirichlet_wall_bc(2, Side::Lo, Potential::Ground);
    poisson->set_dirichlet_wall_bc(2, Side::Hi, Potential::Live);
  }
  

  // Setup plasma engine
  engine->set_verbosity(10);
  engine->set_geom_refinement_depth(2);
  engine->setup_fresh();

  // Poisson solver solves
  poisson->sanity_check();
  poisson->solve();
  

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
