/*!
  @file main.cpp
  @brief Main file for running the chombo-streamer code
  @author Robert Marskar
*/

// Chombo files
#include "ParmParse.H"
#include "EBIndexSpace.H"

#include "sphere_sphere_geometry.H"
#include "plasma_engine.H"
#include "plasma_kinetics.H"
#include "rte_solver.H"
#include "amr_mesh.H"

//
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
  amr->set_max_amr_depth(3);                      // Set max amr depth
  amr->set_ebcf(false);                           // Tell amr to forget about EBCF
  amr->set_refinement_ratio(4);                   // Set refinement ratio
  amr->set_fill_ratio(1.0);                       // Set grid fill ratio
  amr->set_blocking_factor(4);                    // Set blocking factor
  amr->set_buffer_size(1);                        // Set buffer size
  amr->set_max_box_size(16);                      // Set max box size
  amr->set_redist_rad(1);                         // Set redistribution radius
  amr->set_num_ghost(2);                          // Set number of ghost cells (this is overridden inside plasma_engine)
  amr->set_eb_ghost(2);                           // Set EB ghost vectors
  amr->set_physical_domain(physdom);              // Set physical domain

  // Set up plasma engine
  RefCountedPtr<plasma_engine> engine = RefCountedPtr<plasma_engine> (new plasma_engine(physdom,
											compgeom,
											plaskin,
											timestepper,
											amr));
  engine->set_verbosity(10);
  engine->set_geom_refinement_depth(2);
  engine->set_neumann_wall_bc(0,   Side::Lo, 0.0);                  
  engine->set_neumann_wall_bc(0,   Side::Hi, 0.0);
  engine->set_dirichlet_wall_bc(1, Side::Lo, Potential::Ground);
  engine->set_dirichlet_wall_bc(1, Side::Hi, Potential::Live);
  
  engine->setup_fresh();

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
