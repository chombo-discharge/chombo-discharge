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
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new sphere_sphere_geometry());
  RefCountedPtr<physical_domain> physdom         = RefCountedPtr<physical_domain> (new physical_domain(probLo, probHi));
  RefCountedPtr<plasma_kinetics> plaskin         = RefCountedPtr<plasma_kinetics>( NULL);
  RefCountedPtr<plasma_kinetics> timestepper     = RefCountedPtr<time_stepper>( NULL);


  // Set up plasma engine
  RefCountedPtr<plasma_engine> engine = RefCountedPtr<plasma_engine> (new plasma_engine(compgeom, plaskin, timestepper));
  engine->set_verbosity(10);
  engine->set_neumann_wall_bc(0,   Side::Lo, 0.0);
  engine->set_neumann_wall_bc(0,   Side::Hi, 0.0);
  engine->set_dirichlet_wall_bc(1, Side::Lo, Potential::Ground);
  engine->set_dirichlet_wall_bc(1, Side::Hi, Potential::Live);

  engine->set_physical_domain(physdom);
  
  engine->setup_fresh();

  // Real r;
  // ParmParse pp("sphere");
  // pp.get("electrode_radius", r);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
