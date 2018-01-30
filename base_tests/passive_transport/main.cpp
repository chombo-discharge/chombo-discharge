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
  RefCountedPtr<time_stepper> timestepper        = RefCountedPtr<time_stepper>(NULL);
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());
#if CH_SPACEDIM == 2
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new sphere_sphere_geometry());
#else
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new mechanical_shaft());
#endif

  Vector<int> refrat(5);
  refrat[0] = 4;
  refrat[1] = 2;
  refrat[2] = 2;
  refrat[3] = 2;
  refrat[4] = 2;
  
  // Set up the amr strategey
  amr->set_verbosity(10);                         // Set verbosity
  amr->set_coarsest_num_cells(64*IntVect::Unit);  // Set number of cells on coarsest level
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

  // Set up plasma engine. This builds the geometries and initializes amr so it must always be done. 
  RefCountedPtr<plasma_engine> engine = RefCountedPtr<plasma_engine> (new plasma_engine(physdom,
											compgeom,
											plaskin,
											timestepper,
											amr));
  engine->set_verbosity(10);
  engine->set_geom_refinement_depth(-1);
  engine->setup_fresh();

  // Set up a multifluid Poisson solver. This can also be done via a time stepper but we do it directly ehre. 
  RefCountedPtr<poisson_solver> poisson = RefCountedPtr<poisson_solver> (new poisson_multifluid_gmg());
  poisson->set_verbosity(10);
  poisson->set_amr(amr);
  poisson->set_computational_geometry(compgeom);
  poisson->set_physical_domain(physdom);

  if(SpaceDim == 2){
    poisson->set_neumann_wall_bc(0,   Side::Lo, 0.0);                  
    poisson->set_neumann_wall_bc(0,   Side::Hi, 0.0);
    poisson->set_dirichlet_wall_bc(0, Side::Lo, potential::ground);
    poisson->set_dirichlet_wall_bc(0, Side::Hi, potential::ground);
    poisson->set_dirichlet_wall_bc(1, Side::Lo, potential::ground);
    poisson->set_dirichlet_wall_bc(1, Side::Hi, potential::ground);
  }
  else if(SpaceDim == 3){
    poisson->set_neumann_wall_bc(0,   Side::Lo, 0.0);                  
    poisson->set_neumann_wall_bc(0,   Side::Hi, 0.0);
    poisson->set_neumann_wall_bc(1,   Side::Lo, 0.0);                  
    poisson->set_neumann_wall_bc(1,   Side::Hi, 0.0);
    poisson->set_dirichlet_wall_bc(2, Side::Lo, potential::ground);
    poisson->set_dirichlet_wall_bc(2, Side::Hi, potential::live);
  }

  // Poisson solver solves so that we can extract velocities. Compute the electric field. 
  poisson->sanity_check();
  poisson->allocate_internals();
  poisson->solve();
  poisson->write_plot_file();
  MFAMRCellData E;
  EBAMRCellData E_gas;
  amr->allocate(E, SpaceDim);
  amr->allocate_ptr(E_gas);
  amr->compute_gradient(E, poisson->get_state());
  amr->alias(E_gas, phase::gas, E);
  data_ops::scale(E_gas, -1.0);


  // Advance a layout of cdr and rte solvers
  RefCountedPtr<cdr_layout> cdr_solvers = RefCountedPtr<cdr_layout> (new cdr_layout(plaskin));

  cdr_solvers->set_amr(amr);
  cdr_solvers->set_computational_geometry(compgeom);
  cdr_solvers->set_physical_domain(physdom);
  cdr_solvers->sanity_check();
  cdr_solvers->allocate_internals();
  cdr_solvers->initial_data();
  cdr_solvers->set_velocity(E_gas);
  cdr_solvers->set_diffco(0.002);
  cdr_solvers->set_source(0.0);
  cdr_solvers->set_ebflux(0.0);
  cdr_solvers->write_plot_file();

  RefCountedPtr<rte_layout> rte_solvers = RefCountedPtr<rte_layout> (new rte_layout(plaskin));
  rte_solvers->set_amr(amr);
  rte_solvers->set_computational_geometry(compgeom);
  rte_solvers->set_physical_domain(physdom);
  rte_solvers->sanity_check();
  rte_solvers->allocate_internals();
  rte_solvers->set_source(0.0);
  rte_solvers->write_plot_file();

  // Iteration loop. Shut up during advance. 
  int max_steps  = 100;
  const Real cfl = 0.8;

  amr->set_verbosity(-1);
  cdr_solvers->set_verbosity(-1);
  rte_solvers->set_verbosity(-1);
  for (int i = 0; i < max_steps; i++){
    const Real dt_cfl = cfl*cdr_solvers->compute_cfl_dt();
    const Real dt_dif = cfl*cdr_solvers->compute_diffusive_dt();
    const Real dt = Min(dt_cfl, dt_dif);

    // Advance cdr equations
    cdr_solvers->advance(dt);

    // Set source term for RTE equations and solve
    Vector<RefCountedPtr<cdr_solver> > solvers = cdr_solvers->get_solvers();
    rte_solvers->set_source(solvers[0]->get_state());
    rte_solvers->advance(dt);

    if((i+1) % 20 == 0){
      cdr_solvers->write_plot_file();
      rte_solvers->write_plot_file();
    }

    pout() << "step = " << i << "\t cfl dt = " << dt << "\t charge = " << cdr_solvers->compute_Q() << endl;
    
  }

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
