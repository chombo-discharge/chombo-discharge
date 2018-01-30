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
#include "morrow_lowke.H"

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

  // Set up plasma engine
  RefCountedPtr<plasma_engine> engine = RefCountedPtr<plasma_engine> (new plasma_engine(physdom,
											compgeom,
											plaskin,
											timestepper,
											amr));

  // Set up a multifluid Poisson solver
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

      // Set up an eddington sp1 solver
  RefCountedPtr<photon_group> group = RefCountedPtr<photon_group> (new photon_group("photon", 1.0));
  RefCountedPtr<rte_solver> rte = RefCountedPtr<rte_solver> (new eddington_sp1());
  rte->set_verbosity(10);
  rte->set_amr(amr);
  rte->set_computational_geometry(compgeom);
  rte->set_physical_domain(physdom);
  rte->set_photon_group(group);
  rte->sanity_check();
  
  // Setup plasma engine
  engine->set_verbosity(10);
  engine->set_geom_refinement_depth(1);
  engine->setup_fresh();

  // Poisson solver solves
  poisson->sanity_check();
  poisson->allocate_internals();
  poisson->solve();
  poisson->write_plot_file();

  // RTE solves
  rte->allocate_internals();
  rte->advance(0.0);
  rte->write_plot_file(0);

  // Create a time stepper
  time_stepper* ts = static_cast<time_stepper*> (new rk2());

  // New cdr solver
  RefCountedPtr<species> spec = RefCountedPtr<species> (new species());
  cdr_solver* cdr = static_cast<cdr_solver*> (new cdr_sg());
  cdr->set_species(spec);
  cdr->set_verbosity(10);
  cdr->set_amr(amr);
  cdr->set_computational_geometry(compgeom);
  cdr->set_physical_domain(physdom);
  cdr->sanity_check();
  cdr->allocate_internals();
  cdr->initial_data();

  // Compute the electric field from the poisson solver and set that to be the velocity
  MFAMRCellData E;
  EBAMRCellData E_gas;
  amr->allocate(E, SpaceDim);
  amr->allocate_ptr(E_gas);
  amr->compute_gradient(E, poisson->get_state());
  amr->alias(E_gas, phase::gas, E);
  data_ops::scale(E_gas, -1.0);
  cdr->initial_data();
  cdr->set_velocity(E_gas);
  cdr->set_diffco(0.002);
  cdr->set_source(0.0);
  cdr->set_ebflux(0.0);
  cdr->write_plot_file();

  const Real cfl = 0.8;
  cdr->set_verbosity(1);
  amr->set_verbosity(0);
  poisson->set_verbosity(0);
  const Real init_mass = cdr->compute_mass();
  for (int i = 0; i < 100; i++){
    const Real dt_cfl = cfl*cdr->compute_cfl_dt();
    const Real dt_dif = cfl*cdr->compute_diffusive_dt();

    Real dt;
    //    if(dt_cfl < dt_dif){
      dt = dt_cfl;
      dt = dt_dif;
      dt = Min(dt_cfl, dt_dif);
      pout() << "step = " << i << "\t cfl dt = " << dt << "\t mass = " << cdr->compute_mass()/init_mass << endl;
      cdr->advance(dt);
    // }
    // else{
    //   dt = dt_dif;
    //   pout() << "step = " << i << "\t diff dt = " << dt << "\t mass = " << cdr->compute_mass()/init_mass << endl;
    //   cdr->advance(dt);
    // }

    
    if((i+1) % 10 == 0){
      cdr->write_plot_file();
    }
  }

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
