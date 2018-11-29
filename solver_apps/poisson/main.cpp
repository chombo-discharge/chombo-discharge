#include "plasma_engine.H"
#include "geo_coarsener.H"
#include "mechanical_shaft.H"
#include "poisson_multifluid_gmg.H"
#include "plasma_engine.H"
#include "ParmParse.H"

Real potential_curve(const Real a_time){ 
  return 1.0;
}

int main(int argc, char* argv[]){
  // Init MPI
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  char* input_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, input_file);

  // Set up geometry, physical domain, spatial discretization, and a geo_coarsener
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new mechanical_shaft());
  RefCountedPtr<physical_domain> physdom         = RefCountedPtr<physical_domain> (new physical_domain());
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());

  // Build the geometry and create grids
  amr->set_physical_domain(physdom);
  amr->build_domains();
  compgeom->build_geometries(*physdom, amr->get_finest_domain(), amr->get_finest_dx(), amr->get_max_ebis_box_size());
  amr->set_mfis(compgeom->get_mfis());

  // Get irregular tags and grow them
  Vector<IntVectSet> tags = amr->get_irreg_tags();
  for (int lvl = 0; lvl < tags.size(); lvl++){
    tags[lvl].grow(amr->get_irreg_growth());
  }
  amr->regrid(tags);

  // Create volume source term and surface source term. Set them to zero. 
  MFAMRCellData rhs;
  EBAMRIVData sigma;
  amr->allocate(rhs, 1);
  amr->allocate(sigma, phase::gas, 1);
  data_ops::set_value(rhs, 0.0);
  data_ops::set_value(sigma, 0.0);

  // Set up the Poisson equation
  poisson_solver* poisson = (poisson_solver*) (new poisson_multifluid_gmg()); // Instantiate solver
  poisson->set_verbosity(-1);                       // Print messages
  poisson->set_amr(amr);                            // Give the Poisson solver the amr_mesh object
  poisson->set_computational_geometry(compgeom);    // Tell it about electrodes and dielectrics
  poisson->set_physical_domain(physdom);            // Give it information about the physical domain
  poisson->set_potential(potential_curve);
  
  // Solve the Poisson equation and write output
  poisson->allocate_internals();                    // Allocate internals
  poisson->solve(poisson->get_state(), rhs, sigma, true);
  poisson->write_plot_file();

  overallMemoryUsage();
#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}
