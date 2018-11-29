#include "mechanical_shaft.H"
#include "cdr_gdnv.H"
#include "cdr_muscl.H"
#include "ParmParse.H"
#include "data_ops.H"
#include "cdr_sg.H"

// This is the advected species
class adv_spec : public species {
public:
  adv_spec(){
    m_name="adv_spec";
    m_mobile=true;
    m_diffusive=false;
    m_charge=0;
  }
  ~adv_spec(){}
  Real initial_data(const RealVect a_pos, const Real a_time) const {
    return exp(-pow((a_pos.vectorLength()),2)/(1.E-5));
  }
};

int main(int argc, char* argv[]){
  // Init MPI
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  char* input_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, input_file);

  // Set up geometry, physical domain, spatial discretization
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new mechanical_shaft());
  RefCountedPtr<physical_domain> physdom         = RefCountedPtr<physical_domain> (new physical_domain());
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());

  // Build the geometry and create grids
  amr->set_physical_domain(physdom);
  amr->build_domains();
  compgeom->build_geometries(*physdom, amr->get_finest_domain(), amr->get_finest_dx(), amr->get_max_ebis_box_size());
  amr->set_mfis(compgeom->get_mfis());

  // Get irregular tags and grow them
  Vector<IntVectSet> geotags = amr->get_irreg_tags();
  for (int lvl = 0; lvl < geotags.size(); lvl++){
    geotags[lvl].grow(amr->get_irreg_growth());
  }
  amr->regrid(geotags);

  // Create a species which will get advected
  RefCountedPtr<species> spec = RefCountedPtr<species> (new adv_spec());

  // Create a convection-diffusion-reaction solver.
  cdr_solver* cdr = (cdr_solver*) (new cdr_gdnv());
  cdr->set_amr(amr);                            // Give the Poisson solver the amr_mesh object
  cdr->set_computational_geometry(compgeom);    // Tell it about electrodes and dielectrics
  cdr->set_physical_domain(physdom);            // Give it information about the physical domain
  cdr->set_phase(phase::gas);
  cdr->set_verbosity(-1);
  cdr->set_species(spec);
  cdr->allocate_internals();


  // Set velocities, diffusion coefficients, source, and domain bc type. 
  cdr->initial_data();
  cdr->set_velocity(-RealVect(BASISV(2)));
  cdr->set_diffco(1.E-4);
  cdr->set_source(0.0);
  cdr->set_domain_bc(cdr_bc::extrap);


  cdr->write_plot_file();
  for (int k=0; k < 50; k++){
    if(procID() == 0){
      std::cout << "step = " << k << std::endl;
    }
    const Real cfl = cdr->compute_cfl_dt();
    cdr->advance(0.9*cfl);
    cdr->write_plot_file();
  }

  overallMemoryUsage();
#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}
