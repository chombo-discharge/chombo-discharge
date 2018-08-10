/*!
  @file   plasma_engine.cpp
  @brief  Implementation of plasma_engine.H
  @author Robert Marskar
  @date   Nov. 2017
*/

#include "plasma_engine.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"
#include "data_ops.H"
#include "mfalias.H"
#include "tags_factory.H"
#include "units.H"

#include <EBAlias.H>
#include <LevelData.H>
#include <EBCellFAB.H>
#include <EBAMRIO.H>
#include <EBAMRDataOps.H>
#include <ParmParse.H>

Real plasma_engine::s_constant_one(const RealVect a_pos){
  return 1.0;
}

plasma_engine::plasma_engine(){
  CH_TIME("plasma_engine::plasma_engine(weak)");
  if(m_verbosity > 5){
    pout() << "plasma_engine::plasma_engine(weak)" << endl;
  }

  MayDay::Abort("plasma_engine::plasma_engine - weak construction is not allowed (yet)");

}

plasma_engine::plasma_engine(const RefCountedPtr<physical_domain>&        a_physdom,
			     const RefCountedPtr<computational_geometry>& a_compgeom,
			     const RefCountedPtr<plasma_kinetics>&        a_plaskin,
			     const RefCountedPtr<time_stepper>&           a_timestepper,
			     const RefCountedPtr<amr_mesh>&               a_amr,
			     const RefCountedPtr<cell_tagger>&            a_celltagger,
			     const RefCountedPtr<geo_coarsener>&          a_geocoarsen){
  CH_TIME("plasma_engine::plasma_engine(full)");

  this->set_verbosity(1);
  if(m_verbosity > 5){
    pout() << "plasma_engine::plasma_engine(full)" << endl;
  }
  

  this->set_physical_domain(a_physdom);                      // Set physical domain
  this->set_computational_geometry(a_compgeom);              // Set computational geometry
  this->set_plasma_kinetics(a_plaskin);                      // Set plasma kinetics
  this->set_time_stepper(a_timestepper);                     // Set time stepper
  this->set_amr(a_amr);                                      // Set amr
  this->set_cell_tagger(a_celltagger);                       // Set cell tagger
  this->set_geo_coarsen(a_geocoarsen);                       // Set geo coarsener
  this->set_poisson_wall_func(s_constant_one);                       // Set wall function

  this->set_geom_refinement_depth(-1);                       // Set geometric refinement depths
  this->set_verbosity(1);                                    // Set verbosity
  this->set_plot_interval(10);                               // Set plot interval
  this->set_checkpoint_interval(10);                         // Set checkpoint interval
  this->set_num_plot_ghost(0);                               // Set number of ghost cells in plots
  this->set_grow_tags(0);                                    // Grow tagged cells
  this->set_regrid_interval(10);                             // Set regrid interval
  this->set_output_directory("./");                          // Set output directory
  this->set_output_file_names("simulation");                 // Set output file names
  this->set_output_mode(output_mode::full);                  // Set output mode
  this->set_restart_mode(restart_mode::full);                // Set restart mode
  this->set_memory_report_mode(memory_report_mode::overall); // Set memory report mode
  this->set_init_regrids(0);                                 // Number of initial regrids
  this->set_geom_only(false);                                // Only plot geometry
  this->set_dump_mass(false);                                // Dump mass to file
  this->set_dump_charge(false);                              // Dump charges to file
  this->set_ebis_memory_load_balance(false);                 // Set load balance for EBIS generation
  this->set_restart(false);                                  // Restart mode
  this->set_write_ebis(false);                               // Write EBIS
  this->set_read_ebis(false);                                // Read EBIS when restarting simulation
  this->set_restart_step(0);                                 // Restart from this step
  this->set_start_time(0.0);                                 // Start time
  this->set_stop_time(1.0);                                  // Stop time
  this->set_max_steps(100);                                  // Max number of steps
  this->set_max_plot_depth(-1);                              // Set maximum plot depth.
  this->set_max_chk_depth(-1);                               // Set maximum checkpoint depth.

  m_amr->set_physical_domain(m_physdom); // Set physical domain
  m_amr->sanity_check();                 // Sanity check, make sure everything is set up correctly
  m_amr->build_domains();                // Build domains and resolutions, nothing else

  if(!m_celltagger.isNull()){ 
    m_celltagger->define(m_plaskin, m_timestepper, m_amr, m_compgeom, m_physdom);
  }

  m_potential_set = false;
  m_step          = 0;
  m_time          = 0.0;
}

plasma_engine::~plasma_engine(){
  CH_TIME("plasma_engine::~plasma_engine");
}

void plasma_engine::add_potential_to_output(EBAMRCellData& a_output, const int a_cur_var){
  CH_TIME("plasma_engine::add_potential_to_output");
  if(m_verbosity > 5){
    pout() << "plasma_engine::add_potential_to_output" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  const RefCountedPtr<EBIndexSpace> ebis_gas   = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol   = m_mfis->get_ebis(phase::solid);
  const RefCountedPtr<poisson_solver>& poisson = m_timestepper->get_poisson();
  const MFAMRCellData& potential               = poisson->get_state();

  EBAMRCellData scratch;
  m_amr->allocate(scratch, phase::gas, ncomp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    LevelData<EBCellFAB> state_gas, state_sol;

    mfalias::aliasMF(state_gas, phase::gas,   *potential[lvl]);

    state_gas.localCopyTo(*scratch[lvl]);
    
    // Copy all covered cells from the other phase
    if(!ebis_sol.isNull()){
      mfalias::aliasMF(state_sol, phase::solid, *potential[lvl]);
      
      for (DataIterator dit = state_sol.dataIterator(); dit.ok(); ++dit){
	const Box box = state_sol.disjointBoxLayout().get(dit());
	const IntVectSet ivs(box);
	const EBISBox& ebisb_gas = state_gas[dit()].getEBISBox();
	const EBISBox& ebisb_sol = state_sol[dit()].getEBISBox();

	FArrayBox& scratch_gas    = (*scratch[lvl])[dit()].getFArrayBox();
	const FArrayBox& data_gas = state_gas[dit()].getFArrayBox();
	const FArrayBox& data_sol = state_sol[dit()].getFArrayBox();

	for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit){
	  const IntVect iv = ivsit();
	  if(ebisb_gas.isCovered(iv) && !ebisb_sol.isCovered(iv)){ // Regular cells from phase 2
	    scratch_gas(iv, comp) = data_sol(iv, comp);
	  }
	  if(ebisb_sol.isIrregular(iv) && ebisb_gas.isIrregular(iv)){ // Irregular cells
	    //	    scratch_gas(iv, comp) = 0.5*(data_gas(iv, comp) + data_sol(iv, comp));
	    scratch_gas(iv, comp) = data_gas(iv,comp);
	  }
	}
      }
    }
  }

  m_amr->average_down(scratch, phase::gas);
  m_amr->interp_ghost(scratch, phase::gas);
  m_amr->interpolate_to_centroids(scratch, phase::gas);

  const Interval src_interv(comp, comp);
  const Interval dst_interv(a_cur_var, a_cur_var + ncomp -1);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
  }

  // Set covered potential
  //  poisson->set_covered_potential(a_output, a_cur_var, m_time);
}

void plasma_engine::add_electric_field_to_output(EBAMRCellData& a_output, const int a_cur_var){
  CH_TIME("plasma_engine::add_electric_field_to_output");
  if(m_verbosity > 10){
    pout() << "plasma_engine::add_electric_field_to_output" << endl;
  }

  const int ncomp        = SpaceDim;
  const int finest_level = m_amr->get_finest_level();

  const RefCountedPtr<EBIndexSpace> ebis_gas   = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol   = m_mfis->get_ebis(phase::solid);
  const RefCountedPtr<poisson_solver>& poisson = m_timestepper->get_poisson();
  const MFAMRCellData& potential               = poisson->get_state();

  MFAMRCellData E;
  EBAMRCellData scratch;
  m_amr->allocate(scratch, phase::gas, ncomp);
  m_amr->allocate(E, ncomp);

  m_timestepper->compute_E(E, potential);
  m_amr->compute_gradient(E, poisson->get_state());
  data_ops::scale(E, -1.0);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    LevelData<EBCellFAB> E_gas, E_sol;

    mfalias::aliasMF(E_gas, phase::gas,*E[lvl]);

    // Copy all covered cells from the other phase
    if(!ebis_sol.isNull()){
      mfalias::aliasMF(E_sol, phase::solid, *E[lvl]);
      
      for (DataIterator dit = E_sol.dataIterator(); dit.ok(); ++dit){
	const Box box = E_sol.disjointBoxLayout().get(dit());
	const IntVectSet ivs(box);
	const EBISBox& ebisb_gas = E_gas[dit()].getEBISBox();
	const EBISBox& ebisb_sol = E_sol[dit()].getEBISBox();

	FArrayBox& data_gas = E_gas[dit()].getFArrayBox();
	FArrayBox& data_sol = E_sol[dit()].getFArrayBox();

	for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit){
	  const IntVect iv = ivsit();
	  if(ebisb_gas.isCovered(iv) && !ebisb_sol.isCovered(iv)){ // Regular cells from phase 2
	    for (int comp = 0; comp < ncomp; comp++){
	      data_gas(iv, comp) = data_sol(iv,comp);
	    }
	  }
	}
      }
    }

    E_gas.localCopyTo(*scratch[lvl]);
  }

  m_amr->average_down(scratch, phase::gas);
  m_amr->interp_ghost(scratch, phase::gas);
  m_amr->interpolate_to_centroids(scratch, phase::gas);

  const Interval src_interv(0, ncomp - 1);
  const Interval dst_interv(a_cur_var, a_cur_var + ncomp -1);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
  }
}

void plasma_engine::add_space_charge_to_output(EBAMRCellData& a_output, const int a_cur_var){
  CH_TIME("plasma_engine::add_space_charge_to_output");
  if(m_verbosity > 10){
    pout() << "plasma_engine::add_space_charge_to_output" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  EBAMRCellData scratch;
  MFAMRCellData rho;
  m_amr->allocate(rho, ncomp);
  m_amr->allocate(scratch, phase::gas, ncomp);

  RefCountedPtr<cdr_layout>& cdr = m_timestepper->get_cdr();
  Vector<EBAMRCellData*> densities; 
  for (cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    densities.push_back(&(solver->get_state()));
  }

  m_timestepper->compute_rho(rho, densities, centering::cell_center);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    LevelData<EBCellFAB> rho_gas;
    mfalias::aliasMF(rho_gas, phase::gas, *rho[lvl]);

    rho_gas.localCopyTo(*scratch[lvl]);
  }

  m_amr->average_down(scratch, phase::gas);
  m_amr->interp_ghost(scratch, phase::gas);
  m_amr->interpolate_to_centroids(scratch, phase::gas);


  const Interval src_interv(comp, comp);
  const Interval dst_interv(a_cur_var, a_cur_var + ncomp -1);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);

    data_ops::set_covered_value(*a_output[lvl], a_cur_var, 0.0);
  }
}

void plasma_engine::add_surface_charge_to_output(EBAMRCellData& a_output, const int a_cur_var){
  CH_TIME("plasma_engine::add_surface_charge_to_output");
  if(m_verbosity > 10){
    pout() << "plasma_engine::add_surface_charge_to_output" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  RefCountedPtr<sigma_solver>& sigma = m_timestepper->get_sigma();
  EBAMRIVData& data = sigma->get_state();
  sigma->reset_cells(data);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    for (DataIterator dit = data[lvl]->dataIterator(); dit.ok(); ++dit){
      BaseIVFAB<Real>& sig   = (*data[lvl])[dit()];
      const EBGraph& ebgraph = sig.getEBGraph();
      const IntVectSet& ivs  = sig.getIVS();

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	(*a_output[lvl])[dit()](vof, a_cur_var) = sig(vof, comp);
      }
    }
  }
}

void plasma_engine::add_current_density_to_output(EBAMRCellData& a_output, const int a_cur_var){
  CH_TIME("plasma_engine::add_current_density_to_output");
  if(m_verbosity > 10){
    pout() << "plasma_engine::add_current_density_to_output" << endl;
  }

  RefCountedPtr<cdr_layout>& cdr = m_timestepper->get_cdr();

  const int comp         = 0;
  const int ncomp        = SpaceDim;
  const int finest_level = m_amr->get_finest_level();

  // Allocate and compute on centroids
  EBAMRCellData J;
  m_amr->allocate(J, cdr->get_phase(), SpaceDim);
  m_timestepper->compute_J(J);
  m_amr->interpolate_to_centroids(J, cdr->get_phase());


  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Interval src_interv(0, ncomp - 1);
    const Interval dst_interv(a_cur_var, a_cur_var + ncomp -1);
    J[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
  }
}

void plasma_engine::add_cdr_densities_to_output(EBAMRCellData& a_output, const int a_cur_var){
  CH_TIME("plasma_engine::add_cdr_densities_to_output");
  if(m_verbosity > 10){
    pout() << "plasma_engine::add_cdr_densities_to_output" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  RefCountedPtr<cdr_layout>& cdr = m_timestepper->get_cdr();

  EBAMRCellData scratch;
  m_amr->allocate(scratch, cdr->get_phase(), ncomp);

  for (cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    const EBAMRCellData& state        = solver->get_state();
    const int num                     = solver_it.get_solver();

    for (int lvl = 0; lvl <= finest_level; lvl++){
      state[lvl]->localCopyTo(*scratch[lvl]);
    }

    m_amr->average_down(scratch, cdr->get_phase());
    m_amr->interp_ghost(scratch, cdr->get_phase());
    m_amr->interpolate_to_centroids(scratch, cdr->get_phase());

    data_ops::floor(scratch, 0.0);

    const Interval src_interv(comp, comp);
    const Interval dst_interv(a_cur_var + num, a_cur_var + num + ncomp -1);
    for (int lvl = 0; lvl <= finest_level; lvl++){
      scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);

      data_ops::set_covered_value(*a_output[lvl], a_cur_var + num, 0.0);
    }
  }
}

void plasma_engine::add_cdr_velocities_to_output(EBAMRCellData& a_output, const int a_cur_var){
  CH_TIME("plasma_engine::add_cdr_velocities_to_output");
  if(m_verbosity > 10){
    pout() << "plasma_engine::add_cdr_velocities_to_output" << endl;
  }

  m_timestepper->compute_cdr_velocities();

  const int comp         = 0;
  const int ncomp        = SpaceDim;
  const int finest_level = m_amr->get_finest_level();

  RefCountedPtr<cdr_layout>& cdr     = m_timestepper->get_cdr();

  EBAMRCellData scratch;
  m_amr->allocate(scratch, cdr->get_phase(), ncomp);

  for (cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    const EBAMRCellData& velo         = solver->get_velo_cell();
    const int num                     = solver_it.get_solver();

    for (int lvl = 0; lvl <= finest_level; lvl++){
      velo[lvl]->localCopyTo(*scratch[lvl]);
    }

    m_amr->average_down(scratch, cdr->get_phase());
    m_amr->interp_ghost(scratch, cdr->get_phase());
    m_amr->interpolate_to_centroids(scratch, cdr->get_phase());

    const Interval src_interv(0, ncomp -1);
    const Interval dst_interv(a_cur_var + num*ncomp, a_cur_var + num*ncomp + ncomp -1);
    for (int lvl = 0; lvl <= finest_level; lvl++){
      scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);

      for (int comp = 0; comp < SpaceDim; comp ++){
	data_ops::set_covered_value(*a_output[lvl], a_cur_var + num*ncomp + comp, 0.0);
      }
    }
  }
}

void plasma_engine::add_cdr_source_to_output(EBAMRCellData& a_output, const int a_cur_var){
  CH_TIME("plasma_engine::add_cdr_source_to_output");
  if(m_verbosity > 10){
    pout() << "plasma_engine::add_cdr_source_to_output" << endl;
  }

  m_timestepper->compute_cdr_sources();
  
  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  RefCountedPtr<cdr_layout>& cdr     = m_timestepper->get_cdr();

  EBAMRCellData scratch;
  m_amr->allocate(scratch, cdr->get_phase(), ncomp);

  for (cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    const EBAMRCellData& state        = solver->get_source();
    const int num                     = solver_it.get_solver();

    for (int lvl = 0; lvl <= finest_level; lvl++){
      state[lvl]->localCopyTo(*scratch[lvl]);
    }

    m_amr->average_down(scratch, cdr->get_phase());
    m_amr->interp_ghost(scratch, cdr->get_phase());
    m_amr->interpolate_to_centroids(scratch, cdr->get_phase());

    const Interval src_interv(comp, comp);
    const Interval dst_interv(a_cur_var + num, a_cur_var + num + ncomp -1);
    for (int lvl = 0; lvl <= finest_level; lvl++){
      scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);

      data_ops::set_covered_value(*a_output[lvl], a_cur_var + num*ncomp, 0.0);
    }
  }
}

void plasma_engine::add_rte_densities_to_output(EBAMRCellData& a_output, const int a_cur_var){
  CH_TIME("plasma_engine::add_rte_densities_to_output");
  if(m_verbosity > 10){
    pout() << "plasma_engine::add_rte_densities_to_output" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  RefCountedPtr<rte_layout>& rte     = m_timestepper->get_rte();
  const phase::which_phase rte_phase = rte->get_phase();

  EBAMRCellData scratch;
  m_amr->allocate(scratch, rte_phase, 1);
  data_ops::set_value(scratch, 0.0);

  for (rte_iterator solver_it(*rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    const EBAMRCellData& state        = solver->get_state();
    const int num                     = solver_it.get_solver();

    for (int lvl = 0; lvl <= finest_level; lvl++){
      state[lvl]->localCopyTo(*scratch[lvl]);
    }

    m_amr->average_down(scratch, rte_phase);
    m_amr->interp_ghost(scratch, rte_phase);
    m_amr->interpolate_to_centroids(scratch, rte_phase);

    const Interval src_interv(comp, comp);
    const Interval dst_interv(a_cur_var + num, a_cur_var + num + ncomp -1);
    for (int lvl = 0; lvl <= finest_level; lvl++){
      scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);

      data_ops::set_covered_value(*a_output[lvl], a_cur_var + num*ncomp, 0.0);
    }
  }
}

void plasma_engine::add_rte_source_to_output(EBAMRCellData& a_output, const int a_cur_var){
  CH_TIME("plasma_engine::add_rte_source_to_output");
  if(m_verbosity > 10){
    pout() << "plasma_engine::add_rte_source_to_output" << endl;
  }
  
  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  RefCountedPtr<rte_layout>& rte     = m_timestepper->get_rte();
  const phase::which_phase rte_phase = rte->get_phase();

  EBAMRCellData scratch;
  m_amr->allocate(scratch, rte_phase, 1);
  data_ops::set_value(scratch, 0.0);

  for (rte_iterator solver_it(*rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    const EBAMRCellData& state        = solver->get_source();
    const int num                     = solver_it.get_solver();

    for (int lvl = 0; lvl <= finest_level; lvl++){
      state[lvl]->localCopyTo(*scratch[lvl]);
    }

    m_amr->average_down(scratch, rte_phase);
    m_amr->interp_ghost(scratch, rte_phase);
    m_amr->interpolate_to_centroids(scratch, rte_phase);

    const Interval src_interv(comp, comp);
    const Interval dst_interv(a_cur_var + num, a_cur_var + num + ncomp -1);
    for (int lvl = 0; lvl <= finest_level; lvl++){
      scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);

      data_ops::set_covered_value(*a_output[lvl], a_cur_var + num*ncomp, 0.0);
    }
  }
}

void plasma_engine::add_tracer_fields_to_output(EBAMRCellData& a_output, const int a_cur_var){
  CH_TIME("plasma_engine::add_tracer_fields_to_output");
  if(m_verbosity > 10){
    pout() << "plasma_engine::add_tracer_fields_to_output" << endl;
  }

  // Not implemented
}

void plasma_engine::allocate_internals(){
  CH_TIME("plasma_engine::allocate_internals");
  if(m_verbosity > 5){
    pout() << "plasma_engine::allocate_internals" << endl;
  }
  
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();
  const IntVect ghost    = IntVect::Zero;

  m_tags.resize(1 + finest_level);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];

    tags_factory fact = tags_factory();
    m_tags[lvl] = RefCountedPtr<LevelData<tags> > (new LevelData<tags>(dbl, ncomp, ghost, fact));
  }
}

void plasma_engine::cache_tags(const EBAMRTags& a_tags){
  CH_TIME("plasma_engine::cache_tags");
  if(m_verbosity > 5){
    pout() << "plasma_engine::cache_tags" << endl;
  }

  const int ncomp         = 1;
  const int finest_level  = m_amr->get_finest_level();
  const IntVect ghost     = IntVect::Zero;

  m_cached_tags.resize(1 + finest_level);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];

    tags_factory fact = tags_factory();
    m_cached_tags[lvl] = RefCountedPtr<LevelData<tags> > (new LevelData<tags>(dbl, ncomp, ghost, fact));

    a_tags[lvl]->localCopyTo(*m_cached_tags[lvl]);
  }
}

void plasma_engine::deallocate_internals(){
  CH_TIME("plasma_engine::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "plasma_engine::deallocate_internals" << endl;
  }

  //  m_amr->deallocate(m_tags);
}

void plasma_engine::write_ebis(){
  CH_TIME("plasma_engine::write_ebis");
  if(m_verbosity > 5){
    pout() << "plasma_engine::write_ebis" << endl;
  }

  const std::string path_gas = m_output_dir + "/gas_ebis.hdf5";
  const std::string path_sol = m_output_dir + "/sol_ebis.hdf5";
  
  HDF5Handle gas_handle(path_gas.c_str(), HDF5Handle::CREATE);
  HDF5Handle sol_handle(path_sol.c_str(), HDF5Handle::CREATE);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  if(!ebis_gas.isNull()){
    ebis_gas->writeAllLevels(gas_handle);
  }
  if(!ebis_sol.isNull()){
    ebis_sol->writeAllLevels(sol_handle);
  }
}

void plasma_engine::get_geom_tags(){
  CH_TIME("plasma_engine::get_geom_tags");
  if(m_verbosity > 5){
    pout() << "plasma_engine::get_geom_tags" << endl;
  }

  const int maxdepth = m_amr->get_max_amr_depth();

  m_geom_tags.resize(maxdepth);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  CH_assert(ebis_gas != NULL);

  for (int lvl = 0; lvl < maxdepth; lvl++){ // Don't need tags on maxdepth, we will never generate grids below that.
    const ProblemDomain& cur_dom = m_amr->get_domains()[lvl];
    const int which_level = ebis_gas->getLevel(cur_dom);

    IntVectSet cond_tags;
    IntVectSet diel_tags;
    IntVectSet gas_tags;
    IntVectSet solid_tags;
    IntVectSet gas_diel_tags;
    IntVectSet gas_solid_tags;

    // Conductor cells
    if(m_conductor_tag_depth > lvl){ 
      cond_tags = ebis_gas->irregCells(which_level);
      if(!ebis_sol.isNull()){
	cond_tags |= ebis_sol->irregCells(which_level);
	cond_tags -= m_mfis->interface_region(cur_dom);
      }
    }

    // Dielectric cells
    if(m_dielectric_tag_depth > lvl){ 
      if(!ebis_sol.isNull()){
	diel_tags = ebis_sol->irregCells(which_level);
      }
    }

    // Gas-solid interface cells
    if(m_gas_solid_interface_tag_depth > lvl){ 
      if(!ebis_sol.isNull()){
	gas_tags = ebis_gas->irregCells(which_level);
      }
    }

    // Gas-dielectric interface cells
    if(m_gas_dielectric_interface_tag_depth > lvl){
      if(!ebis_sol.isNull()){
	gas_diel_tags = m_mfis->interface_region(cur_dom);
      }
    }

    // Gas-conductor interface cells
    if(m_gas_conductor_interface_tag_depth > lvl){ 
      gas_solid_tags = ebis_gas->irregCells(which_level);
      if(!ebis_sol.isNull()){
	gas_solid_tags -= m_mfis->interface_region(cur_dom);
      }
    }

    // Solid-solid interfaces
    if(m_solid_solid_interface_tag_depth > lvl){ 
      if(!ebis_sol.isNull()){
	solid_tags = ebis_sol->irregCells(which_level);

	// Do the intersection with the conductor cells
	IntVectSet tmp = ebis_gas->irregCells(which_level);
	tmp |= ebis_sol->irregCells(which_level);
	tmp -= m_mfis->interface_region(cur_dom);

	solid_tags &= tmp;
      }
    }

    m_geom_tags[lvl].makeEmpty();
    m_geom_tags[lvl] |= diel_tags;
    m_geom_tags[lvl] |= cond_tags;
    m_geom_tags[lvl] |= gas_diel_tags;
    m_geom_tags[lvl] |= gas_solid_tags;
    m_geom_tags[lvl] |= gas_tags;
    m_geom_tags[lvl] |= solid_tags;
  }

  // Remove tags using the geocoarsener if we have it
  if(!m_geocoarsen.isNull()){
    const Vector<real_box> coarsen_boxes = m_geocoarsen->get_coarsen_boxes();
    const Vector<int>      coarsen_level = m_geocoarsen->get_coarsen_levels();


    if(!(coarsen_boxes.size() == coarsen_level.size())){
      if(m_verbosity > 2){
	pout() << "plasma_engine::get_geom_tags - m_geocoarsen is not well defined. Skipping the coarsening step" << endl;
      }
    }
    else{
      if(coarsen_boxes.size() > 0){

	const RealVect origin = m_physdom->get_prob_lo();
	for (int lvl = 0; lvl < maxdepth; lvl++){
	  const Real dx         = m_amr->get_dx()[lvl];
	  const int num_coarsen = coarsen_boxes.size();

	  const IntVectSet tmp = m_geom_tags[lvl];

	  for (IVSIterator it(tmp); it.ok(); ++it){
	    const IntVect iv   = it();
	    const RealVect pos = origin + RealVect(iv)*dx;
	    for (int ibox = 0; ibox < num_coarsen; ibox++){
	      const RealVect lo = coarsen_boxes[ibox].get_lo();
	      const RealVect hi = coarsen_boxes[ibox].get_hi();


	      if(pos > lo && pos < hi && lvl >= coarsen_level[ibox]){
		m_geom_tags[lvl] -= iv;
	      }
	    }
	  }
	}
      }
    }
  }


  // Grow tags. This is an ad-hoc fix that prevents ugly grid near EBs (i.e. cases where only ghost cells are used
  // for elliptic equations)
  const int growth = Max(1, m_amr->get_irreg_growth());
  for (int lvl = 0; lvl < maxdepth; lvl++){
    m_geom_tags[lvl].grow(growth);
  }
}

void plasma_engine::get_loads_and_boxes(long long& a_myPoints,
					long long& a_myPointsGhosts,
					long long& a_myBoxes,
					long long& a_totalPoints,
					long long& a_totalPointsGhosts,
					long long& a_totalBoxes,
					Vector<long long>& a_my_level_boxes,
					Vector<long long>& a_total_level_boxes,
					const int& a_finestLevel,
					const Vector<DisjointBoxLayout>& a_grids){
  CH_TIME("plasma_engine::get_loads_and_boxes");
  if(m_verbosity > 5){
    pout() << "plasma_engine::get_loads_and_boxes" << endl;
  }

  a_myPoints          = 0;
  a_myPointsGhosts    = 0;
  a_myBoxes           = 0;
  a_totalPoints       = 0;
  a_totalPointsGhosts = 0;
  a_totalBoxes        = 0;

  a_my_level_boxes.resize(1 + a_finestLevel);
  a_total_level_boxes.resize(1 + a_finestLevel);

  const int ghost = m_amr->get_num_ghost();

  for (int lvl = 0; lvl <= a_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = a_grids[lvl];
    const Vector<Box> boxes      = dbl.boxArray();
    const Vector<int> procs      = dbl.procIDs();
    
    // Find the total number of points and boxes for this level
    long long pointsThisLevel       = 0;
    long long pointsThisLevelGhosts = 0;
    long long boxesThisLevel        = 0;
    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit){
      Box box      = dbl[lit()];
      Box grownBox = dbl[lit()];
      grownBox.grow(ghost);
      
      //
      pointsThisLevel       += box.numPts();
      pointsThisLevelGhosts += grownBox.numPts();
      boxesThisLevel        += 1;
    }
    a_total_level_boxes[lvl] = boxesThisLevel;
    

    // Find the total number of points and boxes that this processor owns
    long long myPointsLevel       = 0;
    long long myPointsLevelGhosts = 0;
    long long myBoxesLevel        = 0;
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      Box box      = dbl[dit()];
      Box grownBox = dbl[dit()];
      grownBox.grow(3);
      
      myPointsLevel       += box.numPts();
      myPointsLevelGhosts += grownBox.numPts();
      myBoxesLevel        += 1;
    }
    a_my_level_boxes[lvl] = myBoxesLevel;

    // Total for this level
    a_totalPoints           += pointsThisLevel;
    a_totalPointsGhosts     += pointsThisLevelGhosts;
    a_totalBoxes            += boxesThisLevel;
    a_myPoints              += myPointsLevel;
    a_myPointsGhosts        += myPointsLevelGhosts;
    a_myBoxes               += myBoxesLevel;
    a_my_level_boxes[lvl]    = myBoxesLevel;
    a_total_level_boxes[lvl] = boxesThisLevel;
  }
}

void plasma_engine::grid_report(){
  CH_TIME("plasma_engine::grid_report");
  if(m_verbosity > 5){
    pout() << "plasma_engine::grid_report" << endl;
  }

  pout() << endl;

  const int finest_level                 = m_amr->get_finest_level();
  const Vector<DisjointBoxLayout>& grids = m_amr->get_grids();
  const Vector<ProblemDomain>& domains   = m_amr->get_domains();
  const Vector<Real> dx                  = m_amr->get_dx();

  // Grid stuff goes into here
  long long totPoints;
  long long totPointsGhosts;
  long long totBoxes;
  long long myPoints;
  long long myPointsGhosts;
  long long myBoxes;
  Vector<long long> my_level_boxes;
  Vector<long long> total_level_boxes;

  //
  const long long uniformPoints = (domains[finest_level].domainBox()).numPts();

  // Track memory
#ifdef CH_USE_MEMORY_TRACKING
  int BytesPerMB = 1024*1024;
  long long curMem;
  long long peakMem;
  overallMemoryUsage(curMem, peakMem);

#ifdef CH_MPI
  int unfreed_mem = curMem;
  int peak_mem    = peakMem;

  int max_unfreed_mem;
  int max_peak_mem;

  int result1 = MPI_Allreduce(&unfreed_mem, &max_unfreed_mem, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
  int result2 = MPI_Allreduce(&peak_mem,    &max_peak_mem,    1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
#endif
#endif


  
  // Get loads and boxes
  this->get_loads_and_boxes(myPoints,
			    myPointsGhosts,
			    myBoxes,
			    totPoints,
			    totPointsGhosts,
			    totBoxes,
			    my_level_boxes,
			    total_level_boxes,
			    finest_level,
			    grids);

  // Write a report
  pout() << "-----------------------------------------------------------------------" << endl
	 << "plasma_engine::Grid report - timestep = " << m_step << endl 
	 << "\t\t\t        Finest level          = " << finest_level << endl
    	 << "\t\t\t        Total number boxes    = " << totBoxes << endl
    	 << "\t\t\t        Total number of cells = " << totPoints << " (" << totPointsGhosts << ")" << endl
	 << "\t\t\t        Total boxes per level = " << total_level_boxes << endl
    	 << "\t\t\t        Proc. boxes per level = " << my_level_boxes << endl
	 << "\t\t\t        Grid sparsity         = " << 1.0*totPoints/uniformPoints << endl
	 << "\t\t\t        Finest dx             = " << dx[finest_level] << endl
	 << "\t\t\t        Number of cells       = " << myPoints << " (" << myPointsGhosts << ")" << endl
    	 << "\t\t\t        Num boxes             = " << myBoxes << endl
#ifdef CH_USE_MEMORY_TRACKING
	 << "\t\t\t        Unfreed memory        = " << curMem/BytesPerMB << " (MB)" << endl
    	 << "\t\t\t        Peak memory usage     = " << peakMem/BytesPerMB << " (MB)" << endl
#ifdef CH_MPI
    	 << "\t\t\t        Max unfreed memory    = " << max_unfreed_mem/BytesPerMB << " (MB)" << endl
	 << "\t\t\t        Max peak memory       = " << max_peak_mem/BytesPerMB << " (MB)" << endl
#endif
    	 << "-----------------------------------------------------------------------" << endl
#endif
	 << endl;

  pout() << endl;
}

void plasma_engine::memory_report(const memory_report_mode::which_mode a_mode){
#ifdef CH_USE_MEMORY_TRACKING
  CH_TIME("plasma_engine::grid_report");
  if(m_verbosity > 5){
    pout() << "plasma_engine::grid_report" << endl;
  }

  if(a_mode == memory_report_mode::overall){
    overallMemoryUsage();
  }
  else if(a_mode == memory_report_mode::unfreed){
    ReportUnfreedMemory(pout());
  }
  else if(a_mode == memory_report_mode::allocated){
    ReportAllocatedMemory(pout());
  }
  pout() << endl;
#endif
}

void plasma_engine::read_checkpoint_file(const std::string& a_restart_file){
  CH_TIME("plasma_engine::read_checkpoint_file");
  if(m_verbosity > 3){
    pout() << "plasma_engine::read_checkpoint_file" << endl;
  }

  RefCountedPtr<cdr_layout>& cdr         = m_timestepper->get_cdr();
  RefCountedPtr<rte_layout>& rte         = m_timestepper->get_rte();
  RefCountedPtr<poisson_solver>& poisson = m_timestepper->get_poisson();
  RefCountedPtr<sigma_solver>& sigma     = m_timestepper->get_sigma();

  HDF5Handle handle_in(a_restart_file, HDF5Handle::OPEN_RDONLY);
  HDF5HeaderData header;
  header.readFromFile(handle_in);

  // Base resolution should not have changed. If it did, issue an error
  const Real coarsest_dx = header.m_real["coarsest_dx"];
  if(!coarsest_dx == m_amr->get_dx()[0]){
    MayDay::Abort("plasma_engine::read_checkpoint_file - coarsest_dx != dx[0], did you change the base level resolution?!?");
  }

  m_time       = header.m_real["time"];
  m_dt         = header.m_real["dt"];
  m_step       = header.m_int["step"];
  
  int finest_level = header.m_int["finest_level"];

  // Read in grids
  Vector<Vector<Box> > boxes(1 + finest_level);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    handle_in.setGroupToLevel(lvl);

    const int status = read(handle_in, boxes[lvl]);

    if(status != 0){
      MayDay::Error("plasma_engine::read_checkpoint_file - file has no grids");
    }
  }

  m_amr->set_finest_level(finest_level); // Set finest level
  m_amr->set_grids(boxes);               // Set up amr
  m_timestepper->instantiate_solvers();  // Instantiate solvrs, they can now be filled with data
  this->allocate_internals();            // Allocate internal storage which also needs to be filled


  // Transient storage
  EBAMRCellData sig, tags;
  m_amr->allocate(sig, phase::gas, 1);
  m_amr->allocate(tags, phase::gas, 1);
  data_ops::set_value(sig, 0.0);
  data_ops::set_value(tags,  0.0);

  finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    handle_in.setGroupToLevel(lvl);

    // CDR solver
    for (cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver = solver_it();
      EBAMRCellData& solver_state       = solver->get_state();
      const std::string solver_name     = solver->get_name();

      if(m_restart_mode != restart_mode::surface_charge_only){
	read<EBCellFAB>(handle_in, *solver_state[lvl], solver_name, dbl, Interval(), false);
      }
    }

    // RTE solvers
    if(!m_timestepper->stationary_rte()){
      for (rte_iterator solver_it(*rte); solver_it.ok(); ++solver_it){
	RefCountedPtr<rte_solver>& solver = solver_it();
	EBAMRCellData& solver_state       = solver->get_state();
	const std::string solver_name     = solver->get_name();

	read<EBCellFAB>(handle_in, *solver_state[lvl], solver_name, dbl, Interval(), false);
      }
    }

    // Read in sigma and tags
    read<EBCellFAB>(handle_in, *sig[lvl],  "sigma", dbl, Interval(), false);
    read<EBCellFAB>(handle_in, *tags[lvl], "tags",  dbl, Interval(), false);
  }

  // Copy data to sigma solver
  data_ops::set_value(sigma->get_state(), 0.0);
  if(m_restart_mode != restart_mode::volume_charge_only){
    data_ops::incr(sigma->get_state(), sig, 1.0);
  }
  sigma->reset_cells(sigma->get_state());

  // Instantiate m_tags
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBCellFAB& tmp   = (*tags[lvl])[dit()];
      const EBISBox& ebisbox = tmp.getEBISBox();
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

      DenseIntVectSet& tagged_cells = (*m_tags[lvl])[dit()].get_ivs();

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	if(tmp(vof, 0) >= 0.9999){
	  tagged_cells |= vof.gridIndex();
	}
      }
    }
  }


  handle_in.close();
}

void plasma_engine::regrid(const bool a_use_initial_data){
  CH_TIME("plasma_engine::regrid");
  if(m_verbosity > 2){
    pout() << "plasma_engine::regrid" << endl;
  }

  // We need to be careful with memory allocations here. Therefore we do:
  // --------------------------------------------------------------------
  // 1.  Tag cells, this calls cell_tagger which allocates and deallocate its own storage so 
  //     there's a peak in memory consumption here. We have to eat this one because we
  //     potentially need all the solver data for tagging, so that data can't be touched.
  //     If we don't get new tags, we exit this routine already here. 
  // 2.  Deallocate internal storage for the time_stepper - this frees up a bunch of memory that
  //     we don't need since we won't advance until after regridding anyways. 
  // 3.  Cache tags, this doubles up on the memory for m_tags but that shouldn't matter.
  // 4.  Free up m_tags for safety because it will be regridded anyways. 
  // 5.  Cache solver states
  // 6.  Deallocate internal storage for solver. This should free up a bunch of memory. 
  // 7.  Regrid amr_mesh - this shouldn't cause any extra memory issues
  // 8.  Regrid plasma_engine - this
  // 9.  Regrid the cell tagger. I'm not explicitly releasing storage from here since it's so small.
  // 10. Solve elliptic equations and fill solvers


  Vector<IntVectSet> tags;

  const Real start_time = MPI_Wtime();   // Timer

  const bool got_new_tags = this->tag_cells(tags, m_tags);         // Tag cells

  if(!got_new_tags){
    if(a_use_initial_data){
      m_timestepper->initial_data();
    }

    if(m_verbosity > 1){
      pout() << "\nplasma_engine::regrid - Didn't find any new cell tags. Skipping the regrid step\n" << endl;
    }
    return;
  }

  m_timestepper->deallocate_internals(); // Deallocate internal storage for the time stepper.
  this->cache_tags(m_tags);              // Cache m_tags because after regrid, ownership will change
  this->deallocate_internals();          // Deallocate internal storage for plasma_engine

  m_timestepper->cache_states();                // Cache solver states
  m_timestepper->deallocate_solver_internals(); // Deallocate solver internals
  
  const Real cell_tags = MPI_Wtime();    // Timer

  // Regrid base
  const int old_finest_level = m_amr->get_finest_level();
  m_amr->regrid(tags, old_finest_level + 1);
  const Real base_regrid = MPI_Wtime(); // Base regrid time

  const int new_finest_level = m_amr->get_finest_level();
  this->regrid_internals(old_finest_level, new_finest_level);        // Regrid internals for plasma_engine
  m_timestepper->regrid_solvers(old_finest_level, new_finest_level); // Regrid solvers
  m_timestepper->regrid_internals();                                 // Regrid internal storage for time_stepper
  m_celltagger->regrid();                                            // Regrid cell tagger

  if(a_use_initial_data){
    m_timestepper->initial_data();
  }

  const Real solver_regrid = MPI_Wtime(); // Timer

  // Solve the elliptic parts
  bool converged = m_timestepper->solve_poisson();
  if(!converged){ // If we don't converge, try new solver settings
    if(m_verbosity > 0){
      pout() << "plasma_engine::regrid - Poisson solver failed to converge. Trying to auto-tune new settings." << endl;
    }
	  
    RefCountedPtr<poisson_solver> poisson = m_timestepper->get_poisson();
    poisson->auto_tune();
    converged = m_timestepper->solve_poisson();

    if(!converged){
      if(m_verbosity > 0){
	pout() << "plasma_engine::regrid - Poisson solver fails to converge" << endl;
      }
    }
  }
  if(m_timestepper->stationary_rte()){     // Solve RTE equations by using data that exists inside solvers
    const Real dummy_dt = 0.0;
    m_timestepper->solve_rte(dummy_dt);    // Argument does not matter, it's a stationary solver.
  }

  const Real elliptic_solve = MPI_Wtime(); // Elliptic solve time

  // Fill solvers with important stuff
  m_timestepper->compute_cdr_velocities();
  m_timestepper->compute_cdr_diffusion(); 
  m_timestepper->compute_cdr_sources();
  m_timestepper->compute_rte_sources();

  const Real solver_filling = MPI_Wtime();


  const Real stop_time = MPI_Wtime();

  if(m_verbosity > 1){
    this->regrid_report(stop_time - start_time,
			cell_tags - start_time,
			base_regrid - cell_tags,
			solver_regrid - base_regrid,
			elliptic_solve - solver_regrid,
			solver_filling - elliptic_solve);
  }
}

void plasma_engine::regrid_internals(const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("plasma_engine::regrid_internals");
  if(m_verbosity > 2){
    pout() << "plasma_engine::regrid_internals" << endl;
  }

  this->allocate_internals();

  // Copy cached tags back over to m_tags
  for (int lvl = 0; lvl <= a_old_finest_level; lvl++){
    m_cached_tags[lvl]->copyTo(*m_tags[lvl]);

#if 1 // For some reason this needs to happen, but I don't really see how the DenseIntVectSet leaks over
      // when I copy from the old grids to the new. Anyways, this fixes it (somehow).
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      (*m_tags[lvl])[dit()].get_ivs() &= dbl.get(dit());
    }
#endif
  }
}

void plasma_engine::regrid_report(const Real a_total_time,
				  const Real a_tag_time,
				  const Real a_base_regrid_time,
				  const Real a_solver_regrid_time,
				  const Real a_elliptic_solve_time,
				  const Real a_solver_filling_time){
  CH_TIME("plasma_engine::regrid_report");
  if(m_verbosity > 5){
    pout() << "plasma_engine::regrid_report" << endl;
  }

  const Real elapsed    = a_total_time;
  const int elapsed_hrs = floor(elapsed/3600);
  const int elapsed_min = floor((elapsed - 3600*elapsed_hrs)/60);
  const int elapsed_sec = floor( elapsed - 3600*elapsed_hrs - 60*elapsed_min);
  const int elapsed_ms  = floor((elapsed - 3600*elapsed_hrs - 60*elapsed_min - elapsed_sec)*1000.);

  char metrics[30];
  sprintf(metrics, "%3.3ih %2.2im %2.2is %3.3ims",
	  elapsed_hrs, 
	  elapsed_min, 
	  elapsed_sec, 
	  elapsed_ms);

  pout() << "-----------------------------------------------------------------------" << endl
	 << "plasma_engine::regrid_report breakdown - Time step #" << m_step << endl
	 << "\t\t\t" << "Total regrid time : " << metrics << endl
	 << "\t\t\t" << "Cell tagging      : " << 100.*(a_tag_time/a_total_time) << "%" << endl
    	 << "\t\t\t" << "Base regrid       : " << 100.*(a_base_regrid_time/a_total_time) << "%" << endl
	 << "\t\t\t" << "Solver regrid     : " << 100.*(a_solver_regrid_time/a_total_time) << "%" << endl
    	 << "\t\t\t" << "Elliptic solve    : " << 100.*(a_elliptic_solve_time/a_total_time) << "%" << endl
	 << "\t\t\t" << "Solver filling    : " << 100.*(a_solver_filling_time/a_total_time) << "%" << endl
	 << "-----------------------------------------------------------------------" << endl;
}

void plasma_engine::run(const Real a_start_time, const Real a_end_time, const int a_max_steps){
  CH_TIME("plasma_engine::run");
  if(m_verbosity > 1){
    pout() << "plasma_engine::run" << endl;
  }

  if(m_verbosity > 0){
    pout() << "=================================" << endl;
    if(!m_restart){
      pout() << "plasma_engine::run -- starting run" << endl;
    }
    else{
      pout() << "plasma_engine::run -- restarting run" << endl;
    }
  }

  if(a_max_steps > 0){
    if(!m_restart){
      m_time = a_start_time;
      m_step = 0;
    }

    m_timestepper->compute_dt(m_dt, m_timecode);
    m_timestepper->synchronize_solver_times(m_step, m_time, m_dt);

    bool last_step     = false;
    bool first_step    = true;
    const Real init_dt = m_dt;

    if(m_verbosity > 0){
      this->grid_report();
    }

    m_wallclock_start = MPI_Wtime();

    // This is actually a debugging
    ofstream mass_dump_file;
    ofstream charge_dump_file;
    if(m_dump_mass){
      this->open_mass_dump_file(mass_dump_file);
    }
    if(m_dump_charge){
      this->open_charge_dump_file(charge_dump_file);
    }

    while(m_time < a_end_time && m_step < a_max_steps){

      const int max_sim_depth = m_amr->get_max_sim_depth();
      const int max_amr_depth = m_amr->get_max_amr_depth();     
      if(m_step%m_regrid_interval == 0 && m_regrid_interval > 0 && max_sim_depth > 0 && max_amr_depth > 0){
	if(!first_step){
	  this->regrid(false);
	  if(m_verbosity > 0){
	    this->grid_report();
	  }
	}
      }

      if(m_dump_mass){
	this->dump_mass(mass_dump_file);
      }
      if(m_dump_charge){
	this->dump_charge(charge_dump_file);
      }


      if(!first_step){
	m_timestepper->compute_dt(m_dt, m_timecode);
      }

      if(first_step){
	first_step = false;
      }

      // Did the time step become too small?
      if(m_dt < 1.0E-5*init_dt){
	m_step++;

#ifdef CH_USE_HDF5
	this->write_plot_file();
	this->write_checkpoint_file();
#endif

	MayDay::Abort("plasma_engine::run - the time step became too small");
      }

      // Last time step can be smaller than m_dt so that we end on a_end_time
      if(m_time + m_dt > a_end_time){
	m_dt = a_end_time - m_time;
	last_step = true;
      }

      // Time stepper advances solutions
      m_wallclock1 = MPI_Wtime();
      const Real actual_dt = m_timestepper->advance(m_dt);
      m_wallclock2 = MPI_Wtime();

      // Synchronize times
      m_dt    = actual_dt;
      m_time += actual_dt;
      m_step += 1;
      m_timestepper->synchronize_solver_times(m_step, m_time, m_dt);

      if(m_verbosity > 0){
	this->step_report(a_start_time, a_end_time, a_max_steps);
      }


#ifdef CH_USE_HDF5
      if(m_step%m_plot_interval == 0 && m_plot_interval > 0 || last_step == true && m_plot_interval > 0){
	if(m_verbosity > 2){
	  pout() << "plasma_engine::run -- Writing plot file" << endl;
	}
	this->write_plot_file();
      }

      // Write checkpoint file
      if(m_step % m_chk_interval == 0 && m_chk_interval > 0 || last_step == true && m_chk_interval > 0){
	if(m_verbosity > 2){
	  pout() << "plasma_engine::run -- Writing checkpoint file" << endl;
	}
	this->write_checkpoint_file();
      }
#endif
    }

    if(m_dump_mass){
      this->close_mass_dump_file(mass_dump_file);
    }
    if(m_dump_charge){
      this->close_charge_dump_file(charge_dump_file);
    }
  }

  if(m_verbosity > 0){
    pout() << "==================================" << endl;
    pout() << "plasma_engine::run -- ending run  " << endl;
    pout() << "==================================" << endl;
  }

}

void plasma_engine::setup_and_run(){
  CH_TIME("plasma_engine::setup_and_run");
  if(m_verbosity > 0){
    pout() << "plasma_engine::setup_and_run" << endl;
  }

  char iter_str[100];
  sprintf(iter_str, ".check%07d.%dd.hdf5", m_restart_step, SpaceDim);
  const std::string restart_file = m_output_dir + "/" + m_output_names + std::string(iter_str);

  this->setup(m_init_regrids, m_restart, restart_file);

  if(!m_geometry_only){
    this->run(m_start_time, m_stop_time, m_max_steps);
  }
}

void plasma_engine::set_verbosity(const int a_verbosity){
  CH_TIME("plasma_engine::set_verbosity");

  
  m_verbosity = a_verbosity;

  ParmParse pp("plasma_engine");
  pp.query("verbosity", m_verbosity);

  if(m_verbosity > 5){
    pout() << "plasma_engine::set_verbosity" << endl;
  }
}

void plasma_engine::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
  CH_TIME("plasma_engine::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_computational_geometry" << endl;
  }
  m_compgeom = a_compgeom;
  m_mfis     = a_compgeom->get_mfis();
}

void plasma_engine::set_plasma_kinetics(const RefCountedPtr<plasma_kinetics>& a_plaskin){
  CH_TIME("plasma_engine::set_plasma_kinetics");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_plasma_kinetics" << endl;
  }
  m_plaskin = a_plaskin;
}

void plasma_engine::set_time_stepper(const RefCountedPtr<time_stepper>& a_timestepper){
  CH_TIME("plasma_engine::set_time_stepper");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_time_stepper" << endl;
  }
  m_timestepper = a_timestepper;
}

void plasma_engine::set_cell_tagger(const RefCountedPtr<cell_tagger>& a_celltagger){
  CH_TIME("plasma_engine::set_cell_tagger");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_cell_tagger" << endl;
  }
  m_celltagger = a_celltagger;
}

void plasma_engine::set_geo_coarsen(const RefCountedPtr<geo_coarsener>& a_geocoarsen){
  CH_TIME("plasma_engine::set_geo_coarsen");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_geo_coarsen" << endl;
  }
  m_geocoarsen = a_geocoarsen;
}

void plasma_engine::set_geom_refinement_depth(const int a_depth){
  CH_TIME("plasma_engine::set_geom_refinement_depth");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_geom_refinement_depth" << endl;
  }

  const int max_depth = m_amr->get_max_amr_depth();
  
  int depth     = (a_depth < 0) ? max_depth : a_depth;
  int depth1    = depth;
  int depth2    = depth;
  int depth3    = depth;
  int depth4    = depth;
  int depth5    = depth;
  int depth6    = depth;

  { // Get parameter from input script
    ParmParse pp("plasma_engine");
    pp.query("refine_geometry", depth);
    depth  = (depth < 0) ? max_depth : depth;
    depth1 = depth;
    depth2 = depth;
    depth3 = depth;
    depth4 = depth;
    depth5 = depth;
    depth6 = depth;
  }

  { // Get fine controls from input script
    ParmParse pp("plasma_engine");
    pp.query("refine_electrodes",               depth1);
    pp.query("refine_dielectrics",              depth2);
    pp.query("refine_electrode_gas_interface",  depth3);
    pp.query("refine_dielectric_gas_interface", depth4);
    pp.query("refine_solid_gas_interface",      depth5);
    pp.query("refine_solid_solid_interface",    depth6);

    depth1 = (depth1 < 0) ? depth : depth1;
    depth2 = (depth2 < 0) ? depth : depth2;
    depth3 = (depth3 < 0) ? depth : depth3;
    depth4 = (depth4 < 0) ? depth : depth4;
    depth5 = (depth5 < 0) ? depth : depth5;
    depth6 = (depth6 < 0) ? depth : depth6;
  }
  
  this->set_geom_refinement_depth(depth1, depth2, depth3, depth4, depth5, depth6);
}

void plasma_engine::set_geom_refinement_depth(const int a_depth1,
					      const int a_depth2,
					      const int a_depth3,
					      const int a_depth4,
					      const int a_depth5,
					      const int a_depth6){
  CH_TIME("plasma_engine::set_geom_refinement_depth(full");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_geom_refinement_depth(full)" << endl;
  }
  
  const int max_depth = m_amr->get_max_amr_depth();

  m_conductor_tag_depth                = Min(a_depth1, max_depth);
  m_dielectric_tag_depth               = Min(a_depth2, max_depth);
  m_gas_conductor_interface_tag_depth  = Min(a_depth3, max_depth);
  m_gas_dielectric_interface_tag_depth = Min(a_depth4, max_depth);
  m_gas_solid_interface_tag_depth      = Min(a_depth5, max_depth);
  m_solid_solid_interface_tag_depth    = Min(a_depth6, max_depth);


  m_geom_tag_depth = 0;
  m_geom_tag_depth = Max(m_geom_tag_depth, a_depth1);
  m_geom_tag_depth = Max(m_geom_tag_depth, a_depth2);
  m_geom_tag_depth = Max(m_geom_tag_depth, a_depth3);
  m_geom_tag_depth = Max(m_geom_tag_depth, a_depth4);
  m_geom_tag_depth = Max(m_geom_tag_depth, a_depth5);
  m_geom_tag_depth = Max(m_geom_tag_depth, a_depth6);
}

void plasma_engine::set_output_mode(const output_mode::which_mode a_mode){
  CH_TIME("plasma_engine::set_output_mode");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_output_mode" << endl;
  }

  m_output_mode = a_mode;

  std::string str;
  ParmParse pp("plasma_engine");
  pp.query("output_mode", str);
  if(str == "full"){
    m_output_mode = output_mode::full;
  }
  else if(str == "medium"){
    m_output_mode = output_mode::medium;
  }
  else if(str == "light"){
    m_output_mode = output_mode::light;
  }
  else if(str == "ultra_light"){
    m_output_mode = output_mode::light; // Ultra_light not yet supported
  }
}

void plasma_engine::set_restart_mode(const restart_mode::which_mode a_mode){
  CH_TIME("plasma_engine::set_restart_mode");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_restart_mode" << endl;
  }

  m_restart_mode = a_mode;

  std::string str;
  ParmParse pp("plasma_engine");
  pp.query("restart_mode", str);
  if(str == "full"){
    m_restart_mode = restart_mode::full;
  }
  else if(str == "surface_only"){
    m_restart_mode = restart_mode::surface_charge_only;
  }
  else if(str == "volume_only"){
    m_restart_mode = restart_mode::volume_charge_only; 
  }
}

void plasma_engine::set_memory_report_mode(const memory_report_mode::which_mode a_memory_mode){
  CH_TIME("plasma_engine::set_memory_report_mode");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_memory_report_mode" << endl;
  }

  m_memory_mode = a_memory_mode;

  std::string str;
  ParmParse pp("plasma_engine");
  pp.query("memory_report_mode", str);
  if(str == "overall"){
    m_memory_mode = memory_report_mode::overall;
  }
  else if(str == "unfreed"){
    m_memory_mode = memory_report_mode::unfreed;
  }
  else if(str == "allocated"){
    m_memory_mode = memory_report_mode::allocated;
  }
}

void plasma_engine::set_output_directory(const std::string a_output_dir){
  CH_TIME("plasma_engine::set_output_directory");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_output_directory" << endl;
  }

  m_output_dir = a_output_dir;

  ParmParse pp("plasma_engine");
  pp.query("output_directory", m_output_dir);
}

void plasma_engine::set_output_file_names(const std::string a_output_names){
  CH_TIME("plasma_engine::set_output_names");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_output_names" << endl;
  }

  m_output_names = a_output_names;

  ParmParse pp("plasma_engine");
  pp.query("output_names", m_output_names);
}

void plasma_engine::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("plasma_engine::set_amr");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_amr" << endl;
  }

  m_amr = a_amr;
  m_amr->set_physical_domain(m_physdom);
  m_amr->set_mfis(m_compgeom->get_mfis());
}

void plasma_engine::set_potential(Real (*a_potential)(const Real a_time)){
  CH_TIME("plasma_engine::set_potential");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_potential" << endl;
  }
  
  m_potential     = a_potential;
  m_potential_set = true;
}

void plasma_engine::set_poisson_wall_func(const int a_dir, const Side::LoHiSide a_side, Real (*a_func)(const RealVect a_pos)){
  CH_TIME("plasma_engine::set_poisson_wall_func(dir, side, func)");
  if(m_verbosity > 4){
    pout() << "plasma_engine::set_poisson_wall_func(dir, side, func)" << endl;
  }

  if(a_dir == 0){
    if(a_side == Side::Lo){
      m_wall_func_x_lo = a_func;
    }
    else if(a_side == Side::Hi){
      m_wall_func_x_hi = a_func;
    }
  }
  else if(a_dir == 1){
    if(a_side == Side::Lo){
      m_wall_func_y_lo = a_func;
    }
    else if(a_side == Side::Hi){
      m_wall_func_y_hi = a_func;
    }
  }
#if CH_SPACEDIM==3
  else if(a_dir == 2){
    if(a_side == Side::Lo){
      m_wall_func_z_lo = a_func;
    }
    else if(a_side == Side::Hi){
      m_wall_func_z_hi = a_func;
    }
  }
#endif
}

void plasma_engine::set_poisson_wall_func(Real (*a_func)(const RealVect a_pos)){
  CH_TIME("plasma_engine::set_poisson_wall_func(func)");
  if(m_verbosity > 4){
    pout() << "plasma_engine::set_poisson_wall_func(func)" << endl;
  }

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      this->set_poisson_wall_func(dir, sit(), a_func);
    }
  }
}

void plasma_engine::setup(const int a_init_regrids, const bool a_restart, const std::string a_restart_file){
  CH_TIME("plasma_engine::setup");
  if(m_verbosity > 5){
    pout() << "plasma_engine::setup" << endl;
  }

  if(m_geometry_only){
    this->setup_geometry_only();
  }
  else{
    if(!a_restart){
      this->setup_fresh(a_init_regrids);
#ifdef CH_USE_HDF5
      if(m_plot_interval > 0){
	this->write_plot_file();
      }
#endif
    }
    else{
      this->setup_for_restart(a_init_regrids, a_restart_file);
    }
  }
}

void plasma_engine::setup_geometry_only(){
  CH_TIME("plasma_engine::setup_geometry_only");
  if(m_verbosity > 5){
    pout() << "plasma_engine::setup_geometry_only" << endl;
  }

  this->sanity_check();

  if(m_ebis_memory_load_balance){
    EBIndexSpace::s_useMemoryLoadBalance = true;
  }
  else {
    EBIndexSpace::s_useMemoryLoadBalance = false;
  }

  m_compgeom->build_geometries(*m_physdom,                 // Build the multifluid geometries
			       m_amr->get_finest_domain(),
			       m_amr->get_finest_dx(),
			       m_amr->get_max_ebis_box_size());

  this->get_geom_tags();       // Get geometric tags.
  
  m_amr->set_num_ghost(m_timestepper->query_ghost()); // Query solvers for ghost cells. Give it to amr_mesh before grid gen.
  
  Vector<IntVectSet> tags = m_geom_tags;
  m_amr->build_grids(tags, m_geom_tag_depth);
  m_amr->define_eblevelgrid();
  //  m_amr->regrid(m_geom_tags, m_geom_tag_depth);       // Regrid using geometric tags for now

  if(m_verbosity > 0){
    this->grid_report();
  }

  this->write_geometry();                             // Write geometry only
}

void plasma_engine::setup_fresh(const int a_init_regrids){
  CH_TIME("plasma_engine::setup_fresh");
  if(m_verbosity > 5){
    pout() << "plasma_engine::setup_fresh" << endl;
  }

  this->sanity_check();                                    // Sanity check before doing anything expensive

  if(m_ebis_memory_load_balance){
    EBIndexSpace::s_useMemoryLoadBalance = true;
  }
  else {
    EBIndexSpace::s_useMemoryLoadBalance = false;
  }

  m_compgeom->build_geometries(*m_physdom,                 // Build the multifluid geometries
			       m_amr->get_finest_domain(),
			       m_amr->get_finest_dx(),
			       m_amr->get_max_ebis_box_size());

  if(m_write_ebis){
    this->write_ebis();        // Write EBIndexSpace's for later use
  }
  this->get_geom_tags();       // Get geometric tags.
  
  m_amr->set_num_ghost(m_timestepper->query_ghost()); // Query solvers for ghost cells. Give it to amr_mesh before grid gen.
  m_amr->regrid(m_geom_tags, m_geom_tag_depth);       // Regrid using geometric tags for now

  this->allocate_internals();

  if(m_verbosity > 0){
    this->grid_report();
  }

  m_timestepper->set_amr(m_amr);
  m_timestepper->set_plasma_kinetics(m_plaskin);
  m_timestepper->set_computational_geometry(m_compgeom);       // Set computational geometry
  m_timestepper->set_physical_domain(m_physdom);               // Physical domain
  m_timestepper->set_potential(m_potential);                   // Potential
  m_timestepper->set_poisson_wall_func(0, Side::Lo, m_wall_func_x_lo); // Set function-based Poisson on xlo
  m_timestepper->set_poisson_wall_func(0, Side::Hi, m_wall_func_x_hi); // Set function-based Poisson on xhi
  m_timestepper->set_poisson_wall_func(1, Side::Lo, m_wall_func_y_lo); // Set function-based Poisson on ylo
  m_timestepper->set_poisson_wall_func(1, Side::Hi, m_wall_func_y_hi); // Set function-based Poisson on yhi
#if CH_SPACEDIM==3
  m_timestepper->set_poisson_wall_func(2, Side::Lo, m_wall_func_z_lo); // Set function-based Poisson on zlo
  m_timestepper->set_poisson_wall_func(2, Side::Hi, m_wall_func_z_hi); // Set function-based Poisson on zhi
#endif
  m_timestepper->instantiate_solvers();                   // Instantiate sigma and cdr with initial data (and rte if transient)
  m_timestepper->initial_data();                          // Fill cdr and rte with initial data

  if (a_init_regrids >= 0){
    RefCountedPtr<poisson_solver> poisson = m_timestepper->get_poisson();
    poisson->auto_tune();
    m_timestepper->solve_poisson();                       // Solve Poisson equation by using initial data
    if(m_timestepper->stationary_rte()){                  // Solve RTE equations by using initial data and electric field
      const Real dummy_dt = 0.0;
      m_timestepper->solve_rte(dummy_dt);                 // Argument does not matter, it's a stationary solver.
    }
    
  }

  if(!m_celltagger.isNull()){
    m_celltagger->regrid();
  }
  m_timestepper->regrid_internals();

  // Fill solvers with important stuff
  m_timestepper->compute_cdr_velocities();
  m_timestepper->compute_cdr_diffusion(); 
  m_timestepper->compute_cdr_sources();
  m_timestepper->compute_rte_sources();

  // Initial regrids
  for (int i = 0; i < a_init_regrids; i++){
    if(m_verbosity > 5){
      pout() << "plasma_engine::initial_regrids" << endl;
    }
    
    this->regrid(true);

    if(m_verbosity > 0){
      this->grid_report();
    }
  }
}

void plasma_engine::setup_for_restart(const int a_init_regrids, const std::string a_restart_file){
  CH_TIME("plasma_engine::setup_for_restart");
  if(m_verbosity > 5){
    pout() << "plasma_engine::setup_for_restart" << endl;
  }

  this->sanity_check();                                    // Sanity check before doing anything expensive

  m_compgeom->build_geometries(*m_physdom,                 // Build the multifluid geometries
			       m_amr->get_finest_domain(),
			       m_amr->get_finest_dx(),
			       m_amr->get_max_ebis_box_size());

  this->get_geom_tags();       // Get geometric tags.

  m_timestepper->set_amr(m_amr);                         // Set amr
  m_timestepper->set_plasma_kinetics(m_plaskin);         // Set plasma kinetics
  m_timestepper->set_computational_geometry(m_compgeom); // Set computational geometry
  m_timestepper->set_physical_domain(m_physdom);         // Physical domain
  m_timestepper->set_potential(m_potential);             // Potential
  m_timestepper->set_poisson_wall_func(0, Side::Lo, m_wall_func_x_lo); // Set function-based Poisson on xlo
  m_timestepper->set_poisson_wall_func(0, Side::Hi, m_wall_func_x_hi); // Set function-based Poisson on xhi
  m_timestepper->set_poisson_wall_func(1, Side::Lo, m_wall_func_y_lo); // Set function-based Poisson on ylo
  m_timestepper->set_poisson_wall_func(1, Side::Hi, m_wall_func_y_hi); // Set function-based Poisson on yhi
#if CH_SPACEDIM==3
  m_timestepper->set_poisson_wall_func(2, Side::Lo, m_wall_func_z_lo); // Set function-based Poisson on zlo
  m_timestepper->set_poisson_wall_func(2, Side::Hi, m_wall_func_z_hi); // Set function-based Poisson on zhi
#endif
  m_amr->set_num_ghost(m_timestepper->query_ghost());    // Query solvers for ghost cells. Give it to amr_mesh before grid gen.

  this->read_checkpoint_file(a_restart_file); // Read checkpoint file - this sets up amr, instantiates solvers and fills them

  if(m_restart_mode == restart_mode::surface_charge_only){
    m_timestepper->initial_cdr_data();
    if(!m_timestepper->stationary_rte()){
      m_timestepper->initial_rte_data();
    }
  }

  m_timestepper->solve_poisson();                       // Solve Poisson equation by 
  if(m_timestepper->stationary_rte()){                  // Solve RTE equations if stationary solvers
    const Real dummy_dt = 0.0;
    m_timestepper->solve_rte(dummy_dt);                 // Argument does not matter, it's a stationary solver.
  }
  m_timestepper->regrid_internals(); // Prepare internal storage for time stepper
  if(!m_celltagger.isNull()){
    m_celltagger->regrid();            // Prepare internal storage for cell tagger
  }

  // Fill solvers with important stuff
  m_timestepper->compute_cdr_velocities();
  m_timestepper->compute_cdr_diffusion(); 
  m_timestepper->compute_cdr_sources();
  m_timestepper->compute_rte_sources();


  // Initial regrids
  for (int i = 0; i < a_init_regrids; i++){
    if(m_verbosity > 0){
      pout() << "plasma_engine -- initial regrid # " << i + 1 << endl;
    }
      
    this->regrid(false);

    if(m_verbosity > 0){
      this->grid_report();
    }
  }
}

void plasma_engine::set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom){
  CH_TIME("plasma_engine::set_physical_domain");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_physical_domain" << endl;
  }
  m_physdom = a_physdom;
}

void plasma_engine::set_regrid_interval(const int a_regrid_interval){
  CH_TIME("plasma_engine::set_regrid_interval");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_regrid_interval" << endl;
  }
  
  m_regrid_interval = a_regrid_interval;

  ParmParse pp("plasma_engine");
  pp.query("regrid_interval", m_regrid_interval);
}

void plasma_engine::set_plot_interval(const int a_plot_interval){
  CH_TIME("plasma_engine::set_plot_interval");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_plot_interval" << endl;
  }
  
  m_plot_interval = a_plot_interval;

  ParmParse pp("plasma_engine");
  pp.query("plot_interval", m_plot_interval);
}

void plasma_engine::set_checkpoint_interval(const int a_chk_interval){
  CH_TIME("plasma_engine::set_checkpoint_interval");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_checkpoint_interval" << endl;
  }
  
  m_chk_interval = a_chk_interval;

  ParmParse pp("plasma_engine");
  pp.query("checkpoint_interval", m_chk_interval);
}

void plasma_engine::set_num_plot_ghost(const int a_num_plot_ghost){
  CH_TIME("plasma_engine::set_num_plot_ghost");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_num_plot_ghost" << endl;
  }
  
  m_num_plot_ghost = a_num_plot_ghost;


  ParmParse pp("plasma_engine");
  pp.query("num_plot_ghost", m_num_plot_ghost);

  m_num_plot_ghost = (m_num_plot_ghost < 0) ? 0 : m_num_plot_ghost;
  m_num_plot_ghost = (m_num_plot_ghost > 3) ? 3 : m_num_plot_ghost;
}

void plasma_engine::set_restart(const bool a_restart){
  CH_TIME("plasma_engine::set_restart");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_restart" << endl;
  }
  
  m_restart = a_restart;

  std::string str = "";
  ParmParse pp("plasma_engine");
  pp.query("restart", str);
  if(str == "true"){
    m_restart = true;
  }
  else if(str == "false"){
    m_restart = false;
  }
}

void plasma_engine::set_ebis_memory_load_balance(const bool a_balance){
  CH_TIME("plasma_engine::set_ebis_memory_load_balance");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_ebis_memory_load_balance" << endl;
  }

  m_ebis_memory_load_balance = a_balance;

  { // get parameter from input script
    std::string str;
    ParmParse pp("plasma_engine");
    pp.query("ebis_memory_load_balance", str);
    if(str == "true"){
      m_ebis_memory_load_balance = true;
    }
    else if(str == "false"){
      m_ebis_memory_load_balance = false;
    }
  }
}

void plasma_engine::set_write_ebis(const bool a_write_ebis){
  CH_TIME("plasma_engine::set_write_ebis");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_write_ebis" << endl;
  }

  m_write_ebis = a_write_ebis;

  { // get parameter from input script
    std::string str;
    ParmParse pp("plasma_engine");
    pp.query("write_ebis", str);
    if(str == "true"){
      m_write_ebis = true;
    }
    else if(str == "false"){
      m_write_ebis = false;
    }
  }
}

void plasma_engine::set_read_ebis(const bool a_read_ebis){
  CH_TIME("plasma_engine::set_read_ebis");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_read_ebis" << endl;
  }

  m_read_ebis = a_read_ebis;

  { // get parameter from input script
    std::string str;
    ParmParse pp("plasma_engine");
    pp.query("read_ebis", str);
    if(str == "true"){
      m_read_ebis = true;
    }
    else if(str == "false"){
      m_read_ebis = false;
    }
  }
}

void plasma_engine::set_dump_mass(const bool a_dump_mass){
  CH_TIME("plasma_engine::set_dump_mass");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_dump_mass" << endl;
  }

  m_dump_mass = a_dump_mass;

  { // get parameter from input script
    std::string str;
    ParmParse pp("plasma_engine");
    pp.query("dump_mass", str);
    if(str == "true"){
      m_dump_mass = true;
    }
    else if(str == "false"){
      m_dump_mass = false;
    }
  }
}

void plasma_engine::set_dump_charge(const bool a_dump_charge){
  CH_TIME("plasma_engine::set_dump_charge");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_dump_charge" << endl;
  }

  m_dump_charge = a_dump_charge;

  { // get parameter from input script
    std::string str;
    ParmParse pp("plasma_engine");
    pp.query("dump_charge", str);
    if(str == "true"){
      m_dump_charge = true;
    }
    else if(str == "false"){
      m_dump_charge = false;
    }
  }
}

void plasma_engine::set_geom_only(const bool a_geom_only){
  CH_TIME("plasma_engine::set_geom_only");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_geom_only" << endl;
  }

  m_geometry_only = a_geom_only;

  { // Get parameter from input file
    std::string str;
    ParmParse pp("plasma_engine");
    pp.query("geometry_only", str);
    if(str == "true"){
      m_geometry_only = true;
    }
    else if(str == "false"){
      m_geometry_only = false;
    }
  }
}

void plasma_engine::set_restart_step(const int a_restart_step){
  CH_TIME("plasma_engine::set_restart_step");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_restart_step" << endl;
  }
  
  m_restart_step = a_restart_step;

  int step;
  ParmParse pp("plasma_engine");
  pp.query("restart_step", step);
  if(step >= 0){
    m_restart_step = step;
  }
}

void plasma_engine::set_start_time(const Real a_start_time){
  CH_TIME("plasma_engine::set_start_time");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_start_time" << endl;
  }
  
  m_start_time = a_start_time;

  ParmParse pp("plasma_engine");
  pp.query("start_time", m_start_time);
}

void plasma_engine::set_stop_time(const Real a_stop_time){
  CH_TIME("plasma_engine::set_stop_time");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_stop_time" << endl;
  }
  
  m_stop_time = a_stop_time;
  
  ParmParse pp("plasma_engine");
  pp.query("stop_time", m_stop_time);
}

void plasma_engine::set_max_steps(const int a_max_steps){
  CH_TIME("plasma_engine::set_max_steps");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_max_steps" << endl;
  }
  
  m_max_steps = a_max_steps;

  ParmParse pp("plasma_engine");
  pp.query("max_steps", m_max_steps);
}

void plasma_engine::set_max_plot_depth(const int a_max_plot_depth){
  CH_TIME("plasma_engine::set_max_plot_depth");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_max_plot_depth" << endl;
  }

  m_max_plot_depth = a_max_plot_depth;

  {
    ParmParse pp("plasma_engine");
    pp.query("max_plot_depth", m_max_plot_depth);
  }
}

void plasma_engine::set_max_chk_depth(const int a_max_chk_depth){
  CH_TIME("plasma_engine::set_max_chk_depth");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_max_chk_deepth" << endl;
  }

  m_max_chk_depth = a_max_chk_depth;

  {
    ParmParse pp("plasma_engine");
    pp.query("max_chk_depth", m_max_chk_depth);
  }
}

void plasma_engine::set_init_regrids(const int a_init_regrids){
  CH_TIME("plasma_engine::set_init_regrids");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_init_regrids" << endl;
  }
  
  m_init_regrids = a_init_regrids;


  ParmParse pp("plasma_engine");
  pp.query("initial_regrids", m_init_regrids);

}

void plasma_engine::set_grow_tags(const int a_grow_tags){
  CH_TIME("plasma_engine::set_grow_tags");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_grow_tags" << endl;
  }
  
  m_grow_tags = a_grow_tags;

  ParmParse pp("plasma_engine");
  pp.query("grow_tags", m_grow_tags);

  m_grow_tags = Max(0, m_grow_tags);
}

void plasma_engine::sanity_check(){
  CH_TIME("plasma_engine::sanity_check");
  if(m_verbosity > 4){
    pout() << "plasma_engine::sanity_check" << endl;
  }

  CH_assert(!m_timestepper.isNull());
  CH_assert(m_potential_set);
}

void plasma_engine::step_report(const Real a_start_time, const Real a_end_time, const int a_max_steps){
  CH_TIME("plasma_engine::step_report");
  if(m_verbosity > 5){
    pout() << "plasma_engine::step_report" << endl;
  }

  // Compute the maximum electric field
  Real Emax;
  m_timestepper->compute_Emax(Emax, phase::gas);

  //
  Real nmax;
  std::string solver_max;
  m_timestepper->get_cdr_max(nmax, solver_max);

  pout() << endl;
  std::string str;
  if(m_timecode == time_code::cfl){
    str = " (Restricted by CFL)";
  }
  if(m_timecode == time_code::diffusion){
    str = " (Restricted by diffusion)";
  }
  if(m_timecode == time_code::source){
    str = " (Restricted by source term)";
  }
  if(m_timecode == time_code::relaxation_time){
    str = " (Restricted by relaxation time)";
  }
  if(m_timecode == time_code::restricted){
    str = " (Restricted by time stepper)";
  }
  if(m_timecode == time_code::hardcap){
    str = " (Restricted by a hardcap)";
  }
  pout() << "plasma_engine::Time step report -- Time step #" << m_step << endl
	 << "                                   Time  = " << m_time << endl
	 << "                                   dt    = " << m_dt << str << endl
	 << "                                   Emax  = " << Emax << endl
	 << "                                   n_max = " << nmax << "(" + solver_max + ")" << endl;

  // Get the total number of poitns across all levels
  const int finest_level                 = m_amr->get_finest_level();
  const Vector<DisjointBoxLayout>& grids = m_amr->get_grids();
  const Vector<ProblemDomain>& domains   = m_amr->get_domains();
  const Vector<Real>& dx                 = m_amr->get_dx();
  long long totalPoints = 0;
  long long uniformPoints = (domains[finest_level].domainBox()).numPts();
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    long long pointsThisLevel = 0;
    for (LayoutIterator lit = grids[lvl].layoutIterator(); lit.ok(); ++lit){
      pointsThisLevel += grids[lvl][lit()].numPts();
    }
    totalPoints += pointsThisLevel;
  }

  char metrics[300];
  pout() << endl;

  // Percentage completed of time steps
  const Real percentStep = (1.0*m_step/a_max_steps)*100.;
  sprintf(metrics,"%31c -- %5.2f percentage of time steps completed",' ', percentStep);
  pout() << metrics << endl;

  const Real percentTime = ((m_time - a_start_time)/(a_end_time - a_start_time))*100.;
  sprintf(metrics,"%31c -- %5.2f percentage of simulation time completed",' ', percentTime);
  pout() << metrics << endl;


  // Hours, minutes, seconds and millisecond of the previous iteration
  const Real elapsed   = m_wallclock2 - m_wallclock_start;
  const int elapsedHrs = floor(elapsed/3600);
  const int elapsedMin = floor((elapsed - 3600*elapsedHrs)/60);
  const int elapsedSec = floor( elapsed - 3600*elapsedHrs - 60*elapsedMin);
  const int elapsedMs  = floor((elapsed - 3600*elapsedHrs - 60*elapsedMin - elapsedSec)*1000);

  // Write a string with total elapsed time
  sprintf(metrics, 
	  "%31c -- Elapsed time          : %3.3ih %2.2im %2.2is %3.3ims",
	  ' ',
	  elapsedHrs, 
	  elapsedMin, 
	  elapsedSec, 
	  elapsedMs);
  pout() << metrics << endl;

  // Hours, minutes, seconds and millisecond of the previous iteration
  const Real lastadv = m_wallclock2 - m_wallclock1;
  const int advHrs = floor(lastadv/3600);
  const int advMin = floor((lastadv - 3600*advHrs)/60);
  const int advSec = floor( lastadv - 3600*advHrs - 60*advMin);
  const int advMs  = floor((lastadv - 3600*advHrs - 60*advMin - advSec)*1000);

  // Write a string with the previous iteration metrics
  sprintf(metrics, 
	  "%31c -- Last time step        : %3.3ih %2.2im %2.2is %3.3ims",
	  ' ',
	  advHrs, 
	  advMin, 
	  advSec, 
	  advMs);
  pout() << metrics << endl;

  // This is the time remaining
  const Real maxPercent = Max(percentTime, percentStep);
  const Real remaining  = 100.*elapsed/maxPercent - elapsed;
  const int remHrs = floor(remaining/3600);
  const int remMin = floor((remaining - 3600*remHrs)/60);
  const int remSec = floor( remaining - 3600*remHrs - 60*remMin);
  const int remMs  = floor((remaining - 3600*remHrs - 60*remMin - remSec)*1000);

  // Write a string with the previous iteration metrics
  sprintf(metrics, 
	  "%31c -- Estimated remaining   : %3.3ih %2.2im %2.2is %3.3ims",
	  ' ',
	  remHrs, 
	  remMin, 
	  remSec, 
	  remMs);
  pout() << metrics << endl;

  // Write memory usage
#ifdef CH_USE_MEMORY_TRACKING
  const int BytesPerMB = 1024*1024;
  long long curMem;
  long long peakMem;
  overallMemoryUsage(curMem, peakMem);

  pout() << "                                -- Unfreed memory        : " << curMem/BytesPerMB << "(MB)" << endl;
  pout() << "                                -- Peak memory usage     : " << peakMem/BytesPerMB << "(MB)" << endl;

#ifdef CH_MPI
  int unfreed_mem = curMem;
  int peak_mem    = peakMem;

  int max_unfreed_mem;
  int max_peak_mem;

  int result1 = MPI_Allreduce(&unfreed_mem, &max_unfreed_mem, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
  int result2 = MPI_Allreduce(&peak_mem,    &max_peak_mem,    1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
  pout() << "                                -- Max unfreed memory    : " << max_unfreed_mem/BytesPerMB << "(MB)" << endl;
  pout() << "                                -- Max peak memory usage : " << max_peak_mem/BytesPerMB << "(MB)" << endl;
#endif
#endif

}

bool plasma_engine::tag_cells(Vector<IntVectSet>& a_all_tags, EBAMRTags& a_cell_tags){
  CH_TIME("plasma_engine::tag_cells");
  if(m_verbosity > 5){
    pout() << "plasma_engine::tag_cells" << endl;
  }

  bool got_new_tags = false;

  const int finest_level  = m_amr->get_finest_level();
  a_all_tags.resize(1 + finest_level, IntVectSet());

  if(!m_celltagger.isNull()){
    got_new_tags = m_celltagger->tag_cells(a_cell_tags);
  }

  // Gather tags from a_tags
  for (int lvl = 0; lvl <= finest_level; lvl++){
    for (DataIterator dit = a_cell_tags[lvl]->dataIterator(); dit.ok(); ++dit){
      a_all_tags[lvl] |= IntVectSet((*a_cell_tags[lvl])[dit()].get_ivs());
    }
  }

  // Add geometric tags. 
  for (int lvl = 0; lvl < finest_level; lvl++){
    a_all_tags[lvl].grow(m_grow_tags);
    a_all_tags[lvl] |= m_geom_tags[lvl];
  }

#if 0 // Debug - if this fails, you have tags on m_amr->m_max_amr_depth and something has gone wrong. 
  if(finest_level == m_amr->get_max_amr_depth()){
    for (int lvl = 0; lvl <= finest_level; lvl++){
      pout() << "level = " << lvl << "\t num_pts = " << a_all_tags[lvl].numPts() << endl;
    }
    CH_assert(a_all_tags[finest_level].isEmpty());
  }
#endif

  return got_new_tags;
}

void plasma_engine::write_geometry(){
  CH_TIME("plasma_engine::write_geometry");
  if(m_verbosity > 3){
    pout() << "plasma_engine::write_geometry" << endl;
  }

  EBAMRCellData output;
  m_amr->allocate(output, phase::gas, 1);
  data_ops::set_value(output, 0.0);
  Vector<std::string> names(1, "dummy_data");

  const int finest_level                 = m_amr->get_finest_level();
  const Vector<DisjointBoxLayout>& grids = m_amr->get_grids();
  const Vector<ProblemDomain>& domains   = m_amr->get_domains();
  const Vector<Real>& dx                 = m_amr->get_dx();
  const Vector<int>& ref_rat             = m_amr->get_ref_rat();

  bool replace_covered = false;
  Vector<Real> covered_values;

  Vector<LevelData<EBCellFAB>*> output_ptr(1 + finest_level);
  m_amr->alias(output_ptr, output);

  // Dummy file name
  char file_char[1000];
  const std::string prefix = m_output_dir + "/" + m_output_names;
  sprintf(file_char, "%s.geometry.%dd.hdf5", prefix.c_str(), SpaceDim);
  string fname(file_char);

  writeEBHDF5(fname, 
	      grids,
	      output_ptr,
	      names, 
	      domains[0],
	      dx[0], 
	      m_dt,
	      m_time,
	      ref_rat,
	      finest_level + 1,
	      replace_covered,
	      covered_values,
	      m_num_plot_ghost*IntVect::Unit);
}

void plasma_engine::write_plot_file(){
  CH_TIME("plasma_engine::write_plot_file");
  if(m_verbosity > 3){
    pout() << "plasma_engine::write_plot_file" << endl;
  }
  
  Vector<string> names = this->get_output_variable_names();
  const int num_output = names.size();

  EBAMRCellData output;
  m_amr->allocate(output, phase::gas, num_output);
  data_ops::set_value(output, 0.0);

  // Add data to output
  int cur_var = 0.0;
  this->add_potential_to_output(output,      cur_var); cur_var += 1;
  this->add_electric_field_to_output(output, cur_var); cur_var += SpaceDim;
  this->add_space_charge_to_output(output,   cur_var); cur_var += 1;
  this->add_surface_charge_to_output(output, cur_var); cur_var += 1;
  if(m_output_mode == output_mode::full){
    this->add_current_density_to_output(output, cur_var); cur_var += SpaceDim;
  }
  this->add_cdr_densities_to_output(output,  cur_var); cur_var += m_plaskin->get_num_species();;
  if(m_output_mode == output_mode::full){
    this->add_cdr_velocities_to_output(output, cur_var); cur_var += m_plaskin->get_num_species()*SpaceDim;
  }
  if(m_output_mode == output_mode::medium || m_output_mode == output_mode::full){
    this->add_cdr_source_to_output(output,     cur_var); cur_var += m_plaskin->get_num_species();
  }
  this->add_rte_densities_to_output(output, cur_var); cur_var += m_plaskin->get_num_photons();
  if(m_output_mode == output_mode::medium || m_output_mode == output_mode::full){
    this->add_rte_source_to_output(output, cur_var); cur_var += m_plaskin->get_num_photons();
  }
  if(m_output_mode == output_mode::full && !m_celltagger.isNull()){
    this->add_tracer_fields_to_output(output, cur_var); cur_var += m_celltagger->get_num_tracers();
  }

  // Filename
  char file_char[1000];
  const std::string prefix = m_output_dir + "/" + m_output_names;
  sprintf(file_char, "%s.step%07d.%dd.hdf5", prefix.c_str(), m_step, SpaceDim);
  string fname(file_char);

  // Write data file
  const int finest_level                 = m_amr->get_finest_level();
  const Vector<DisjointBoxLayout>& grids = m_amr->get_grids();
  const Vector<ProblemDomain>& domains   = m_amr->get_domains();
  const Vector<Real>& dx                 = m_amr->get_dx();
  const Vector<int>& ref_rat             = m_amr->get_ref_rat();

  bool replace_covered = false;
  Vector<Real> covered_values;

  Vector<LevelData<EBCellFAB>* > output_ptr(1 + finest_level);
  m_amr->alias(output_ptr, output);

  // Write HDF5 file
  if(m_verbosity > 3){
    pout() << "plasma_engine::write_plot_file - writing plot file..." << endl;
  }

  int plot_depth = Min(m_max_plot_depth, finest_level);
  if(m_max_plot_depth < 0){
    plot_depth = finest_level;
  }

  writeEBHDF5(fname, 
	      grids,
	      output_ptr,
	      names, 
	      domains[0],
	      dx[0], 
	      m_dt,
	      m_time,
	      ref_rat,
	      plot_depth + 1,
	      replace_covered,
	      covered_values,
	      m_num_plot_ghost*IntVect::Unit);
  
  if(m_verbosity > 3){
    pout() << "plasma_engine::write_plot_file - writing plot file... DONE!" << endl;
  }
}

void plasma_engine::write_checkpoint_file(){
  CH_TIME("plasma_engine::write_checkpoint_file");
  if(m_verbosity > 3){
    pout() << "plasma_engine::write_checkpoint_file" << endl;
  }
  
  const int finest_level = m_amr->get_finest_level();
  int finest_chk_level  = Min(m_max_chk_depth, finest_level);
  if(m_max_chk_depth < 0){
    finest_chk_level = finest_level;
  }

  const phase::which_phase cur_phase = phase::gas;
  
  RefCountedPtr<cdr_layout>& cdr   = m_timestepper->get_cdr();
  RefCountedPtr<rte_layout>& rte   = m_timestepper->get_rte();
  RefCountedPtr<sigma_solver>& sig = m_timestepper->get_sigma();

  HDF5HeaderData header;
  header.m_real["coarsest_dx"] = m_amr->get_dx()[0];
  header.m_real["time"]        = m_time;
  header.m_real["dt"]          = m_dt;
  header.m_int["step"]         = m_step;
  header.m_int["finest_level"] = finest_level;

  // Storage for sigma and m_tags - these are things that can't be written directly
  EBAMRCellData sigma, tags;
  m_amr->allocate(sigma, cur_phase, 1);
  m_amr->allocate(tags,  cur_phase, 1);

  // Copy data from sigma solver to sigma
  data_ops::set_value(sigma, 0.0);
  data_ops::incr(sigma, sig->get_state(), 1.0);



  // Set tagged cells = 1
  for (int lvl = 0; lvl <= finest_chk_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(cur_phase)[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs   = IntVectSet((*m_tags[lvl])[dit()].get_ivs());

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	(*tags[lvl])[dit()](vof, 0) = 1.;
      }
    }
  }

  // Output file name
  char str[100];
  const std::string prefix = m_output_dir + "/" + m_output_names;
  sprintf(str, "%s.check%07d.%dd.hdf5", prefix.c_str(), m_step, SpaceDim);


  HDF5Handle handle_out(str, HDF5Handle::CREATE);
  header.writeToFile(handle_out);

  for (int lvl = 0; lvl <= finest_chk_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    handle_out.setGroupToLevel(lvl);
    write(handle_out, dbl);

    for (cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
      const RefCountedPtr<cdr_solver>& solver = solver_it();
      const EBAMRCellData& state              = solver->get_state();
      const std::string name                  = solver->get_name();

      write(handle_out, *state[lvl], name);
    }

    if(!m_timestepper->stationary_rte()){ // Must write RTE data if the solvers are transient
      for (rte_iterator solver_it(*rte); solver_it.ok(); ++solver_it){
	const RefCountedPtr<rte_solver>& solver = solver_it();
	const EBAMRCellData& state              = solver->get_state();
	const std::string name                  = solver->get_name();
	
	write(handle_out, *state[lvl], name);
      }
    }

    write(handle_out, *sigma[lvl], "sigma");
    write(handle_out, *tags[lvl],  "tags");
  }

  
  handle_out.close();
}

void plasma_engine::open_mass_dump_file(ofstream& a_file){
  if(procID() == 0){
    const std::string prefix = m_output_dir + "/" + "mass_dump.txt";
    a_file.open(prefix);

    const Vector<std::string> names = m_timestepper->get_cdr()->get_names();

    a_file << "# time step" << "\t" << "time";
    for (int i = 0; i < names.size(); i++){
      a_file << "\t" << names[i];
    }
    a_file << endl;
  }
}

void plasma_engine::open_charge_dump_file(ofstream& a_file){
  if(procID() == 0){
    const std::string prefix = m_output_dir + "/" + "charge_dump.txt";
    a_file.open(prefix);

    const Vector<std::string> names = m_timestepper->get_cdr()->get_names();

    a_file << "# time step" << "\t" << "time";
    for (int i = 0; i < names.size(); i++){
      a_file << "\t" << names[i];
    }
    a_file << "\t" << "surface charge" << endl;
    a_file << endl;
  }
}

void plasma_engine::dump_mass(ofstream& a_file){
  const Vector<Real> masses = m_timestepper->get_cdr()->compute_mass();
  if(procID() == 0){
    a_file << m_step << "\t" << m_time;
    for (int i = 0; i < masses.size(); i++){
      a_file << "\t" << masses[i];
    }
    a_file << endl;
  }
}

void plasma_engine::dump_charge(ofstream& a_file){
  const Real surface_charge  = m_timestepper->get_sigma()->compute_charge();
  const Vector<Real> charges = m_timestepper->get_cdr()->compute_charge();
  if(procID() == 0){
    a_file << m_step << "\t" << m_time;
    for (int i = 0; i < charges.size(); i++){
      a_file << "\t" << charges[i]/units::s_Qe;
    }
    a_file << "\t" << surface_charge/units::s_Qe;
    a_file << endl;
  }
}

void plasma_engine::close_mass_dump_file(ofstream& a_file){
  if(procID() == 0){
    a_file.close();
  }
}

void plasma_engine::close_charge_dump_file(ofstream& a_file){
  if(procID() == 0){
    a_file.close();
  }
}

Vector<string> plasma_engine::get_output_variable_names(){
  CH_TIME("plasma_engine::get_output_variable_names");
  if(m_verbosity > 5){
    pout() << "plasma_engine::get_output_variable_names" << endl;
  }

  // Handle to solvers
  RefCountedPtr<cdr_layout>& cdr         = m_timestepper->get_cdr();
  RefCountedPtr<poisson_solver>& poisson = m_timestepper->get_poisson();
  RefCountedPtr<rte_layout>& rte         = m_timestepper->get_rte();
  RefCountedPtr<sigma_solver>& sigma     = m_timestepper->get_sigma();

  int cur_name = 0;
  int num_vars = 0;

  Vector<string> names;

  names.clear();

  // Things related to the field solver. Potential, Electric field, and surface charge density => (2 + SpaceDim)
  names.push_back("Potential"); cur_name++; num_vars++;
  names.push_back("x-Electric field"); cur_name++; num_vars++;
  names.push_back("y-Electric field"); cur_name++; num_vars++;
  if(SpaceDim == 3){
    names.push_back("z-Electric field"); cur_name++; num_vars++;
  }
  names.push_back("Space charge density"); cur_name++; num_vars++;
  names.push_back("Surface charge density"); cur_name++; num_vars++;
  if(m_output_mode == output_mode::full){
    names.push_back("x-Current density"); cur_name++; num_vars++;
    names.push_back("y-Current density"); cur_name++; num_vars++;
    if(SpaceDim == 3){
      names.push_back("z-Current density"); cur_name++; num_vars++;
    }
  }

  // Ion densities
  for (cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    RefCountedPtr<species>& spec      = solver_it.get_species();
    const std::string name = spec->get_name();
    const std::string unit = spec->get_unit();

    names.push_back(name + " ." + unit + "."); cur_name++; num_vars++;
  }

  // Ion velocities
  if(m_output_mode == output_mode::full){
    for (cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver = solver_it();
      RefCountedPtr<species>& spec      = solver_it.get_species();
      const std::string& name           = spec->get_name();
    
      names.push_back("x-Velocity " + name);   cur_name++; num_vars++;
      names.push_back("y-Velocity " + name);   cur_name++; num_vars++;
      if(SpaceDim == 3){
	names.push_back("z-Velocity " + name); cur_name++; num_vars++;
      }
    }
  }

  // Ion source terms
  if(m_output_mode == output_mode::full || m_output_mode == output_mode::medium){
    for (cdr_iterator solver_it(*cdr); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver = solver_it();
      RefCountedPtr<species>& spec      = solver_it.get_species();
      const std::string& name           = spec->get_name();
    
      names.push_back(name + " source");  cur_name++; num_vars++;
    }
  }

  // Photon solvers' stuff. Photon density and source term => (2 + SpaceDim)
  for (rte_iterator solver_it(*rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<photon_group>& photon = solver_it.get_photon();
    const std::string& name = photon->get_name();

    names.push_back(name + " density"); cur_name++; num_vars++;
  }

  // Photon solvers' stuff. Photon density and source term => (2 + SpaceDim)
  if(m_output_mode == output_mode::full || m_output_mode == output_mode::medium){
    for (rte_iterator solver_it(*rte); solver_it.ok(); ++solver_it){
      RefCountedPtr<rte_solver>& solver   = solver_it();
      RefCountedPtr<photon_group>& photon = solver_it.get_photon();
      const std::string& name = photon->get_name();

      names.push_back(name + " source"); cur_name++; num_vars++;
    }
  }

  // Tracer field(s)
  if(m_output_mode == output_mode::full){
    if(!m_celltagger.isNull()){
      const Vector<EBAMRCellData>& tracer = m_celltagger->get_tracer_fields();
      for (int i = 0; i < tracer.size(); i++){
	std::string one = "Tracer field-";

	long int j = i;
	char s[2]; 
	sprintf(s,"%ld", j);

	std::string two(s);
	names.push_back(one + two); cur_name++; num_vars++;
      }
    }
  }

  return names;
}

void plasma_engine::compute_norm(std::string a_chk_coarse, std::string a_chk_fine){
  CH_TIME("plasma_engine::compute_norm");
  if(m_verbosity > 5){
    pout() << "plasma_engine::compute_norm" << endl;
  }

  // TLDR: This routine is long and ugly. In short, it instantiates a all solvers and read the fine-level checkpoint file
  // into those solvers.
  RefCountedPtr<cdr_layout>& cdr         = m_timestepper->get_cdr();

  this->read_checkpoint_file(a_chk_fine); // This reads the fine-level data

  Vector<EBISLayout>& fine_ebisl        = m_amr->get_ebisl(phase::gas);
  Vector<DisjointBoxLayout>& fine_grids = m_amr->get_grids();
  Vector<ProblemDomain>& fine_domains   = m_amr->get_domains();

  // Read the second file header
  HDF5Handle handle_in(a_chk_coarse, HDF5Handle::OPEN_RDONLY);
  HDF5HeaderData header;
  header.readFromFile(handle_in);
  int finest_coar_level = header.m_int["finest_level"];


  // Read coarse data into these data holders
  Vector<EBAMRCellData> cdr_densities(1 + m_plaskin->get_num_species());
  for (int i = 0; i < cdr_densities.size(); i++){
    cdr_densities[i].resize(1+finest_coar_level);
  }

  // Read cdr data
  for (int lvl = 0; lvl <= finest_coar_level; lvl++){
    handle_in.setGroupToLevel(lvl);
    
    for (cdr_iterator solver_it = cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      RefCountedPtr<cdr_solver>& solver = solver_it();
      const std::string solver_name = solver->get_name();
      
      EBCellFactory cellfact(fine_ebisl[lvl]);
      cdr_densities[idx][lvl] = RefCountedPtr<LevelData<EBCellFAB> > 
	(new LevelData<EBCellFAB>(fine_grids[lvl], 1, 3*IntVect::Unit, cellfact));
      
      read<EBCellFAB>(handle_in, *cdr_densities[idx][lvl], solver_name, fine_grids[lvl], Interval(), false);
    }
  }


  // We will now compute n_2h - A(n_h)
  pout() << "files are '"  << a_chk_coarse << "' and '" << a_chk_fine << "'" << endl;
  EBCellFactory cellfact(fine_ebisl[finest_coar_level]);
  LevelData<EBCellFAB> diff(fine_grids[finest_coar_level], 1, 0*IntVect::Unit, cellfact);
  for (cdr_iterator solver_it = cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<cdr_solver>& solver = solver_it();
    const std::string solver_name = solver->get_name();

    LevelData<EBCellFAB>& fine_sol = *solver->get_state()[finest_coar_level];
    LevelData<EBCellFAB>& coar_sol = *cdr_densities[idx][finest_coar_level];

    data_ops::set_value(diff, 0.0);
    data_ops::incr(diff, coar_sol,   1.0);
    data_ops::incr(diff, fine_sol,  -1.0);
    data_ops::scale(diff, 1.E-18);

    writeEBLevelname (&diff, "diff.hdf5");
    Real volume;

    Real Linf;
    Real L1;
    Real L2;

    Linf = EBLevelDataOps::kappaNorm(volume, diff, EBLEVELDATAOPS_ALLVOFS, fine_domains[finest_coar_level], 0);
    L1   = EBLevelDataOps::kappaNorm(volume, diff, EBLEVELDATAOPS_ALLVOFS, fine_domains[finest_coar_level], 1);
    L2   = EBLevelDataOps::kappaNorm(volume, diff, EBLEVELDATAOPS_ALLVOFS, fine_domains[finest_coar_level], 2);
    
    pout() << idx << "\t Linf = " << Linf << "\t L1 = " << L1 << "\t L2 = " << L2 << endl;
  }
  
  handle_in.close();
}
