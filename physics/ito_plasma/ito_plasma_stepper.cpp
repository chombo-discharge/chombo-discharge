/*!
  @file   ito_plasma_stepper.cpp
  @brief  Implementation of ito_plasma_stepper.H
  @author Robert Marskar
  @date   May 2020
*/

#include "ito_plasma_stepper.H"
#include "data_ops.H"
#include "units.H"

#include <EBArith.H>
#include <PolyGeom.H>

using namespace physics::ito_plasma;

ito_plasma_stepper::ito_plasma_stepper(){
  m_verbosity = -1;
  m_name      = "ito_plasma_stepper";
  m_phase     = phase::gas;

  m_dt   = 0.0;
  m_time = 0.0;

  m_regrid_superparticles = false;

  m_fluid_realm    = realm::primal;
  m_particle_realm = realm::primal;
}

ito_plasma_stepper::ito_plasma_stepper(RefCountedPtr<ito_plasma_physics>& a_physics) : ito_plasma_stepper(){
  m_physics   = a_physics;
}

ito_plasma_stepper::~ito_plasma_stepper(){
}

void ito_plasma_stepper::set_verbosity(const int a_verbosity){
  CH_TIME("ito_plasma_stepper::set_verbosity");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_verbosity" << endl;
  }
  m_verbosity = a_verbosity;
}

void ito_plasma_stepper::setup_solvers(){
  CH_TIME("ito_plasma_stepper::setup_solver");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_solvers" << endl;
  }

  // Parse class options
  this->parse_options();

  // Set up solvers
  this->setup_ito();
  this->setup_poisson();
  this->setup_rte();
  this->setup_sigma();

  // Allocate internal stuff
  this->allocate_internals();
}

void ito_plasma_stepper::setup_ito(){
  CH_TIME("ito_plasma_stepper::setup_ito");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_ito" << endl;
  }

  m_ito->set_verbosity(m_verbosity);
  m_ito->parse_options();
  m_ito->set_amr(m_amr);
  m_ito->set_phase(m_phase);
  m_ito->set_computational_geometry(m_compgeom);
  m_ito->set_realm(m_particle_realm);
}

void ito_plasma_stepper::setup_poisson(){
  CH_TIME("ito_plasma_stepper::setup_poisson");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_poisson" << endl;
  }

  m_poisson->set_verbosity(m_verbosity);
  m_poisson->parse_options();
  m_poisson->set_amr(m_amr);
  m_poisson->set_computational_geometry(m_compgeom);
  m_poisson->set_potential(m_potential); // Needs to happen AFTER set_poisson_wall_func
  m_poisson->set_realm(m_fluid_realm);
  m_poisson->sanity_check();
}

void ito_plasma_stepper::setup_rte(){
  CH_TIME("ito_plasma_stepper::setup_rte");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_rte" << endl;
  }

  m_rte->set_verbosity(m_verbosity);
  m_rte->parse_options();
  m_rte->set_phase(m_phase);
  m_rte->set_amr(m_amr);
  m_rte->set_computational_geometry(m_compgeom);
  m_rte->set_realm(m_particle_realm);
  m_rte->sanity_check();
}

void ito_plasma_stepper::setup_sigma(){
  CH_TIME("ito_plasma_stepper::setup_sigma");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_sigma" << endl;
  }

  m_sigma = RefCountedPtr<sigma_solver> (new sigma_solver());
  m_sigma->set_amr(m_amr);
  m_sigma->set_verbosity(m_verbosity);
  m_sigma->set_computational_geometry(m_compgeom);
  m_sigma->set_realm(m_fluid_realm);
}

void ito_plasma_stepper::allocate() {
  CH_TIME("ito_plasma_stepper::allocate");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::allocate" << endl;
  }

  m_ito->allocate_internals();
  m_rte->allocate_internals();
  m_poisson->allocate_internals();
  m_sigma->allocate_internals();
}

void ito_plasma_stepper::post_initialize(){

}

void ito_plasma_stepper::initial_data(){
  CH_TIME("ito_plasma_stepper::initial_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::initial_data" << endl;
  }

  m_ito->initial_data(); // This deposits, of course. 
  m_rte->initial_data();
  this->initial_sigma();

  m_ito->sort_particles_by_cell();
  m_ito->make_superparticles(m_ppc);
  m_ito->sort_particles_by_patch();
  
  // Solve Poisson equation and compute the E-field
  this->solve_poisson();

  // Fill solvers with velocities and diffusion coefficients
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
}

void ito_plasma_stepper::initial_sigma(){
  CH_TIME("ito_plasma_stepper::initial_sigma");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::initial_sigma" << endl;
  }

  const RealVect origin  = m_amr->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  EBAMRIVData& sigma = m_sigma->get_state();
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_fluid_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_fluid_realm, phase::gas)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      BaseIVFAB<Real>& state = (*sigma[lvl])[dit()];

      const EBISBox& ebisbox = ebisl[dit()];
      const IntVectSet& ivs  = state.getIVS();
      const EBGraph& ebgraph = state.getEBGraph();
      
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = origin + vof.gridIndex()*dx + 0.5*ebisbox.bndryCentroid(vof)*dx;
	
	for (int comp = 0; comp < state.nComp(); comp++){
	  state(vof, comp) = m_physics->initial_sigma(m_time, pos);
	}
      }
    }
  }

  m_amr->average_down(sigma, m_fluid_realm, phase::gas);
  m_sigma->reset_cells(sigma);
}

void ito_plasma_stepper::post_checkpoint_setup(){
  CH_TIME("ito_plasma_stepper::post_checkpoint_setup");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::post_checkpoint_setup" << endl;
  }

  // Recompute poisson
  this->solve_poisson();
  this->allocate_internals();

  
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
}

void ito_plasma_stepper::write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const {
  CH_TIME("ito_plasma_stepper::write_checkpoint_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::write_checkpoint_data" << endl;
  }

  for (ito_iterator<ito_solver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<ito_solver>& solver = solver_it();
    solver->write_checkpoint_level(a_handle, a_lvl);
  }

  for (rte_iterator<mc_photo> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<mc_photo>& solver = solver_it();
    solver->write_checkpoint_level(a_handle, a_lvl);
  }

  m_poisson->write_checkpoint_level(a_handle, a_lvl);
  m_sigma->write_checkpoint_level(a_handle, a_lvl);
}

void ito_plasma_stepper::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl){
  CH_TIME("ito_plasma_stepper::read_checkpoint_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::read_checkpoint_data" << endl;
  }

  for (ito_iterator<ito_solver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    solver->read_checkpoint_level(a_handle, a_lvl);
  }

  for (rte_iterator<mc_photo> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<mc_photo>& solver = solver_it();
    solver->read_checkpoint_level(a_handle, a_lvl);
  }

  m_poisson->read_checkpoint_level(a_handle, a_lvl);
  m_sigma->read_checkpoint_level(a_handle, a_lvl);
}

void ito_plasma_stepper::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const {
  CH_TIME("ito_plasma_stepper::write_plot_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::write_plot_data" << endl;
  }

  // Poisson solver copies over its output data
  a_plotvar_names.append(m_poisson->get_plotvar_names());
  m_poisson->write_plot_data(a_output, a_icomp);

  // Surface charge solver writes
  a_plotvar_names.append(m_sigma->get_plotvar_names());
  m_sigma->write_plot_data(a_output, a_icomp);

  // Ito solvers copy their output data
  for (ito_iterator<ito_solver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    a_plotvar_names.append(solver->get_plotvar_names());
    solver->write_plot_data(a_output, a_icomp);
  }

  // RTE solvers copy their output data
  for (rte_iterator<mc_photo> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<mc_photo>& solver = solver_it();
    a_plotvar_names.append(solver->get_plotvar_names());
    solver->write_plot_data(a_output, a_icomp);
  }

  // Write the current to the output
  this->write_J(a_output, a_icomp);
  a_plotvar_names.push_back("x-J");
  a_plotvar_names.push_back("y-J");
  if(SpaceDim == 3){
    a_plotvar_names.push_back("z-J");
  }
}

void ito_plasma_stepper::write_J(EBAMRCellData& a_output, int& a_icomp) const{
  CH_TIME("ito_plasma_stepper::write_J");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::write_J" << endl;
  }

  const Interval src_interv(0, SpaceDim-1);
  const Interval dst_interv(a_icomp, a_icomp + SpaceDim -1);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    if(m_J.get_realm() == a_output.get_realm()){
      m_J[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else {
      m_J[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }
  a_icomp += SpaceDim;
}

void ito_plasma_stepper::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt){
  CH_TIME("ito_plasma_stepper::synchronize_solver_times");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::synchronize_solver_times" << endl;
  }

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;

  m_ito->set_time(a_step,     a_time, a_dt);
  m_poisson->set_time(a_step, a_time, a_dt);
  m_rte->set_time(a_step,     a_time, a_dt);
  m_sigma->set_time(a_step,   a_time, a_dt);
}

void ito_plasma_stepper::print_step_report(){
  CH_TIME("ito_plasma_stepper::print_step_report");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::print_step_report" << endl;
  }

  const Real Emax = this->compute_Emax(m_phase);
  
  const size_t l_particles         = m_ito->get_num_particles(true);
  const size_t g_particles        = m_ito->get_num_particles(false);
  
  const size_t l_eb_particles      = m_ito->get_num_eb_particles(true);
  const size_t g_eb_particles     = m_ito->get_num_eb_particles(false);
  
  const size_t l_domain_particles  = m_ito->get_num_domain_particles(true);
  const size_t g_domain_particles = m_ito->get_num_domain_particles(false);

  const size_t l_source_particles  = m_ito->get_num_source_particles(true);
  const size_t g_source_particles = m_ito->get_num_source_particles(false);


  pout() << "                                   Emax      = " << Emax << endl
	 << "                                   #part     = " << l_particles << " (" << g_particles << ")" << endl
	 << "                                   #eb part  = " << l_eb_particles << " (" << g_eb_particles << ")" << endl
	 << "                                   #dom part = " << l_domain_particles << " (" << g_domain_particles << ")" << endl
    	 << "                                   #src part = " << l_source_particles << " (" << g_source_particles << ")" << endl;
}

void ito_plasma_stepper::compute_dt(Real& a_dt, time_code& a_timecode){
  CH_TIME("ito_plasma_stepper::compute_dt");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_dt" << endl;
  }
  
  a_dt = m_ito->compute_dt();
  a_dt = a_dt*m_max_cells_hop;
  a_timecode = time_code::cfl;
  
}

void ito_plasma_stepper::register_realms(){
  CH_TIME("ito_plasma_stepper::register_realms");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::register_realms" << endl;
  }

  m_amr->register_realm(m_fluid_realm);
  m_amr->register_realm(m_particle_realm);
}

void ito_plasma_stepper::register_operators(){
  CH_TIME("ito_plasma_stepper::register_operators");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::register_operators" << endl;
  }

  m_ito->register_operators();
  m_poisson->register_operators();
  m_rte->register_operators();
  m_sigma->register_operators();
}
  
void ito_plasma_stepper::pre_regrid(const int a_lmin, const int a_old_finest_level){
  CH_TIME("ito_plasma_stepper::pre_regrid");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::pre_regrid" << endl;
  }

  m_ito->pre_regrid(a_lmin,     a_old_finest_level);
  m_poisson->pre_regrid(a_lmin, a_old_finest_level);
  m_rte->pre_regrid(a_lmin,     a_old_finest_level);
  m_sigma->pre_regrid(a_lmin,   a_old_finest_level);
}

void ito_plasma_stepper::deallocate(){
  CH_TIME("ito_plasma_stepper::deallocate");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::deallocate" << endl;
  }

  // Don't deallocate anything. 
}

void ito_plasma_stepper::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("ito_plasma_stepper::regrid");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::regrid" << endl;
  }

  // Allocate new memory
  this->allocate_internals();

  // Regrid solvers
  m_ito->regrid(a_lmin,     a_old_finest_level, a_new_finest_level);
  m_poisson->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  m_rte->regrid(a_lmin,     a_old_finest_level, a_new_finest_level);
  m_sigma->regrid(a_lmin,   a_old_finest_level, a_new_finest_level);

  if(m_regrid_superparticles){
    m_ito->sort_particles_by_cell();
    m_ito->make_superparticles(m_ppc);
    m_ito->sort_particles_by_patch();
  }

  // Redeposit particles
  m_ito->deposit_particles();

  // Recompute the electric field
  const bool converged = this->solve_poisson();
  if(!converged){
    MayDay::Abort("ito_plasma_stepper::regrid - Poisson solve did not converge after regrid!!!");
  }

  // Recompute new velocities and diffusion coefficients
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
}

void ito_plasma_stepper::post_regrid(){

}

int ito_plasma_stepper::get_num_plot_vars() const {
  CH_TIME("ito_plasma_stepper::get_num_plot_vars");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::get_num_plot_vars" << endl;
  }

  int ncomp = 0;
  
  for (ito_iterator<ito_solver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    ncomp += solver->get_num_plotvars();
  }
  
  for (rte_iterator<mc_photo> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<mc_photo>& solver = solver_it();
    ncomp += solver->get_num_plotvars();
  }

  ncomp += m_poisson->get_num_plotvars();
  ncomp += m_sigma->get_num_plotvars();
  ncomp += SpaceDim; // For plotting the current density

  return ncomp;
}

void ito_plasma_stepper::set_ito(RefCountedPtr<ito_layout<ito_solver> >& a_ito){
  CH_TIME("ito_plasma_stepper::set_ito");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_ito" << endl;
  }

  m_ito = a_ito;
}

void ito_plasma_stepper::set_poisson(RefCountedPtr<poisson_solver>& a_poisson){
  CH_TIME("ito_plasma_stepper::set_poisson");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_poisson" << endl;
  }

  m_poisson = a_poisson;
}

void ito_plasma_stepper::set_rte(RefCountedPtr<rte_layout<mc_photo> >& a_rte){
  CH_TIME("ito_plasma_stepper::set_rte");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_rte" << endl;
  }
  
  m_rte = a_rte;
}

void ito_plasma_stepper::set_potential(Real (*a_potential)(const Real a_time)){
  CH_TIME("ito_plasma_stepper::set_potential");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_potential" << endl;
  }

  m_potential = a_potential;
}

Real ito_plasma_stepper::compute_Emax(const phase::which_phase a_phase) {
  CH_TIME("ito_plasma_stepper::compute_Emax");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_Emax" << endl;
  }

  // Get a handle to the E-field
  EBAMRCellData Ephase;
  m_amr->allocate_ptr(Ephase);
  m_amr->alias(Ephase, m_phase, m_poisson->get_E());

  // Interpolate to centroids
  EBAMRCellData E;
  m_amr->allocate(E, m_fluid_realm, m_phase, SpaceDim);
  data_ops::copy(E, Ephase);
  m_amr->interpolate_to_centroids(E, m_fluid_realm, m_phase);

  Real max, min;
  data_ops::get_max_min_norm(max, min, E);

  return max;
}

Real ito_plasma_stepper::get_time() const{
  CH_TIME("ito_plasma_stepper::get_time");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::get_time" << endl;
  }

  return m_time;
}

void ito_plasma_stepper::compute_E(MFAMRCellData& a_E, const MFAMRCellData& a_potential){
  CH_TIME("ito_plasma_stepper::compute_E(mfamrcell,mfamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(mfamrcell, mfamrcell" << endl;
  }

  m_amr->compute_gradient(a_E, a_potential, m_fluid_realm);
  data_ops::scale(a_E, -1.0);

  m_amr->average_down(a_E, m_fluid_realm);
  m_amr->interp_ghost(a_E, m_fluid_realm);

}

void ito_plasma_stepper::compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase){
  CH_TIME("ito_plasma_stepper::compute_E(ebamrcell, phase)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(ebamrcell, phase" << endl;
  }

  this->compute_E(a_E, a_phase, m_poisson->get_state());
}

void ito_plasma_stepper::compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase, const MFAMRCellData& a_potential){
  CH_TIME("ito_plasma_stepper::compute_E(ebamrcell, phase, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(ebamrcell, phase mfamrcell" << endl;
  }
  
  EBAMRCellData pot_gas;
  m_amr->allocate_ptr(pot_gas);
  m_amr->alias(pot_gas, a_phase, a_potential);

  m_amr->compute_gradient(a_E, pot_gas, m_fluid_realm, a_phase);
  data_ops::scale(a_E, -1.0);

  m_amr->average_down(a_E, m_fluid_realm, a_phase);
  m_amr->interp_ghost(a_E, m_fluid_realm, a_phase);
}

void ito_plasma_stepper::compute_E(EBAMRFluxData& a_E_face, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell){
  CH_TIME("ito_plasma_stepper::compute_E(ebamrflux, phase, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(ebamrflux, phase mfamrcell" << endl;
  }

  CH_assert(a_E_face[0]->nComp() == SpaceDim);
  CH_assert(a_E_cell[0]->nComp() == SpaceDim);

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){

    const DisjointBoxLayout& dbl = m_amr->get_grids(m_fluid_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_fluid_realm, a_phase)[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBCellFAB& E_cell = (*a_E_cell[lvl])[dit()];
      const EBISBox& ebisbox  = ebisl[dit()];
      const EBGraph& ebgraph  = ebisbox.getEBGraph();
      const Box& box          = dbl.get(dit());
      
      for (int dir = 0; dir < SpaceDim; dir++){
      	EBFaceFAB& E_face = (*a_E_face[lvl])[dit()][dir];
	E_face.setVal(0.0);

      	EBLevelDataOps::averageCellToFace(E_face,
      					  E_cell,
      					  ebgraph,
      					  box,
      					  0,
      					  dir,
      					  domain,
      					  dir,
      					  dir);
      }

    }
    a_E_face[lvl]->exchange();
  }
}

void ito_plasma_stepper::compute_E(EBAMRIVData& a_E_eb,  const phase::which_phase a_phase, const EBAMRCellData& a_E_cell){
  CH_TIME("ito_plasma_stepper::compute_E(ebamriv, phase, ebamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(ebamriv, phase ebamrcell)" << endl;
  }

  CH_assert(a_E_eb[0]->nComp()   == SpaceDim);
  CH_assert(a_E_cell[0]->nComp() == SpaceDim);

  const irreg_amr_stencil<eb_centroid_interp>& interp_stencil = m_amr->get_eb_centroid_interp_stencils(m_fluid_realm, a_phase);
  interp_stencil.apply(a_E_eb, a_E_cell);
}

void ito_plasma_stepper::compute_rho(){
  CH_TIME("ito_plasma_stepper::compute_rho()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_rho()" << endl;
  }
  
  this->compute_rho(m_poisson->get_source(), m_ito->get_densities());
}

void ito_plasma_stepper::compute_rho(MFAMRCellData& a_rho, const Vector<EBAMRCellData*>&  a_densities){
  CH_TIME("ito_plasma_stepper::compute_rho(rho, densities)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_rho(rho, densities)" << endl;
  }

  // TLDR: a_densities is from the ito solvers so it is defined over the particle realm. But a_rho is defined over
  //       the fluid realm so we need scratch storage we can copy into. We use m_fluid_scratch1 for that. 

  // Reset
  data_ops::set_value(a_rho, 0.0);

  // Make alias
  EBAMRCellData rhoPhase;
  m_amr->allocate_ptr(rhoPhase);
  m_amr->alias(rhoPhase, m_phase, a_rho);

  // Increment each solver
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<ito_solver>& solver   = solver_it();
    const RefCountedPtr<ito_species>& species = solver->get_species();
    const int idx = solver_it.get_solver();
    const int q   = species->get_charge();

    if(species->get_charge() != 0){
      m_fluid_scratch1.copy(*a_densities[idx]);
      data_ops::incr(rhoPhase, m_fluid_scratch1, q);
    }
  }

  data_ops::scale(a_rho, units::s_Qe);

  m_amr->average_down(a_rho, m_fluid_realm);
  m_amr->interp_ghost(a_rho, m_fluid_realm);

  // Interpolate to centroids
  m_amr->interpolate_to_centroids(rhoPhase, m_fluid_realm, m_phase);
}

void ito_plasma_stepper::compute_J(EBAMRCellData& a_J, const Real a_dt){
  CH_TIME("ito_plasma_stepper::compute_J(J)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_J(J)" << endl;
  }

  // TLDR: a_J is defined over the fluid realm but the computation takes place on the particle realm.
  //       If the realms are different we compute on a scratch storage instead

  data_ops::set_value(a_J, 0.0);

  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    LevelData<EBCellFAB>* data;
    if(lvl > 0){ 
      RefCountedPtr<EBMGInterp>& interp = m_amr->get_eb_mg_interp(m_particle_realm, m_phase)[lvl];
      interp->pwcInterp(*m_particle_scratchD[lvl], *m_particle_scratchD[lvl-1], Interval(0, SpaceDim-1));
    }
    this->compute_J(*m_particle_scratchD[lvl], lvl, a_dt);
  }

  // Particle realm to fluid realm copy
  a_J.copy(m_particle_scratchD);
  
  m_amr->average_down(a_J, m_fluid_realm, m_phase);
  m_amr->interp_ghost(a_J, m_fluid_realm, m_phase);

  m_amr->interpolate_to_centroids(a_J, m_fluid_realm, m_phase);
}

void ito_plasma_stepper::compute_J(LevelData<EBCellFAB>& a_J, const int a_level, const Real a_dt){
  CH_TIME("ito_plasma_stepper::compute_J(J, level)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_J(J, level)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box = dbl.get(dit());

    // Add current on present level. 
    this->compute_J(a_J[dit()], a_level, dit(), box, a_dt);
  }
}

void ito_plasma_stepper::compute_J(EBCellFAB& a_J, const int a_level, const DataIndex a_dit, const Box& a_box, const Real a_dt){
  CH_TIME("ito_plasma_stepper::compute_J(J, level, dit)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_J(J, level, dit)" << endl;
  }

  const int vofId = 0;
  const RealVect prob_lo = m_amr->get_prob_lo();
  const RealVect dx = m_amr->get_dx()[a_level]*RealVect::Unit;
  const Real idV = 1./pow(m_amr->get_dx()[a_level], SpaceDim);

  // TLDR: This code computes q*(Xnew - Xold)/dt for all charged particles
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<ito_solver>& solver   = solver_it();
    const RefCountedPtr<ito_species>& species = solver->get_species();

    if(solver->is_mobile() || solver->is_diffusive() && species->get_charge() != 0){
      const List<ito_particle>& particles = solver->get_particles()[a_level][a_dit].listItems();
      
      for (ListIterator<ito_particle> lit(particles); lit.ok(); ++lit){
	const ito_particle& p = lit();
	const RealVect rv = (p.position() - prob_lo)/dx;
	const IntVect iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

	const VolIndex vof(iv, vofId);

	const Real weight = p.mass();
	const RealVect v  = (p.position() - p.oldPosition())/a_dt;

	for (int dir = 0; dir < SpaceDim; dir++){
	  a_J(vof, dir) += units::s_Qe*weight*v[dir]*idV;
	}
      }
    }
  }
}

Real ito_plasma_stepper::compute_relaxation_time(){
  CH_TIME("ito_plasma_stepper::compute_relaxation_time()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_relaxation_time()" << endl;
  }

  Real dt = 1.E99;
  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const Real thisDt = this->compute_relaxation_time(lvl);

    dt = Min(dt, thisDt);
  }

  return dt;
}

Real ito_plasma_stepper::compute_relaxation_time(const int a_level){
  CH_TIME("ito_plasma_stepper::compute_relaxation_time(level)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_relaxation_time(level)" << endl;
  }

  Real dt = 1.E99;

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_fluid_realm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Real thisDt = this->compute_relaxation_time(a_level, dit());

    dt = Min(dt, thisDt);
  }

#ifdef CH_MPI
  Real tmp = dt;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("cdr_solver::compute_cfl_dt() - communication error on norm");
  }
  dt = tmp;
#endif

  return dt;
}

Real ito_plasma_stepper::compute_relaxation_time(const int a_level, const DataIndex a_dit){
  CH_TIME("ito_plasma_stepper::compute_relaxation_time(level, dit)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_relaxation_time(level, dit)" << endl;
  }

  const int comp    = 0;
  const Real SAFETY = 1.E-10;

  const Box box = m_amr->get_grids(m_fluid_realm)[a_level].get(a_dit);
  const EBISBox& ebisbox = m_amr->get_ebisl(m_fluid_realm, m_phase)[a_level][a_dit];

  // Get a handle to the E-field
  EBAMRCellData amrE;
  m_amr->allocate_ptr(amrE);
  m_amr->alias(amrE, m_phase, m_poisson->get_E());
  
  const EBCellFAB& E = (*amrE[a_level])[a_dit];
  const EBCellFAB& J = (*m_J[a_level])[a_dit];

  EBCellFAB dt(ebisbox, box, 1);
  EBCellFAB e_magnitude(ebisbox, box, 1);
  EBCellFAB j_magnitude(ebisbox, box, 1);

  e_magnitude.setVal(0.0);
  j_magnitude.setVal(0.0);

  data_ops::vector_length(e_magnitude, E, box);
  data_ops::vector_length(j_magnitude, J, box);
  j_magnitude += SAFETY;

  dt.setVal(units::s_eps0);
  dt *= e_magnitude;
  dt /= j_magnitude;
  
  return dt.min(comp);
}

bool ito_plasma_stepper::solve_poisson(){
  CH_TIME("ito_plasma_stepper::solve_poisson()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::solve_poisson()" << endl;
  }

  this->compute_rho(); // This computes rho onto m_poisson->get_source()
  const bool converged = m_poisson->solve(m_poisson->get_state(),
					  m_poisson->get_source(),
					  m_sigma->get_state(),
					  false);

  m_poisson->compute_E();
    
  return converged;
}

bool ito_plasma_stepper::solve_poisson(MFAMRCellData&                a_potential,
				       MFAMRCellData&                a_rho,
				       const Vector<EBAMRCellData*>& a_densities,
				       const EBAMRIVData&            a_sigma){
  CH_TIME("ito_plasma_stepper::solve_poisson(phi, rho, densities, sigma)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::solve_poisson(phi, rho, densities, sigma)" << endl;
  }

  this->compute_rho(a_rho, a_densities);

  const bool converged = m_poisson->solve(a_potential,
					  a_rho,
					  a_sigma,
					  false);
  m_poisson->compute_E();
  
  return converged;
}

void ito_plasma_stepper::deposit_particles(){
  CH_TIME("ito_plasma_stepper::deposit_particles");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::deposit_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->deposit_particles();
  }
}

void ito_plasma_stepper::deposit_mobile_particles(){
  CH_TIME("ito_plasma_stepper::deposit_mobile_particles");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::deposit_mobile_particles" << endl;
  }
  
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    if(solver_it()->is_mobile()){
      solver_it()->deposit_particles();
    }
  }
}

void ito_plasma_stepper::deposit_diffusive_particles(){
  CH_TIME("ito_plasma_stepper::deposit_diffusive_particles");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::deposit_diffusive_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    if(solver_it()->is_diffusive()){
      solver_it()->deposit_particles();
    }
  }
}

void ito_plasma_stepper::deposit_mobile_or_diffusive_particles(){
  CH_TIME("ito_plasma_stepper::deposit_mobile_or_diffusive_particles");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::deposit_mobile_or_diffusive_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();
    
    if(mobile || diffusive){
      solver->deposit_particles();
    }
  }
}

void ito_plasma_stepper::deposit_stationary_particles(){
  CH_TIME("ito_plasma_stepper::deposit_stationary_particles");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::deposit_stationary_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();
    
    if(!mobile && !diffusive){
      solver->deposit_particles();
    }
  }
}

void ito_plasma_stepper::remap_particles(){
  CH_TIME("ito_plasma_stepper::remap_particles");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::remap_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    solver->remap();
  }
}

void ito_plasma_stepper::remap_mobile_particles(){
  CH_TIME("ito_plasma_stepper::remap_mobile_particles");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::remap_mobile_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    if(solver->is_mobile()){
      solver->remap();
    }
  }
}

void ito_plasma_stepper::remap_diffusive_particles(){
  CH_TIME("ito_plasma_stepper::remap_diffusive_particles");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::remap_diffusive_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    if(solver->is_diffusive()){
      solver->remap();
    }
  }
}

void ito_plasma_stepper::remap_mobile_or_diffusive_particles(){
  CH_TIME("ito_plasma_stepper::remap_mobile_or_diffusive_particles");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::remap_mobile_or_diffusive_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();
    
    if(mobile || diffusive){
      solver->remap();
    }
  }
}

void ito_plasma_stepper::compute_ito_velocities(){
  CH_TIME("ito_plasma_stepper::compute_ito_velocities()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_velocities()" << endl;
  }
  
  EBAMRCellData E;
  m_amr->allocate_ptr(E);
  m_amr->alias(E, m_phase, m_poisson->get_E());

  Vector<EBAMRCellData*> velocities = m_ito->get_velocities();
  Vector<EBAMRCellData*> densities  = m_ito->get_densities();

  this->compute_ito_velocities(velocities, densities, E, m_time);
}

void ito_plasma_stepper::compute_ito_velocities(Vector<EBAMRCellData*>&       a_velo,
						const Vector<EBAMRCellData*>& a_densities,
						const EBAMRCellData&          a_E,
						const Real                    a_time){
  CH_TIME("ito_plasma_stepper::compute_ito_velocities(velo, E, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_velocities(velo, E, time)" << endl;
  }

  // TLDR: The a_E that comes in is fluid_realm. We copy that to equivalent storage on the particle realm
  m_particle_scratchD.copy(a_E);
  m_amr->average_down(m_particle_scratchD, m_particle_realm, m_phase);
  m_amr->interp_ghost(m_particle_scratchD, m_particle_realm, m_phase);

  // Interpolate to centroids
  //  m_amr->interpolate_to_centroids(m_particle_scratchD, m_particle_realm, m_phase);

  const int num_ito_species = m_physics->get_num_ito_species();
  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){

    Vector<LevelData<EBCellFAB>* > velocities(num_ito_species);
    Vector<LevelData<EBCellFAB>* > densities(num_ito_species);
    
    for (int idx = 0; idx < a_velo.size(); idx++){
      velocities[idx] = (*a_velo[idx])[lvl];
      densities[idx]  = (*a_densities[idx])[lvl];
    }

    this->compute_ito_velocities(velocities, densities, *m_particle_scratchD[lvl], lvl, a_time);
  }

  // Average down, interpolate ghost cells, and then interpolate to particle positions
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<ito_solver>& solver = solver_it();

    if(solver->is_mobile()){
      m_amr->average_down(*a_velo[idx], m_particle_realm, m_phase);
      m_amr->interp_ghost(*a_velo[idx], m_particle_realm, m_phase);

      solver->interpolate_velocities();
    }
  }
}

void ito_plasma_stepper::compute_ito_velocities(Vector<LevelData<EBCellFAB>* >&       a_velo,
						const Vector<LevelData<EBCellFAB>* >& a_densities,
						const LevelData<EBCellFAB>&           a_E,
						const int                             a_level,
						const Real                            a_time){
  CH_TIME("ito_plasma_stepper::compute_ito_velocities(velo, E, level, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_velocities(velo, E, level, time)" << endl;
  }

  const int num_ito_species = m_physics->get_num_ito_species();

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());

    Vector<EBCellFAB*> velocities(num_ito_species);
    Vector<EBCellFAB*> densities(num_ito_species);;
    
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();

      if(solver_it()->is_mobile()){
	velocities[idx] = &(*a_velo[idx])[dit()];
      }
      densities[idx] = &(*a_densities[idx])[dit()];
    }

    this->compute_ito_velocities(velocities, densities, a_E[dit()], a_level, dit(), box, a_time);
  }
}

void ito_plasma_stepper::compute_ito_velocities(Vector<EBCellFAB*>&       a_velo,
						const Vector<EBCellFAB*>& a_densities,
						const EBCellFAB&          a_E,
						const int                 a_level,
						const DataIndex           a_dit,
						const Box                 a_box,
						const Real                a_time){
  CH_TIME("ito_plasma_stepper::compute_ito_velocities(velo, E, level, dit, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_velocities(velo, E, level, dit, time)" << endl;
  }

  const int comp         = 0;
  const Real dx          = m_amr->get_dx()[a_level];
  const RealVect prob_lo = m_amr->get_prob_lo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();

  // Do regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv   = bit();
    const RealVect pos = m_amr->get_prob_lo() + dx*(RealVect(iv) + 0.5*RealVect::Unit);
    const RealVect e   = RealVect(D_DECL(E(iv,0), E(iv, 1), E(iv, 2)));

    // Make grid densities
    Vector<Real> densities;
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      densities.push_back((*a_densities[idx]).getSingleValuedFAB()(iv, comp));
    }

    // Call ito_physics and compute velocities for each particle species
    Vector<RealVect> velocities = m_physics->compute_ito_velocities(a_time, pos, e, densities);
    
    // Put velocities where they belong
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<ito_solver>& solver = solver_it();
      if(solver->is_mobile()){
	const int idx = solver_it.get_solver();
	for (int dir = 0; dir < SpaceDim; dir++){
	  (*a_velo[idx]).getSingleValuedFAB()(iv, dir) = velocities[idx][dir];
	}
      }
    }
  }

  // Do irregular cells
  VoFIterator& vofit = (*m_amr->get_vofit(m_particle_realm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect e    = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, prob_lo);
    
    // Get densities
    Vector<Real> densities;
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      densities.push_back((*a_densities[idx])(vof, comp));
    }
    
    // Compute velocities
    Vector<RealVect> velocities = m_physics->compute_ito_velocities(a_time, pos, e, densities);

    // Put velocities in the appropriate place.
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      if(solver_it()->is_mobile()){
	const int idx = solver_it.get_solver();
	for (int dir = 0; dir < SpaceDim; dir++){
	  (*a_velo[idx])(vof, dir) = velocities[idx][dir];
	}
      }
    }
  }

  // Covered is bogus.
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    if(solver_it()->is_mobile()){
      const int idx = solver_it.get_solver();
      for (int dir = 0; dir < SpaceDim; dir++){
	a_velo[idx]->setCoveredCellVal(0.0, dir);
      }
    }
  }
}

void ito_plasma_stepper::compute_ito_diffusion(){
  CH_TIME("ito_plasma_stepper::compute_ito_diffusion()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_diffusion()" << endl;
  }

  EBAMRCellData E;
  m_amr->allocate_ptr(E);
  m_amr->alias(E, m_phase, m_poisson->get_E());

  Vector<EBAMRCellData*> diffusion = m_ito->get_diffusion();
  Vector<EBAMRCellData*> densities = m_ito->get_densities();

  this->compute_ito_diffusion(diffusion, densities, E, m_time);
}

void ito_plasma_stepper::compute_ito_diffusion(Vector<EBAMRCellData*>&       a_diffco,
						const Vector<EBAMRCellData*>& a_densities,
						const EBAMRCellData&          a_E,
						const Real                    a_time){
  CH_TIME("ito_plasma_stepper::compute_ito_diffusion(velo, E, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_diffusion(velo, E, time)" << endl;
  }

  // TLDR: The a_E that comes in is defined on the fluid realm. Copy it to the particle realm.
  m_particle_scratchD.copy(a_E);
  m_amr->average_down(m_particle_scratchD, m_particle_realm, m_phase);
  m_amr->interp_ghost(m_particle_scratchD, m_particle_realm, m_phase);

  //  m_amr->interpolate_to_centroids(m_particle_scratchD, m_fluid_realm, m_phase);

  const int num_ito_species = m_physics->get_num_ito_species();
  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){

    Vector<LevelData<EBCellFAB>* > diffusion(num_ito_species);
    Vector<LevelData<EBCellFAB>* > densities(num_ito_species);
    
    for (int idx = 0; idx < a_diffco.size(); idx++){
      diffusion[idx] = (*a_diffco[idx])[lvl];
      densities[idx] = (*a_densities[idx])[lvl];
    }

    this->compute_ito_diffusion(diffusion, densities, *m_particle_scratchD[lvl], lvl, a_time);
  }

  // Average down, interpolate ghost cells, and then interpolate to particle positions
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<ito_solver>& solver = solver_it();

    if(solver->is_diffusive()){
      m_amr->average_down(*a_diffco[idx], m_particle_realm, m_phase);
      m_amr->interp_ghost(*a_diffco[idx], m_particle_realm, m_phase);

      solver->interpolate_diffusion();
    }
  }
}

void ito_plasma_stepper::compute_ito_diffusion(Vector<LevelData<EBCellFAB>* >&       a_diffco,
					       const Vector<LevelData<EBCellFAB>* >& a_densities,
					       const LevelData<EBCellFAB>&           a_E,
					       const int                             a_level,
					       const Real                            a_time){
  CH_TIME("ito_plasma_stepper::compute_ito_diffusion(velo, E, level, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_diffusion(velo, E, level, time)" << endl;
  }

  const int num_ito_species = m_physics->get_num_ito_species();

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());

    Vector<EBCellFAB*> diffusion(num_ito_species);
    Vector<EBCellFAB*> densities(num_ito_species);;
    
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();

      if(solver_it()->is_diffusive()){
	diffusion[idx] = &(*a_diffco[idx])[dit()];
      }
      densities[idx] = &(*a_densities[idx])[dit()];
    }

    this->compute_ito_diffusion(diffusion, densities, a_E[dit()], a_level, dit(), box, a_time);
  }
}

void ito_plasma_stepper::compute_ito_diffusion(Vector<EBCellFAB*>&       a_diffco,
					       const Vector<EBCellFAB*>& a_densities,
					       const EBCellFAB&          a_E,
					       const int                 a_level,
					       const DataIndex           a_dit,
					       const Box                 a_box,
					       const Real                a_time){
  CH_TIME("ito_plasma_stepper::compute_ito_diffusion(velo, E, level, dit, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_diffusion(velo, E, level, dit, time)" << endl;
  }

  const int comp         = 0;
  const Real dx          = m_amr->get_dx()[a_level];
  const RealVect prob_lo = m_amr->get_prob_lo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();

  // Do regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv   = bit();
    const RealVect pos = m_amr->get_prob_lo() + dx*(RealVect(iv) + 0.5*RealVect::Unit);
    const RealVect e   = RealVect(D_DECL(E(iv,0), E(iv, 1), E(iv, 2)));

    // Make grid densities
    Vector<Real> densities;
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      densities.push_back((*a_densities[idx]).getSingleValuedFAB()(iv, comp));
    }

    // Call ito_physics and compute diffusion for each particle species
    Vector<Real> diffusion = m_physics->compute_ito_diffusion(a_time, pos, e, densities);
    
    // Put diffusion where they belong
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<ito_solver>& solver = solver_it();
      if(solver->is_diffusive()){
	const int idx = solver_it.get_solver();
	(*a_diffco[idx]).getSingleValuedFAB()(iv, comp) = diffusion[idx];
      }
    }
  }

  // Do irregular cells
  VoFIterator& vofit = (*m_amr->get_vofit(m_particle_realm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect e    = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, prob_lo);
    
    // Get densities
    Vector<Real> densities;
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      densities.push_back((*a_densities[idx])(vof, comp));
    }
    
    // Compute diffusion
    Vector<Real> diffusion = m_physics->compute_ito_diffusion(a_time, pos, e, densities);

    // Put diffusion in the appropriate place.
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      if(solver_it()->is_diffusive()){
	const int idx = solver_it.get_solver();
	(*a_diffco[idx])(vof, comp) = diffusion[idx];
      }
    }
  }

  // Covered is bogus.
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    if(solver_it()->is_diffusive()){
      const int idx = solver_it.get_solver();
      a_diffco[idx]->setCoveredCellVal(0.0, comp);
    }
  }
}

void ito_plasma_stepper::advance_reaction_network(const Real a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network(a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network(a_dt)" << endl;
  }
  const int num_ito_species = m_physics->get_num_ito_species();
  const int num_rte_species = m_physics->get_num_rte_species();
  
  Vector<particle_container<ito_particle>* > particles(num_ito_species);  // Current particles. 
  Vector<particle_container<photon>* > bulk_photons(num_rte_species);     // Photons absorbed on mesh
  Vector<particle_container<photon>* > new_photons(num_rte_species);      // Produced photons go here.

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    particles[solver_it.get_solver()] = &(solver_it()->get_particles());
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    bulk_photons[solver_it.get_solver()] = &(solver_it()->get_bulk_photons());
    new_photons[solver_it.get_solver()] = &(solver_it()->get_source_photons());
  }

  // Make alias
  EBAMRCellData E;
  m_amr->allocate_ptr(E);
  m_amr->alias(E, m_phase, m_poisson->get_E());

  this->advance_reaction_network(particles, bulk_photons, new_photons, E, a_dt);
}

void ito_plasma_stepper::advance_reaction_network(Vector<particle_container<ito_particle>* >& a_particles,
						  Vector<particle_container<photon>* >&       a_photons,
						  Vector<particle_container<photon>* >&       a_newPhotons,
						  const EBAMRCellData&                        a_E,
						  const Real                                  a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network(container x3, a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network(container x3, a_dt)" << endl;
  }

  // TLDR: The a_E that comes in is defined on the fluid realm. Copy it to the particle relam
  m_particle_scratchD.copy(a_E);
  m_amr->average_down(m_particle_scratchD, m_particle_realm, m_phase);
  m_amr->interp_ghost(m_particle_scratchD, m_particle_realm, m_phase);

  m_amr->interpolate_to_centroids(m_particle_scratchD, m_particle_realm, m_phase);

  const int num_ito_species = m_physics->get_num_ito_species();
  const int num_rte_species = m_physics->get_num_rte_species();

  Vector<AMRCellParticles<ito_particle>* > particles(num_ito_species);
  Vector<AMRCellParticles<photon>* >       photons(num_ito_species);
  Vector<AMRCellParticles<photon>* >       newPhotons(num_ito_species);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    particles[idx] = &(a_particles[idx]->get_cell_particles());
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    photons[idx]    = &(a_photons[idx]->get_cell_particles());
    newPhotons[idx] = &(a_newPhotons[idx]->get_cell_particles());
  }
				
  //Advance reaction network
  this->advance_reaction_network(particles, photons, newPhotons, m_particle_scratchD, a_dt);

  // Discard shit that is under the PVR
#if 0
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    a_particles[solver_it.get_solver()]->discard_invalid_particles();
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    a_photons[solver_it.get_solver()]->discard_invalid_particles();
    a_newPhotons[solver_it.get_solver()]->discard_invalid_particles();
  }
#endif
}

void ito_plasma_stepper::advance_reaction_network(Vector<AMRCellParticles<ito_particle>* >& a_particles,
						  Vector<AMRCellParticles<photon>* >&       a_photons,
						  Vector<AMRCellParticles<photon>* >&       a_newPhotons,
						  const EBAMRCellData&                      a_E,
						  const Real                                a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network(AMRCellParticles x3, lvl, a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network(AMRCellParticles x3, lvl, a_dt)" << endl;
  }

  const int num_ito_species = m_physics->get_num_ito_species();
  const int num_rte_species = m_physics->get_num_rte_species();

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    Vector<LayoutData<BinFab<ito_particle> >* > particles(num_ito_species);
    Vector<LayoutData<BinFab<photon> >* >       photons(num_rte_species);
    Vector<LayoutData<BinFab<photon> >* >       newPhotons(num_rte_species);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      particles[idx] = &(*(*a_particles[idx])[lvl]);
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      photons[idx]    = &(*(*a_photons[idx])[lvl]);
      newPhotons[idx] = &(*(*a_newPhotons[idx])[lvl]);
    }
    
    this->advance_reaction_network(particles, photons, newPhotons, *a_E[lvl], lvl, a_dt);
  }
}

void ito_plasma_stepper::advance_reaction_network(Vector<LayoutData<BinFab<ito_particle> >* >& a_particles,
						  Vector<LayoutData<BinFab<photon> >* >&       a_photons,
						  Vector<LayoutData<BinFab<photon> >* >&       a_newPhotons,
						  const LevelData<EBCellFAB>&                  a_E,
						  const int                                    a_lvl,
						  const Real                                   a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network(LayoutData<BinFab> x3, lvl, dit, box, a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network(LayoutData<BinFab> x3, lvl, dit, box, a_dt)" << endl;
  }

  const int num_ito_species = m_physics->get_num_ito_species();
  const int num_rte_species = m_physics->get_num_rte_species();

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[a_lvl];
  const Real dx = m_amr->get_dx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());

    Vector<BinFab<ito_particle>* > particles(num_ito_species);
    Vector<BinFab<photon>* > photons(num_rte_species);;
    Vector<BinFab<photon>* > newPhotons(num_rte_species);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      particles[idx] = &((*a_particles[idx])[dit()]);
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      photons[idx]    = &((*a_photons[idx])[dit()]);
      newPhotons[idx] = &((*a_newPhotons[idx])[dit()]);
    }

    this->advance_reaction_network(particles, photons, newPhotons, a_E[dit()], a_lvl, dit(), box, dx, a_dt);
  }
}

void ito_plasma_stepper::advance_reaction_network(Vector<BinFab<ito_particle>* >& a_particles,
						  Vector<BinFab<photon>* >&       a_photons,
						  Vector<BinFab<photon>* >&       a_newPhotons,
						  const EBCellFAB&                a_E,
						  const int                       a_lvl,
						  const DataIndex                 a_dit,
						  const Box                       a_box,
						  const Real                      a_dx,
						  const Real                      a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network(BinFab x3, lvl, dit, box, a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network(BinFab x3, lvl, dit, box, a_dt)" << endl;
  }

  const int comp = 0;

  const int num_ito_species = m_physics->get_num_ito_species();
  const int num_rte_species = m_physics->get_num_rte_species();

  const RealVect prob_lo = m_amr->get_prob_lo();
  const RealVect dx      = a_dx*RealVect::Unit;

  const EBISBox& ebisbox = m_amr->get_ebisl(m_particle_realm, m_phase)[a_lvl][a_dit];
  const EBISBox& ebgraph = m_amr->get_ebisl(m_particle_realm, m_phase)[a_lvl][a_dit];

  const BaseFab<Real>& Efab = a_E.getSingleValuedFAB();

  // Regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    
    if(ebisbox.isRegular(iv)){
      const Real kappa   = 1.0;
      const RealVect pos = prob_lo + a_dx*(RealVect(iv) + 0.5*RealVect::Unit);
      const RealVect e   = RealVect(D_DECL(Efab(iv, 0), Efab(iv, 1), Efab(iv, 2)));
      
      Vector<List<ito_particle>* > particles(num_ito_species);
      Vector<List<photon>* >       photons(num_rte_species);
      Vector<List<photon>* >       newPhotons(num_rte_species);

      for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.get_solver();
      
	List<ito_particle>& bp = (*a_particles[idx])(iv, comp);
	particles[idx] = &bp;
      }

      for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.get_solver();
      
	List<photon>& bp    = (*a_photons[idx])(iv, comp);
	List<photon>& bpNew = (*a_newPhotons[idx])(iv, comp);
      
	photons[idx]    = &bp;
	newPhotons[idx] = &bpNew;
      }

      // Dummy stuff for regular cells
      const RealVect lo = -0.5*RealVect::Unit;
      const RealVect hi =  0.5*RealVect::Unit;
      const RealVect n  = RealVect::Zero;
      const RealVect c  = RealVect::Zero;

      // Update reaction rates
      m_physics->update_reaction_rates(e, a_dx, kappa);
      
      // Advance reactions
      m_physics->advance_reaction_network(particles, photons, newPhotons, e, pos, c, c, n, lo, hi, a_dx, kappa, a_dt);
    }
  }

  // Now do the irregular cells
  VoFIterator& vofit = (*m_amr->get_vofit(m_particle_realm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex vof = vofit();
    const IntVect iv   = vof.gridIndex();
    const RealVect pos = prob_lo + a_dx*(RealVect(iv) + 0.5*RealVect::Unit);
    const RealVect cen = ebisbox.centroid(vof);
    const Real kappa   = ebisbox.volFrac(vof);
    const RealVect e   = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect n   = ebisbox.normal(vof);
    const RealVect ebc = ebisbox.bndryCentroid(vof);


    // Compute a small box that encloses the cut-cell volume
    RealVect lo = -0.5*RealVect::Unit;
    RealVect hi =  0.5*RealVect::Unit;
    if(kappa < 1.0){
      data_ops::compute_min_valid_box(lo, hi, n, ebc);
    }

    Vector<List<ito_particle>* > particles(num_ito_species);
    Vector<List<photon>* >       photons(num_rte_species);
    Vector<List<photon>* >       newPhotons(num_rte_species);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      
      List<ito_particle>& bp = (*a_particles[idx])(iv, comp);
      particles[idx] = &bp;
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      
      List<photon>& bp    = (*a_photons[idx])(iv, comp);
      List<photon>& bpNew = (*a_newPhotons[idx])(iv, comp);
      
      photons[idx]    = &bp;
      newPhotons[idx] = &bpNew;
    }

    // Update reaction rates
    m_physics->update_reaction_rates(e, a_dx, kappa);

    // Advance reactions
    m_physics->advance_reaction_network(particles, photons, newPhotons, e, pos, cen, ebc, n, lo, hi, a_dx, kappa, a_dt);
  }
}

void ito_plasma_stepper::advance_photons(const Real a_dt){
  CH_TIME("ito_plasma_stepper::advance_photons(a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_advance_photons(a_dt)" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<mc_photo>& solver = solver_it();
    
    // Add source photons and move the photons
    particle_container<photon>& photons        = solver->get_photons();
    particle_container<photon>& bulkPhotons    = solver->get_bulk_photons();
    particle_container<photon>& ebPhotons      = solver->get_eb_photons();
    particle_container<photon>& domainPhotons  = solver->get_domain_photons();
    particle_container<photon>& sourcePhotons  = solver->get_source_photons();

    if(solver->is_instantaneous()){
      solver->clear(photons);

      // Add source photons
      photons.add_particles(sourcePhotons);
      solver->clear(sourcePhotons);

      // Instantaneous advance
      solver->advance_photons_stationary(bulkPhotons, ebPhotons, domainPhotons, photons);
    }
    else{
      // Add source photons
      photons.add_particles(sourcePhotons);
      solver->clear(sourcePhotons);

      // Stationary advance
      solver->advance_photons_transient(bulkPhotons, ebPhotons, domainPhotons, photons, a_dt);
    }
  }
}

void ito_plasma_stepper::sort_photons_by_cell(){
  CH_TIME("ito_plasma_stepper::sort_photons_by_cell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_photons_by_cell()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_photons_by_cell();
  }
}

void ito_plasma_stepper::sort_photons_by_patch(){
  CH_TIME("ito_plasma_stepper::sort_photons_by_patch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_photons_by_patch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_photons_by_patch();
  }
}

void ito_plasma_stepper::sort_source_photons_by_cell(){
  CH_TIME("ito_plasma_stepper::sort_source_photons_by_cell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_source_photons_by_cell()" << endl;
  }
  
  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_source_photons_by_cell();
  }
}

void ito_plasma_stepper::sort_source_photons_by_patch(){
  CH_TIME("ito_plasma_stepper::sort_source_photons_by_patch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_source_photons_by_patch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_source_photons_by_patch();
  }
}

void ito_plasma_stepper::sort_bulk_photons_by_cell(){
  CH_TIME("ito_plasma_stepper::sort_bulk_photons_by_cell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_bulk_photons_by_cell()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_bulk_photons_by_cell();
  }
}

void ito_plasma_stepper::sort_bulk_photons_by_patch(){
  CH_TIME("ito_plasma_stepper::sort_bulk_photons_by_patch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_bulk_photons_by_patch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_bulk_photons_by_patch();
  }
}

void ito_plasma_stepper::sort_eb_photons_by_cell(){
  CH_TIME("ito_plasma_stepper::sort_eb_photons_by_cell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_eb_photons_by_cell()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_eb_photons_by_cell();
  }
}

void ito_plasma_stepper::sort_eb_photons_by_patch(){
  CH_TIME("ito_plasma_stepper::sort_eb_photons_by_patch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_eb_photons_by_patch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_eb_photons_by_patch();
  }
}

void ito_plasma_stepper::sort_domain_photons_by_cell(){
  CH_TIME("ito_plasma_stepper::sort_domain_photons_by_cell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_domain_photons_by_cell()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_domain_photons_by_cell();
  }
}

void ito_plasma_stepper::sort_domain_photons_by_patch(){
  CH_TIME("ito_plasma_stepper::sort_domain_photons_by_patch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_domain_photons_by_patch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_domain_photons_by_patch();
  }

}

bool ito_plasma_stepper::load_balance(Vector<Vector<int> >&            a_procs,
					   Vector<Vector<Box> >&            a_boxes,
					   const std::string                a_realm,
					   const Vector<DisjointBoxLayout>& a_grids,
					   const int                        a_lmin,
					   const int                        a_finest_level){
  CH_TIME("ito_plasma_stepper::load_balance");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper_stepper::load_balance" << endl;
  }

  bool ret;

  if(!m_load_balance){
    ret = false;
  }
  else{
    if(a_realm == m_particle_realm){
      ret = this->load_balance_particle_realm(a_procs, a_boxes, a_realm, a_grids, a_lmin, a_finest_level);
    }
    else{
      ret = false;
    }
  }

  return ret;
}

bool ito_plasma_stepper::load_balance_particle_realm(Vector<Vector<int> >&            a_procs,
						     Vector<Vector<Box> >&            a_boxes,
						     const std::string                a_realm,
						     const Vector<DisjointBoxLayout>& a_grids,
						     const int                        a_lmin,
						     const int                        a_finest_level){
  
  bool ret = false;
  
  if(m_load_balance){
    RefCountedPtr<ito_solver>& solver           = m_ito->get_solvers()[0];
    particle_container<ito_particle>& particles = solver->get_particles();
  
    particles.regrid(a_grids, m_amr->get_domains(), m_amr->get_dx(), m_amr->get_ref_rat(), a_lmin, a_finest_level);

    a_procs.resize(1 + a_finest_level);
    a_boxes.resize(1 + a_finest_level);
  
    // Compute loads on each level
    for (int lvl = 0; lvl < a_lmin; lvl++){
      a_procs[lvl] = a_grids[lvl].procIDs();
      a_boxes[lvl] = a_grids[lvl].boxArray();
    }

    for (int lvl = a_lmin; lvl <= a_finest_level; lvl++){
      Vector<long int> loads;
      a_boxes[lvl] = a_grids[lvl].boxArray();
    
      solver->compute_loads(loads, a_grids[lvl], lvl);

#ifdef CH_MPI
      int count = loads.size();
      Vector<long int> tmp(count);
      MPI_Allreduce(&(loads[0]),&(tmp[0]), count, MPI_LONG, MPI_SUM, Chombo_MPI::comm);
      loads = tmp;
#endif

      LoadBalance(a_procs[lvl], loads, a_boxes[lvl]);
    }

    // Put particles back
    particles.pre_regrid(a_lmin);

    ret = true;
  }

  return ret;
}
