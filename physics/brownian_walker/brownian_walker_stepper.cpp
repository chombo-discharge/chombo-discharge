/*!
  @file   brownian_walker_stepper.cpp
  @brief  Implementation of brownian_walker_stepper.H
  @author Robert Marskar
  @data   March 2020
*/

#include "brownian_walker_stepper.H"
#include "brownian_walker_species.H"
#include "poly.H"

#include <ParmParse.H>
#include <PolyGeom.H>
#include <BinFab.H>

#include "CD_NamespaceHeader.H"
using namespace physics::brownian_walker;

brownian_walker_stepper::brownian_walker_stepper(){
  ParmParse pp("brownian_walker");

  m_phase = phase::gas;

  pp.get("realm",          m_realm);
  pp.get("diffco",         m_diffco);
  pp.get("omega",          m_omega);
  pp.get("verbosity",      m_verbosity);
  pp.get("ppc",            m_ppc);
  pp.get("max_cells_hop",  m_max_cells_hop);
  pp.get("load_balance",   m_load_balance);
}

brownian_walker_stepper::brownian_walker_stepper(RefCountedPtr<ito_solver>& a_solver) : brownian_walker_stepper() {
  m_solver = a_solver;
}

brownian_walker_stepper::~brownian_walker_stepper(){

}

void brownian_walker_stepper::parse_runtime_options() {

  ParmParse pp("brownian_walker");
  
  pp.get("verbosity",      m_verbosity);
  pp.get("ppc",            m_ppc);
  pp.get("max_cells_hop",  m_max_cells_hop);
  pp.get("load_balance",   m_load_balance);
  
  m_solver->parse_runtime_options();
}

void brownian_walker_stepper::initial_data(){
  CH_TIME("brownian_walker_stepper::initial_data");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::initial_data" << endl;
  }
  
  m_solver->initial_data();
  if(m_ppc > 0){
    m_solver->sort_particles_by_cell(ito_solver::which_container::bulk);
    m_solver->make_superparticles(ito_solver::which_container::bulk,m_ppc);
    m_solver->sort_particles_by_patch(ito_solver::which_container::bulk);
  }

  if(m_solver->is_diffusive()){
    m_solver->set_diffco_func(m_diffco);
  }
  if(m_solver->is_mobile()){
    this->set_velocity();
  }
}

void brownian_walker_stepper::post_initialize(){
  CH_TIME("brownian_walker_stepper::post_initialize");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::post_initialize" << endl;
  }
}

void brownian_walker_stepper::set_velocity(){
  CH_TIME("brownian_walker_stepper::set_velocity");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::set_velocity" << endl;
  }
  m_solver->set_mobility(1.0);
  
  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    this->set_velocity(lvl);
  }

  EBAMRCellData& vel = m_solver->get_velo_func();
  m_amr->averageDown(vel, m_realm, m_phase);
  m_amr->interpGhost(vel, m_realm, m_phase);
}

bool brownian_walker_stepper::load_balance_realm(const std::string a_realm) const {
  CH_TIME("brownian_walker_stepper::load_balance_realm");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::load_balance_realm" << endl;
  }

  bool ret = false;

  if(m_load_balance && a_realm == m_realm){
    ret = true;
  }

  return ret;
}

void brownian_walker_stepper::load_balance_boxes(Vector<Vector<int> >&            a_procs,
						 Vector<Vector<Box> >&            a_boxes,
						 const std::string                a_realm,
						 const Vector<DisjointBoxLayout>& a_grids,
						 const int                        a_lmin,
						 const int                        a_finest_level){
  CH_TIME("brownian_walker_stepper::load_balance_boxes");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::load_balance_boxes" << endl;
  }
  
  if(m_load_balance && a_realm == m_realm){
    particle_container<ito_particle>& particles = m_solver->get_particles(ito_solver::which_container::bulk);
  
    particles.regrid(a_grids, m_amr->getDomains(), m_amr->getDx(), m_amr->getRefinementRatios(), a_lmin, a_finest_level);

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
    
      m_solver->compute_loads(loads, a_grids[lvl], lvl);

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
  }
  else{
    MayDay::Abort("brownian_walker_stepper::load_balance_boxes - logic bust");
  }
}

void brownian_walker_stepper::set_velocity(const int a_level){
  CH_TIME("brownian_walker_stepper::set_velocity(level)");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::set_velocity(level)" << endl;
  }

  // TLDR: This code goes down to each cell on grid level a_level and sets the velocity to omega*r
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box = dbl.get(dit());

    EBCellFAB& vel = (*(m_solver->get_velo_func())[a_level])[dit()];
    BaseFab<Real>& vel_reg = vel.getSingleValuedFAB();

    vel.setVal(0.0);
      
    // Regular cells
    for (BoxIterator bit(box); bit.ok(); ++bit){
      const IntVect iv = bit();
      const RealVect pos = m_amr->getProbLo() + (RealVect(iv) + 0.5*RealVect::Unit)*m_amr->getDx()[a_level];

      const Real r     = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
      const Real theta = atan2(pos[1],pos[0]);

      vel_reg(iv,0) = -r*m_omega*sin(theta);
      vel_reg(iv,1) =  r*m_omega*cos(theta);
    }

    // Irregular and multicells
    const EBISBox& ebisbox = vel.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_level])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){

      const VolIndex vof = vofit();
      const IntVect iv   = vof.gridIndex();
      const RealVect pos = m_amr->getProbLo() + (RealVect(iv) + 0.5*RealVect::Unit)*m_amr->getDx()[a_level];

      const Real r     = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
      const Real theta = atan2(pos[1],pos[0]);

      vel(vof,0) = -r*m_omega*sin(theta);
      vel(vof,1) =  r*m_omega*cos(theta);
    }

    // Now set the mobility for all the particles
    List<ito_particle>& particles = m_solver->get_particles(ito_solver::which_container::bulk)[a_level][dit()].listItems();
    for (ListIterator<ito_particle> lit(particles); lit.ok(); ++lit){
      lit().mobility() = 1.0;
    }
  }
}

void brownian_walker_stepper::write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const{
  CH_TIME("brownian_walker_stepper::write_checkpoint_data");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::write_checkpoint_data" << endl;
  }

  m_solver->write_checkpoint_level(a_handle, a_lvl);
}

void brownian_walker_stepper::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) {
  CH_TIME("brownian_walker_stepper::read_checkpoint_data");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::read_checkpoint_data" << endl;
  }
  
  m_solver->read_checkpoint_level(a_handle, a_lvl);
}

void brownian_walker_stepper::post_checkpoint_setup() {
  CH_TIME("brownian_walker_stepper::post_checkpoint_setup");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::post_checkpoint_setup" << endl;
  }

  m_solver->remap();
  if(m_ppc > 0){
    m_solver->sort_particles_by_cell(ito_solver::which_container::bulk);
    m_solver->make_superparticles(ito_solver::which_container::bulk, m_ppc);
    m_solver->sort_particles_by_patch(ito_solver::which_container::bulk);
  }
  m_solver->deposit_particles();

  if(m_solver->is_diffusive()){
    m_solver->set_diffco_func(m_diffco);
  }
  if(m_solver->is_mobile()){
    this->set_velocity();
  }
}

int brownian_walker_stepper::get_num_plot_vars() const {
  CH_TIME("brownian_walker_stepper::get_num_plot_vars");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::get_num_plot_vars" << endl;
  }

  return m_solver->get_num_plotvars();
}

void brownian_walker_stepper::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const {
  CH_TIME("brownian_walker_stepper::write_plot_data");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::write_plot_data" << endl;
  }
  a_plotvar_names.append(m_solver->get_plotvar_names());
  m_solver->write_plot_data(a_output, a_icomp);
}

void brownian_walker_stepper::compute_dt(Real& a_dt, time_code& a_timecode) {
  CH_TIME("brownian_walker_stepper::compute_dt");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::compute_dt" << endl;
  }

  m_solver->set_mobility(1.0); 
  m_solver->interpolate_velocities();
  m_solver->interpolate_diffusion();

  a_dt = m_max_cells_hop*m_solver->compute_dt();
}

void brownian_walker_stepper::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("brownian_walker_stepper::synchronize_solver_times");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::synchronize_solver_times" << endl;
  }
  
  m_solver->set_time(a_step, a_time, a_dt);

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void brownian_walker_stepper::print_step_report() {
  CH_TIME("brownian_walker_stepper::print_step_report");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::print_step_report" << endl;
  }

  // Do nothing
  const size_t local_particles  = m_solver->get_num_particles(ito_solver::which_container::bulk, true);
  const size_t global_particles = m_solver->get_num_particles(ito_solver::which_container::bulk, false);

  pout() << "                                   #part = " << local_particles << " (" << global_particles << ")" << endl;
  
}

bool brownian_walker_stepper::need_to_regrid() {
  return false;
}

void brownian_walker_stepper::pre_regrid(const int a_lbase, const int a_old_finest_level){
  CH_TIME("brownian_walker_stepper::pre_regrid");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::pre_regrid" << endl;
  }

  // TLDR: base is the finest level that DOES NOT CHANGE
  const int base = 0;
  const int finest_level = m_amr->getFinestLevel();
  
  m_solver->pre_regrid(base, finest_level);
}

void brownian_walker_stepper::deallocate() {
  CH_TIME("brownian_walker_stepper::deallocate");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::deallocate" << endl;
  }
}

void brownian_walker_stepper::setup_solvers() {
  CH_TIME("brownian_walker_stepper::setup_solvers");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::setup_solvers" << endl;
  }

  m_species = RefCountedPtr<ito_species> (new brownian_walker_species());

  m_solver->set_verbosity(m_verbosity);
  m_solver->parse_options();
  m_solver->set_amr(m_amr);
  m_solver->set_species(m_species);
  m_solver->set_phase(m_phase);
  m_solver->set_computational_geometry(m_compgeom);
  m_solver->set_realm(m_realm);
}

void brownian_walker_stepper::registerRealms() {
  CH_TIME("brownian_walker_stepper::registerRealms");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::registerRealms" << endl;
  }

  m_amr->registerRealm(m_realm);
}

void brownian_walker_stepper::registerOperators() {
  CH_TIME("brownian_walker_stepper::registerOperators");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::registerOperators" << endl;
  }

  m_solver->registerOperators();
}

void brownian_walker_stepper::allocate() {
  m_solver->allocate_internals(); // Allocate some internal storage
}

Real brownian_walker_stepper::advance(const Real a_dt) {
  CH_TIME("brownian_walker_stepper::advance");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::advance" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();
  const RealVect origin  = m_amr->getProbLo();
  
  //  m_solver->remap();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const RealVect dx                     = m_amr->getDx()[lvl]*RealVect::Unit;
    const DisjointBoxLayout& dbl          = m_amr->getGrids(m_realm)[lvl];
    ParticleData<ito_particle>& particles = m_solver->get_particles(ito_solver::which_container::bulk)[lvl];

    const EBISLayout& ebisl = m_amr->getEBISLayout(m_realm, m_solver->get_phase())[lvl];

    if(m_solver->is_mobile() || m_solver->is_diffusive()){
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	// Create a copy. 
	List<ito_particle>& particleList = particles[dit()].listItems();
	List<ito_particle>  particleCopy = List<ito_particle>(particleList);

	// The list iterator is NOT an indexing iterator but iterates over the list given
	// in the constructor. So, we need one for velocities and one for the copy
	ListIterator<ito_particle> lit(particleList);
	ListIterator<ito_particle> litC(particleCopy);

	// Half Euler step and evaluate velocity at half step
	if(m_solver->is_mobile()){
	  for (lit.rewind(); lit; ++lit){ 
	    ito_particle& p = particleList[lit];
	    p.position() += 0.5*p.velocity()*a_dt;
	  }

	  m_solver->interpolate_velocities(lvl, dit());

	  // Final stage
	  for (lit.rewind(), litC.rewind(); lit, litC; ++lit, ++litC){
	    ito_particle& p    = particleList[lit];
	    ito_particle& oldP = particleCopy[litC];
	    p.position() = oldP.position() + p.velocity()*a_dt;
	  }
	}

	// Add a diffusion hop
	if(m_solver->is_diffusive()){
	  for (lit.rewind(); lit; ++lit){ 
	    ito_particle& p = particleList[lit];
	    const RealVect ran = m_solver->random_gaussian();
	    const RealVect hop = ran*sqrt(2.0*p.diffusion()*a_dt);
	    p.position() += hop;
	  }
	}

	// Do particle bounceback on the EB
	const EBISBox& ebisbox = ebisl[dit()];
	if(!ebisbox.isAllRegular() || !ebisbox.isAllCovered()){ 
	  
	  // This is the implicit function
	  RefCountedPtr<BaseIF> func;
	  if(m_solver->get_phase() == phase::gas){
	    func = m_compgeom->get_gas_if();
	  }
	  else {
	    func = m_compgeom->get_sol_if();
	  }

	  for (lit.rewind(), litC.rewind(); lit, litC; ++lit, ++litC){
	    ito_particle& p    = particleList[lit];
	    ito_particle& oldP = particleCopy[litC];

	    const RealVect& oldPos = oldP.position();
	    const RealVect& newPos = p.position();

	    const Real fOld = func->value(oldPos);
	    const Real fNew = func->value(newPos);

	    // If the particle crossed the boundary, 
	    if(fOld*fNew <= 0.0){
	      const RealVect xb = poly::brent_root_finder(func, oldPos, newPos);
	      const IntVect iv = locateBin(xb, dx, origin);
	      const VolIndex vof(iv, 0);

	      // Do a dummy check
	      //	      const Box b(iv-IntVect::Unit, iv+IntVect::Unit);
	      if(!ebisbox.isIrregular(iv)){
		MayDay::Abort("brownian_walker_stepper::advance - logic bust in EB bounce-back code");
	      }
	      else{
		const RealVect n     = ebisbox.normal(vof);
		const RealVect bback = 2.0*n*PolyGeom::dot((newPos - xb), n);
		p.position() -= bback;
	      }
	    }
	  }
	}
      }
    }
  }

  // Remap and deposit particles
  m_solver->remap();
  m_solver->deposit_particles();

  return a_dt;
}

void brownian_walker_stepper::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) {
  CH_TIME("brownian_walker_stepper::regrid");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::regrid" << endl;
  }

  m_solver->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  if(m_solver->is_diffusive()){
    m_solver->set_diffco_func(m_diffco);
  }
  if(m_solver->is_mobile()){
    this->set_velocity();
  }

  if(m_ppc > 0){
    m_solver->sort_particles_by_cell(ito_solver::which_container::bulk);
    m_solver->make_superparticles(ito_solver::which_container::bulk, m_ppc);
    m_solver->sort_particles_by_patch(ito_solver::which_container::bulk);
    
    m_solver->set_mobility(1.0); // Superparticle algorithm only conserves mass, energy. Diffusion and mobility needs to be reset.
  }
}

void brownian_walker_stepper::post_regrid(){
  CH_TIME("brownian_walker_stepper::post_regrid");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::post_regrid" << endl;
  }
}
#include "CD_NamespaceFooter.H"
