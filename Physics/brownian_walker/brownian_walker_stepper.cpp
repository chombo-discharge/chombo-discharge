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

  pp.get("Realm",          m_realm);
  pp.get("diffco",         m_faceCenteredDiffusionCoefficient);
  pp.get("omega",          m_omega);
  pp.get("verbosity",      m_verbosity);
  pp.get("ppc",            m_ppc);
  pp.get("max_cells_hop",  m_max_cells_hop);
  pp.get("LoadBalancing",   m_LoadBalancing);
}

brownian_walker_stepper::brownian_walker_stepper(RefCountedPtr<ito_solver>& a_solver) : brownian_walker_stepper() {
  m_solver = a_solver;
}

brownian_walker_stepper::~brownian_walker_stepper(){

}

void brownian_walker_stepper::parseRuntimeOptions() {

  ParmParse pp("brownian_walker");
  
  pp.get("verbosity",      m_verbosity);
  pp.get("ppc",            m_ppc);
  pp.get("max_cells_hop",  m_max_cells_hop);
  pp.get("LoadBalancing",   m_LoadBalancing);
  
  m_solver->parseRuntimeOptions();
}

void brownian_walker_stepper::initialData(){
  CH_TIME("brownian_walker_stepper::initialData");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::initialData" << endl;
  }
  
  m_solver->initialData();
  if(m_ppc > 0){
    m_solver->sortParticlesByCell(ito_solver::which_container::bulk);
    m_solver->make_superparticles(ito_solver::which_container::bulk,m_ppc);
    m_solver->sortParticlesByPatch(ito_solver::which_container::bulk);
  }

  if(m_solver->isDiffusive()){
    m_solver->setDiffusionCoefficient_func(m_faceCenteredDiffusionCoefficient);
  }
  if(m_solver->isMobile()){
    this->setVelocity();
  }
}

void brownian_walker_stepper::postInitialize(){
  CH_TIME("brownian_walker_stepper::postInitialize");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::postInitialize" << endl;
  }
}

void brownian_walker_stepper::setVelocity(){
  CH_TIME("brownian_walker_stepper::setVelocity");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::setVelocity" << endl;
  }
  m_solver->set_mobility(1.0);
  
  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    this->setVelocity(lvl);
  }

  EBAMRCellData& vel = m_solver->get_velo_func();
  m_amr->averageDown(vel, m_realm, m_phase);
  m_amr->interpGhost(vel, m_realm, m_phase);
}

bool brownian_walker_stepper::loadBalanceThisRealm(const std::string a_realm) const {
  CH_TIME("brownian_walker_stepper::loadBalanceThisRealm");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::loadBalanceThisRealm" << endl;
  }

  bool ret = false;

  if(m_LoadBalancing && a_realm == m_realm){
    ret = true;
  }

  return ret;
}

void brownian_walker_stepper::loadBalanceBoxes(Vector<Vector<int> >&            a_procs,
						 Vector<Vector<Box> >&            a_boxes,
						 const std::string                a_realm,
						 const Vector<DisjointBoxLayout>& a_grids,
						 const int                        a_lmin,
						 const int                        a_finestLevel){
  CH_TIME("brownian_walker_stepper::loadBalanceBoxes");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::loadBalanceBoxes" << endl;
  }
  
  if(m_LoadBalancing && a_realm == m_realm){
    ParticleContainer<ito_particle>& particles = m_solver->getParticles(ito_solver::which_container::bulk);
  
    particles.regrid(a_grids, m_amr->getDomains(), m_amr->getDx(), m_amr->getRefinementRatios(), a_lmin, a_finestLevel);

    a_procs.resize(1 + a_finestLevel);
    a_boxes.resize(1 + a_finestLevel);
  
    // Compute loads on each level
    for (int lvl = 0; lvl < a_lmin; lvl++){
      a_procs[lvl] = a_grids[lvl].procIDs();
      a_boxes[lvl] = a_grids[lvl].boxArray();
    }

    for (int lvl = a_lmin; lvl <= a_finestLevel; lvl++){
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
    particles.preRegrid(a_lmin);
  }
  else{
    MayDay::Abort("brownian_walker_stepper::loadBalanceBoxes - logic bust");
  }
}

void brownian_walker_stepper::setVelocity(const int a_level){
  CH_TIME("brownian_walker_stepper::setVelocity(level)");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::setVelocity(level)" << endl;
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
    List<ito_particle>& particles = m_solver->getParticles(ito_solver::which_container::bulk)[a_level][dit()].listItems();
    for (ListIterator<ito_particle> lit(particles); lit.ok(); ++lit){
      lit().mobility() = 1.0;
    }
  }
}

void brownian_walker_stepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const{
  CH_TIME("brownian_walker_stepper::writeCheckpointData");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::writeCheckpointData" << endl;
  }

  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}

void brownian_walker_stepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl) {
  CH_TIME("brownian_walker_stepper::readCheckpointData");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::readCheckpointData" << endl;
  }
  
  m_solver->readCheckpointLevel(a_handle, a_lvl);
}

void brownian_walker_stepper::postCheckpointSetup() {
  CH_TIME("brownian_walker_stepper::postCheckpointSetup");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::postCheckpointSetup" << endl;
  }

  m_solver->remap();
  if(m_ppc > 0){
    m_solver->sortParticlesByCell(ito_solver::which_container::bulk);
    m_solver->make_superparticles(ito_solver::which_container::bulk, m_ppc);
    m_solver->sortParticlesByPatch(ito_solver::which_container::bulk);
  }
  m_solver->deposit_particles();

  if(m_solver->isDiffusive()){
    m_solver->setDiffusionCoefficient_func(m_faceCenteredDiffusionCoefficient);
  }
  if(m_solver->isMobile()){
    this->setVelocity();
  }
}

int brownian_walker_stepper::getNumberOfPlotVariables() const {
  CH_TIME("brownian_walker_stepper::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::getNumberOfPlotVariables" << endl;
  }

  return m_solver->getNumberOfPlotVariables();
}

void brownian_walker_stepper::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const {
  CH_TIME("brownian_walker_stepper::writePlotData");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::writePlotData" << endl;
  }
  a_plotVariableNames.append(m_solver->getPlotVariableNames());
  m_solver->writePlotData(a_output, a_icomp);
}

void brownian_walker_stepper::computeDt(Real& a_dt, TimeCode& a_timeCode) {
  CH_TIME("brownian_walker_stepper::computeDt");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::computeDt" << endl;
  }

  m_solver->set_mobility(1.0); 
  m_solver->interpolate_velocities();
  m_solver->interpolate_diffusion();

  a_dt = m_max_cells_hop*m_solver->computeDt();
}

void brownian_walker_stepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("brownian_walker_stepper::synchronizeSolverTimes");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::synchronizeSolverTimes" << endl;
  }
  
  m_solver->setTime(a_step, a_time, a_dt);

  m_timeStep = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void brownian_walker_stepper::printStepReport() {
  CH_TIME("brownian_walker_stepper::printStepReport");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::printStepReport" << endl;
  }

  // Do nothing
  const size_t local_particles  = m_solver->get_num_particles(ito_solver::which_container::bulk, true);
  const size_t global_particles = m_solver->get_num_particles(ito_solver::which_container::bulk, false);

  pout() << "                                   #part = " << local_particles << " (" << global_particles << ")" << endl;
  
}

bool brownian_walker_stepper::needToRegrid() {
  return false;
}

void brownian_walker_stepper::preRegrid(const int a_lbase, const int a_oldFinestLevel){
  CH_TIME("brownian_walker_stepper::preRegrid");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::preRegrid" << endl;
  }

  // TLDR: base is the finest level that DOES NOT CHANGE
  const int base = 0;
  const int finest_level = m_amr->getFinestLevel();
  
  m_solver->preRegrid(base, finest_level);
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

  m_solver->setVerbosity(m_verbosity);
  m_solver->parseOptions();
  m_solver->setAmr(m_amr);
  m_solver->setSpecies(m_species);
  m_solver->setPhase(m_phase);
  m_solver->setComputationalGeometry(m_computationalGeometry);
  m_solver->setRealm(m_realm);
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
  m_solver->allocateInternals(); // Allocate some internal storage
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
    ParticleData<ito_particle>& particles = m_solver->getParticles(ito_solver::which_container::bulk)[lvl];

    const EBISLayout& ebisl = m_amr->getEBISLayout(m_realm, m_solver->get_phase())[lvl];

    if(m_solver->isMobile() || m_solver->isDiffusive()){
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	// Create a copy. 
	List<ito_particle>& particleList = particles[dit()].listItems();
	List<ito_particle>  particleCopy = List<ito_particle>(particleList);

	// The list iterator is NOT an indexing iterator but iterates over the list given
	// in the constructor. So, we need one for velocities and one for the copy
	ListIterator<ito_particle> lit(particleList);
	ListIterator<ito_particle> litC(particleCopy);

	// Half Euler step and evaluate velocity at half step
	if(m_solver->isMobile()){
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
	if(m_solver->isDiffusive()){
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
	    func = m_computationalGeometry->get_gas_if();
	  }
	  else {
	    func = m_computationalGeometry->get_sol_if();
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

void brownian_walker_stepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) {
  CH_TIME("brownian_walker_stepper::regrid");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::regrid" << endl;
  }

  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  if(m_solver->isDiffusive()){
    m_solver->setDiffusionCoefficient_func(m_faceCenteredDiffusionCoefficient);
  }
  if(m_solver->isMobile()){
    this->setVelocity();
  }

  if(m_ppc > 0){
    m_solver->sortParticlesByCell(ito_solver::which_container::bulk);
    m_solver->make_superparticles(ito_solver::which_container::bulk, m_ppc);
    m_solver->sortParticlesByPatch(ito_solver::which_container::bulk);
    
    m_solver->set_mobility(1.0); // Superparticle algorithm only conserves mass, energy. Diffusion and mobility needs to be reset.
  }
}

void brownian_walker_stepper::postRegrid(){
  CH_TIME("brownian_walker_stepper::postRegrid");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::postRegrid" << endl;
  }
}
#include "CD_NamespaceFooter.H"
