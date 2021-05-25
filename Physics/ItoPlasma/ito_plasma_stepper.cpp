/*!
  @file   ito_plasmpa_stepper.cpp
  @brief  Implementation of ito_plasma_stepper.H
  @author Robert Marskar
  @date   May 2020
*/

#include "ito_plasma_stepper.H"
#include <CD_DataOps.H>
#include <CD_Units.H>
#include "CD_FieldSolverMultigrid.H"

#include <EBArith.H>
#include <PolyGeom.H>
#include <ParmParse.H>

#include "CD_NamespaceHeader.H"
using namespace Physics::ItoPlasma;

ito_plasma_stepper::ito_plasma_stepper(){
  m_verbosity = -1;
  m_name      = "ito_plasma_stepper";
  m_phase     = phase::gas;

  m_dt   = 0.0;
  m_time = 0.0;
  
  m_halo_buffer = 1;
  m_pvr_buffer  = 0;
  m_load_ppc    = 1.0;

  m_nwo_reactions         = true;
  m_regrid_superparticles = true;

  m_fluid_Realm    = Realm::Primal;
  m_particleRealm = Realm::Primal;
}

ito_plasma_stepper::ito_plasma_stepper(RefCountedPtr<ItoPlasmaPhysics>& a_physics) : ito_plasma_stepper(){
  m_physics   = a_physics;
}

ito_plasma_stepper::~ito_plasma_stepper(){
}

void ito_plasma_stepper::setVerbosity(const int a_verbosity){
  CH_TIME("ito_plasma_stepper::setVerbosity");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setVerbosity" << endl;
  }
  m_verbosity = a_verbosity;
}

void ito_plasma_stepper::setupSolvers(){
  CH_TIME("ito_plasma_stepper::setup_solver");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setupSolvers" << endl;
  }

  // Parse class options
  this->parseOptions();

  // Set up solvers
  this->setup_ito();
  this->setupPoisson();
  this->setupRadiativeTransfer();
  this->setupSigma();

  // Allocate internal stuff
  this->allocateInternals();

  // Set the particle buffers for the Ito solver
  this->set_particle_buffers();
}

void ito_plasma_stepper::setup_ito(){
  CH_TIME("ito_plasma_stepper::setup_ito");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_ito" << endl;
  }

  m_ito->setVerbosity(m_verbosity);
  m_ito->parseOptions();
  m_ito->setAmr(m_amr);
  m_ito->setPhase(m_phase);
  m_ito->setComputationalGeometry(m_computationalGeometry);
  m_ito->setRealm(m_particleRealm);
}

void ito_plasma_stepper::setupPoisson(){
  CH_TIME("ito_plasma_stepper::setupPoisson");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setupPoisson" << endl;
  }

  m_fieldSolver->setVerbosity(m_verbosity);
  m_fieldSolver->parseOptions();
  m_fieldSolver->setAmr(m_amr);
  m_fieldSolver->setComputationalGeometry(m_computationalGeometry);
  m_fieldSolver->setVoltage(m_potential); 
  m_fieldSolver->setRealm(m_fluid_Realm);
}

void ito_plasma_stepper::setupRadiativeTransfer(){
  CH_TIME("ito_plasma_stepper::setupRadiativeTransfer");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setupRadiativeTransfer" << endl;
  }

  m_rte->setVerbosity(m_verbosity);
  m_rte->parseOptions();
  m_rte->setPhase(m_phase);
  m_rte->setAmr(m_amr);
  m_rte->setComputationalGeometry(m_computationalGeometry);
  m_rte->setRealm(m_particleRealm);
  m_rte->sanityCheck();
}

void ito_plasma_stepper::setupSigma(){
  CH_TIME("ito_plasma_stepper::setupSigma");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setupSigma" << endl;
  }

  m_sigma = RefCountedPtr<SigmaSolver> (new SigmaSolver());
  m_sigma->setAmr(m_amr);
  m_sigma->setVerbosity(m_verbosity);
  m_sigma->setComputationalGeometry(m_computationalGeometry);
  m_sigma->setRealm(m_fluid_Realm);
}

void ito_plasma_stepper::set_particle_buffers(){
  CH_TIME("ito_plasma_stepper::set_particle_buffers");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_particle_buffers" << endl;
  }
  
  m_ito->setHalobuffer(m_halo_buffer);
  m_ito->setPVRBuffer(m_pvr_buffer);

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<McPhoto>& solver = solver_it();

    solver->setHalobuffer(m_halo_buffer);
    solver->setPVRBuffer(m_pvr_buffer);
  }
}

void ito_plasma_stepper::allocate() {
  CH_TIME("ito_plasma_stepper::allocate");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::allocate" << endl;
  }

  m_ito->allocateInternals();
  m_rte->allocateInternals();
  m_fieldSolver->allocateInternals();
  m_sigma->allocateInternals();
}

void ito_plasma_stepper::postInitialize(){

}

void ito_plasma_stepper::initialData(){
  CH_TIME("ito_plasma_stepper::initialData");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::initialData" << endl;
  }

  m_ito->initialData(); // This deposits, of course. 
  m_rte->initialData();
  this->initialSigma();

  m_ito->sortParticlesByCell( ItoSolver::WhichContainer::bulk);
  m_ito->makeSuperparticles(    ItoSolver::WhichContainer::bulk, m_ppc);
  m_ito->sortParticlesByPatch(ItoSolver::WhichContainer::bulk);
  
  // Solve Poisson equation and compute the E-field
  this->solvePoisson();

  // Fill solvers with velocities and diffusion coefficients
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
}

void ito_plasma_stepper::initialSigma(){
  CH_TIME("ito_plasma_stepper::initialSigma");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::initialSigma" << endl;
  }

  const RealVect origin  = m_amr->getProbLo();
  const int finest_level = m_amr->getFinestLevel();

  EBAMRIVData& sigma = m_sigma->getPhi();
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluid_Realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_fluid_Realm, phase::gas)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      BaseIVFAB<Real>& state = (*sigma[lvl])[dit()];

      const EBISBox& ebisbox = ebisl[dit()];
      const IntVectSet& ivs  = state.getIVS();
      const EBGraph& ebgraph = state.getEBGraph();
      
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = origin + vof.gridIndex()*dx + 0.5*ebisbox.bndryCentroid(vof)*dx;
	
	for (int comp = 0; comp < state.nComp(); comp++){
	  state(vof, comp) = m_physics->initialSigma(m_time, pos);
	}
      }
    }
  }

  m_amr->averageDown(sigma, m_fluid_Realm, phase::gas);
  m_sigma->resetCells(sigma);
}

void ito_plasma_stepper::postCheckpointSetup(){
  CH_TIME("ito_plasma_stepper::postCheckpointSetup");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::postCheckpointSetup" << endl;
  }

  //this->solvePoisson();
  this->allocateInternals();

  m_ito->remap();
  
  this->post_checkpoint_poisson();
  
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
}

void ito_plasma_stepper::post_checkpoint_poisson(){
  CH_TIME("ito_plasma_stepper::post_checkpoint_poisson");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::post_checkpoint_poisson" << endl;
  }

  // Do some post checkpointing stuff. 
  m_fieldSolver->postCheckpoint();
  
  // Do ghost cells and then compute E
  MFAMRCellData& state = m_fieldSolver->getPotential();

  m_amr->averageDown(state, m_fluid_Realm);
  m_amr->interpGhost(state, m_fluid_Realm);

  m_fieldSolver->computeElectricField();    // Solver checkpoints the potential. Now compute the field.

  // Interpolate the fields to centroids
  EBAMRCellData E;
  m_amr->allocatePointer(E);
  m_amr->alias(E, m_phase, m_fieldSolver->getElectricField());
  
  // Fluid Realm
  m_fluid_E.copy(E);
  m_amr->averageDown(m_fluid_E,             m_fluid_Realm, m_phase);
  m_amr->interpGhostPwl(m_fluid_E,         m_fluid_Realm, m_phase);
  m_amr->interpToCentroids(m_fluid_E, m_fluid_Realm, m_phase);

  // Particle Realm
  m_particle_E.copy(E);
  m_amr->averageDown(m_particle_E,             m_particleRealm, m_phase);
  m_amr->interpGhostPwl(m_particle_E,         m_particleRealm, m_phase);
  m_amr->interpToCentroids(m_particle_E, m_particleRealm, m_phase);

  // Compute maximum E
  // const Real Emax = this->computeMaxElectricField(m_phase);
  // std::cout << Emax << std::endl;
}

void ito_plasma_stepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const {
  CH_TIME("ito_plasma_stepper::writeCheckpointData");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::writeCheckpointData" << endl;
  }

  for (ItoIterator<ItoSolver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<ItoSolver>& solver = solver_it();
    solver->writeCheckpointLevel(a_handle, a_lvl);
  }

  for (RtIterator<McPhoto> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<McPhoto>& solver = solver_it();
    solver->writeCheckpointLevel(a_handle, a_lvl);
  }

  m_fieldSolver->writeCheckpointLevel(a_handle, a_lvl);
  m_sigma->writeCheckpointLevel(a_handle, a_lvl);
}

void ito_plasma_stepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl){
  CH_TIME("ito_plasma_stepper::readCheckpointData");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::readCheckpointData" << endl;
  }

  for (ItoIterator<ItoSolver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>& solver = solver_it();
    solver->readCheckpointLevel(a_handle, a_lvl);
  }

  for (RtIterator<McPhoto> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<McPhoto>& solver = solver_it();
    solver->readCheckpointLevel(a_handle, a_lvl);
  }

  m_fieldSolver->readCheckpointLevel(a_handle, a_lvl);
  m_sigma->readCheckpointLevel(a_handle, a_lvl);
}

void ito_plasma_stepper::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const {
  CH_TIME("ito_plasma_stepper::writePlotData");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::writePlotData" << endl;
  }

  // Poisson solver copies over its output data
  a_plotVariableNames.append(m_fieldSolver->getPlotVariableNames());
  m_fieldSolver->writePlotData(a_output, a_icomp);

  // Surface charge solver writes
  a_plotVariableNames.append(m_sigma->getPlotVariableNames());
  m_sigma->writePlotData(a_output, a_icomp);

  // Ito solvers copy their output data
  for (ItoIterator<ItoSolver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>& solver = solver_it();
    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // RTE solvers copy their output data
  for (RtIterator<McPhoto> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<McPhoto>& solver = solver_it();
    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // Write the current to the output
  this->writeJ(a_output, a_icomp);
  a_plotVariableNames.push_back("x-J");
  a_plotVariableNames.push_back("y-J");
  if(SpaceDim == 3){
    a_plotVariableNames.push_back("z-J");
  }

  // Write the number of particles per patch
  this->write_num_particles_per_patch(a_output, a_icomp);
  a_plotVariableNames.push_back("particles_per_patch");
}

void ito_plasma_stepper::writeJ(EBAMRCellData& a_output, int& a_icomp) const{
  CH_TIME("ito_plasma_stepper::writeJ");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::writeJ" << endl;
  }

  const Interval src_interv(0, SpaceDim-1);
  const Interval dst_interv(a_icomp, a_icomp + SpaceDim -1);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    if(m_J.getRealm() == a_output.getRealm()){
      m_J[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else {
      m_J[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }
  a_icomp += SpaceDim;
}

void ito_plasma_stepper::write_num_particles_per_patch(EBAMRCellData& a_output, int& a_icomp) const {
  CH_TIME("ito_plasma_stepper::write_num_particles_per_patch");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::write_num_particles_per_patch" << endl;
  }

  const Interval src_interv(0, 0);
  const Interval dst_interv(a_icomp, a_icomp);

  DataOps::setValue(m_particle_scratch1, 0.0);
  
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const ParticleContainer<ItoParticle>& particles = solver_it()->getParticles(ItoSolver::WhichContainer::bulk);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	(*m_particle_scratch1[lvl])[dit()] += particles[lvl][dit].numItems();
      }
    }
  }

  // Copy to output holder
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    if(m_particle_scratch1.getRealm() == a_output.getRealm()){
      m_particle_scratch1[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else {
      m_particle_scratch1[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }
  
  a_icomp += 1;
}

void ito_plasma_stepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt){
  CH_TIME("ito_plasma_stepper::synchronizeSolverTimes");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::synchronizeSolverTimes" << endl;
  }

  m_timeStep = a_step;
  m_time = a_time;
  m_dt   = a_dt;

  m_ito->setTime(a_step,     a_time, a_dt);
  m_fieldSolver->setTime(a_step, a_time, a_dt);
  m_rte->setTime(a_step,     a_time, a_dt);
  m_sigma->setTime(a_step,   a_time, a_dt);
}

void ito_plasma_stepper::printStepReport(){
  CH_TIME("ito_plasma_stepper::printStepReport");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::printStepReport" << endl;
  }

  const Real Emax = this->computeMaxElectricField(m_phase);
  
  const size_t l_particles        = m_ito->getNumParticles(ItoSolver::WhichContainer::bulk, true);
  const size_t g_particles        = m_ito->getNumParticles(ItoSolver::WhichContainer::bulk, false);
  
  const size_t l_eb_particles     = m_ito->getNumParticles(ItoSolver::WhichContainer::eb, true);
  const size_t g_eb_particles     = m_ito->getNumParticles(ItoSolver::WhichContainer::eb, false);
  
  const size_t l_domain_particles = m_ito->getNumParticles(ItoSolver::WhichContainer::domain, true);
  const size_t g_domain_particles = m_ito->getNumParticles(ItoSolver::WhichContainer::domain, false);

  const size_t l_source_particles = m_ito->getNumParticles(ItoSolver::WhichContainer::source, true);
  const size_t g_source_particles = m_ito->getNumParticles(ItoSolver::WhichContainer::source, false);

  Real avg;
  Real sigma;
  
  int minRank;
  int maxRank;
  
  size_t minParticles;
  size_t maxParticles;
  
  // Compute some particle statistics
  this->get_particle_statistics(avg, sigma, minParticles, maxParticles, minRank, maxRank);

  // How was the time step restricted
  std::string str;
  switch(m_timeCode){
  case TimeCode::Physics:
    str = "dt restricted by 'physics'";
    break;
  case TimeCode::Advection:
    str = "dt restricted by 'advection'";
    break;
  case TimeCode::RelaxationTime:
    str = "dt restricted by 'relaxation time'";
    break;
  case TimeCode::Hardcap:
    str = "dt restricted by 'hardcap'";
    break;
  default:
    str = "dt restricted by 'unspecified'";
    break;
  }
  pout() << "                                   " + str << endl;
  pout() << "                                   Emax      = " << Emax << endl
	 << "                                   #part     = " << l_particles << " (" << g_particles << ")" << endl
	 << "                                   #eb part  = " << l_eb_particles << " (" << g_eb_particles << ")" << endl
	 << "                                   #dom part = " << l_domain_particles << " (" << g_domain_particles << ")" << endl
	 << "                                   #src part = " << l_source_particles << " (" << g_source_particles << ")" << endl
	 << "                                   #part min = " << minParticles << " (on rank = " << minRank << ")" << endl
	 << "                                   #part max = " << maxParticles << " (on rank = " << maxRank << ")" << endl
	 << "                                   #part avg = " << avg << endl
	 << "                                   #part dev = " << sigma << " (" << 100.*sigma/avg << "%)" << endl;
}

void ito_plasma_stepper::get_particle_statistics(Real& a_avg, Real& a_sigma, size_t& a_minPart, size_t& a_maxPart, int& a_minRank, int& a_maxRank){
  CH_TIME("ito_plasma_stepper::get_particle_statistics");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::get_particle_statistics" << endl;
  }

  const int srcProc = 0; 
  const int nProc   = numProc();

  const size_t numLocal = m_ito->getNumParticles(ItoSolver::WhichContainer::bulk, true);

  // Gather on source proc
  Vector<size_t> allCounts(nProc);
  gather(allCounts, numLocal, srcProc);


  // Compute and broadcast the average and the standard deviation
  if(procID() == srcProc){

    a_avg = 0.0;
    for (int i = 0; i < nProc; i++){
      a_avg += 1.0*allCounts[i];
    }
    a_avg *= 1./nProc;

    a_sigma = 0.0;
    for (int i = 0; i < nProc; i++){
      a_sigma += std::pow(1.0*allCounts[i]-a_avg,2);
    }
    a_sigma = sqrt(a_sigma/nProc);
  }
  
  broadcast(a_avg,   srcProc);
  broadcast(a_sigma, srcProc);

  
  // Get the minimum/maximum number of particles
  a_minRank = srcProc;
  a_maxRank = srcProc;
  
  a_minPart = std::numeric_limits<size_t>::max();
  a_maxPart = 0;

  if(procID() == srcProc){
    for (int i = 0; i < nProc; i++){
      if(allCounts[i] < a_minPart){
	a_minPart = allCounts[i];
	a_minRank = i;
      }

      if(allCounts[i] > a_maxPart){
	a_maxPart = allCounts[i];
	a_maxRank = i;
      }
    }
  }

  broadcast(a_minRank, srcProc);
  broadcast(a_maxRank, srcProc);
  
  broadcast(a_minPart, srcProc);
  broadcast(a_maxPart, srcProc);
}

void ito_plasma_stepper::print_timer_diagnostics(Real& a_timer, const std::string a_prefix){
  CH_TIME("ito_plasma_stepper::print_timer_diagnostics");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::print_timer_diagnostics" << endl;
  }

  const int srcProc = 0; 
  const int nProc   = numProc();

  const size_t numLocal = m_ito->getNumParticles(ItoSolver::WhichContainer::bulk, true);

  // Gather all timers on source proc
  Vector<Real> allTimers(nProc);
  gather(allTimers, a_timer, srcProc);

  Real averageTime;
  Real sigmaTime;

  Real minTime;
  Real maxTime;

  int minRank;
  int maxRank;

  // Compute and broadcast the average and the standard deviation
  if(procID() == srcProc){

    averageTime = 0.0;
    for (int i = 0; i < nProc; i++){
      averageTime += allTimers[i];
    }
    averageTime *= 1./nProc;

    sigmaTime = 0.0;
    for (int i = 0; i < nProc; i++){
      sigmaTime += std::pow(1.0*allTimers[i]-averageTime,2);
    }
    sigmaTime = sqrt(sigmaTime/nProc);
  }
  
  broadcast(averageTime, srcProc);
  broadcast(sigmaTime,   srcProc);
  
  // Get the minimum/maximum number of particles
  minRank = srcProc;
  maxRank = srcProc;
  
  minTime = std::numeric_limits<Real>::max();
  maxTime = 0;

  if(procID() == srcProc){
    for (int i = 0; i < nProc; i++){
      if(allTimers[i] < minTime){
	minTime = allTimers[i];
	minRank = i;
      }

      if(allTimers[i] > maxTime){
	maxTime = allTimers[i];
	maxRank = i;
      }
    }
  }

  broadcast(minRank, srcProc);
  broadcast(maxRank, srcProc);
  
  broadcast(minTime, srcProc);
  broadcast(maxTime, srcProc);

  // Fix formatting for the various fields
  std::stringstream ss_locTime;
  std::stringstream ss_minTime;
  std::stringstream ss_maxTime;
  std::stringstream ss_avgTime;
  std::stringstream ss_sigTime;

  ss_locTime << std::fixed << std::setprecision(2) << a_timer;
  ss_minTime << std::fixed << std::setprecision(2) << minTime;
  ss_maxTime << std::fixed << std::setprecision(2) << maxTime;
  ss_avgTime << std::fixed << std::setprecision(2) << averageTime;
  ss_sigTime << std::fixed << std::setprecision(2) << sigmaTime;


  pout() << std::left << std::setw(25) << a_prefix
	 << " | " << std::right << std::setw(8) << ss_locTime.str()
	 << " | " << std::right << std::setw(8) << ss_minTime.str()
	 << " | " << std::right << std::setw(8) << ss_maxTime.str()
	 << " | " << std::right << std::setw(8) << ss_avgTime.str()
	 << " | " << std::right << std::setw(8) << ss_sigTime.str()
	 << " | " << std::right << std::setw(8) << minRank
	 << " | " << std::right << std::setw(8) << maxRank
	 << " | " << endl;
}

void ito_plasma_stepper::print_timer_head(){
  pout() << "--------------------------------------------------------------------------------------------------------"
	 << endl
	 << std::left << std::setw(25) << "Kernel"
	 << " | " << std::right << std::setw(8) << "Loc."
	 << " | " << std::right << std::setw(8) << "Min."
	 << " | " << std::right << std::setw(8) << "Max."
	 << " | " << std::right << std::setw(8) << "Avg." 
	 << " | " << std::right << std::setw(8) << "Dev."
	 << " | " << std::right << std::setw(8)  << "Min rank"
	 << " | " << std::right << std::setw(8)  << "Max rank"
	 << " | " << endl
	 << "-------------------------------------------------------------------------------------------------------|"
	 << endl;
}

void ito_plasma_stepper::print_timer_tail(){
  pout() << "--------------------------------------------------------------------------------------------------------\n";
}

void ito_plasma_stepper::parse_filters(){
  CH_TIME("ito_plasma_stepper::computeDt");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::parse_filters" << endl;
  }

  ParmParse pp(m_name.c_str());

  // Build filters. Always uses a compensation step.
  for (int i = 0; i < 100; i++){

    const int ndigits = round(log10(1.0 + 1.0*i));
    char* cstr = new char[1+ndigits];
    sprintf(cstr, "%d", 1 + i);

    const std::string str = "filter_" + std::string(cstr);

    if(pp.contains(str.c_str())){
      Real alpha;
      int stride;
      int N;
      bool comp;
    
      pp.get(str.c_str(), alpha,    0);
      pp.get(str.c_str(), stride,   1);
      pp.get(str.c_str(), N,        2);
      pp.get(str.c_str(), comp,     3);

      m_filters.emplace_front(alpha, stride, N);
      if(comp){
	const Real alphaC = (N+1) - N*alpha;
	m_filters.emplace_front(alphaC, stride, 1);
      }
    }

    delete cstr;
  }
}

void ito_plasma_stepper::computeDt(Real& a_dt, TimeCode& a_timeCode){
  CH_TIME("ito_plasma_stepper::computeDt");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeDt" << endl;
  }
  
  a_dt = m_ito->computeDt();
  a_dt = a_dt*m_max_cells_hop;
  a_timeCode = TimeCode::Advection;
  
}

void ito_plasma_stepper::registerRealms(){
  CH_TIME("ito_plasma_stepper::registerRealms");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::registerRealms" << endl;
  }

  m_amr->registerRealm(m_fluid_Realm);
  m_amr->registerRealm(m_particleRealm);
}

void ito_plasma_stepper::registerOperators(){
  CH_TIME("ito_plasma_stepper::registerOperators");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::registerOperators" << endl;
  }

  m_ito->registerOperators();
  m_fieldSolver->registerOperators();
  m_rte->registerOperators();
  m_sigma->registerOperators();

  m_amr->registerMask("particle_halo", m_halo_buffer, m_particleRealm);
}
  
void ito_plasma_stepper::preRegrid(const int a_lmin, const int a_oldFinestLevel){
  CH_TIME("ito_plasma_stepper::preRegrid");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::preRegrid" << endl;
  }

  m_ito->preRegrid(a_lmin,     a_oldFinestLevel);
  m_fieldSolver->preRegrid(a_lmin, a_oldFinestLevel);
  m_rte->preRegrid(a_lmin,     a_oldFinestLevel);
  m_sigma->preRegrid(a_lmin,   a_oldFinestLevel);
}

void ito_plasma_stepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("ito_plasma_stepper::regrid");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::regrid" << endl;
  }

  // Allocate new memory
  this->allocateInternals();

  // Regrid solvers
  m_ito->regrid(a_lmin,     a_oldFinestLevel, a_newFinestLevel);
  m_fieldSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_rte->regrid(a_lmin,     a_oldFinestLevel, a_newFinestLevel);
  m_sigma->regrid(a_lmin,   a_oldFinestLevel, a_newFinestLevel);

  if(m_regrid_superparticles){
    m_ito->sortParticlesByCell( ItoSolver::WhichContainer::bulk);
    m_ito->makeSuperparticles(    ItoSolver::WhichContainer::bulk, m_ppc);
    m_ito->sortParticlesByPatch(ItoSolver::WhichContainer::bulk);
  }

  // Redeposit particles
  m_ito->depositParticles();

  // Recompute the electric field
  const bool converged = this->solvePoisson();
  if(!converged){
    MayDay::Abort("ito_plasma_stepper::regrid - Poisson solve did not converge after regrid!!!");
  }

  // Recompute new velocities and diffusion coefficients
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
}

void ito_plasma_stepper::postRegrid(){

}

int ito_plasma_stepper::getNumberOfPlotVariables() const {
  CH_TIME("ito_plasma_stepper::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::getNumberOfPlotVariables" << endl;
  }

  int ncomp = 0;
  
  for (ItoIterator<ItoSolver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>& solver = solver_it();
    ncomp += solver->getNumberOfPlotVariables();
  }
  
  for (RtIterator<McPhoto> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<McPhoto>& solver = solver_it();
    ncomp += solver->getNumberOfPlotVariables();
  }

  ncomp += m_fieldSolver->getNumberOfPlotVariables();
  ncomp += m_sigma->getNumberOfPlotVariables();
  ncomp += SpaceDim; // For plotting the current density
  ncomp += 1;        // For plotting the number of particles per cell

  return ncomp;
}

void ito_plasma_stepper::set_ito(RefCountedPtr<ItoLayout<ItoSolver> >& a_ito){
  CH_TIME("ito_plasma_stepper::set_ito");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_ito" << endl;
  }

  m_ito = a_ito;
}

void ito_plasma_stepper::setFieldSolver(RefCountedPtr<FieldSolver>& a_poisson){
  CH_TIME("ito_plasma_stepper::setFieldSolver");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setFieldSolver" << endl;
  }

  m_fieldSolver = a_poisson;
}

void ito_plasma_stepper::setRadiativeTransferSolvers(RefCountedPtr<RtLayout<McPhoto> >& a_rte){
  CH_TIME("ito_plasma_stepper::setRadiativeTransferSolvers");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setRadiativeTransferSolvers" << endl;
  }
  
  m_rte = a_rte;
}

void ito_plasma_stepper::setVoltage(const std::function<Real(const Real a_time)>& a_potential){
  CH_TIME("ito_plasma_stepper::setVoltage");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setVoltage" << endl;
  }

  m_potential = a_potential;
}

Real ito_plasma_stepper::computeMaxElectricField(const phase::which_phase a_phase) {
  CH_TIME("ito_plasma_stepper::computeMaxElectricField");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeMaxElectricField" << endl;
  }

  // Get a handle to the E-field
  EBAMRCellData Ephase;
  m_amr->allocatePointer(Ephase);
  m_amr->alias(Ephase, m_phase, m_fieldSolver->getElectricField());

  // Interpolate to centroids
  EBAMRCellData E;
  m_amr->allocate(E, m_fluid_Realm, m_phase, SpaceDim);
  DataOps::copy(E, Ephase);
  m_amr->interpToCentroids(E, m_fluid_Realm, m_phase);

  Real max, min;
  DataOps::getMaxMinNorm(max, min, E);

  return max;
}

Real ito_plasma_stepper::getTime() const{
  CH_TIME("ito_plasma_stepper::getTime");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::getTime" << endl;
  }

  return m_time;
}

void ito_plasma_stepper::computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential){
  CH_TIME("ito_plasma_stepper::computeElectricField(mfamrcell,mfamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeElectricField(mfamrcell, mfamrcell" << endl;
  }

  m_amr->computeGradient(a_E, a_potential, m_fluid_Realm);
  DataOps::scale(a_E, -1.0);

  m_amr->averageDown(a_E, m_fluid_Realm);
  m_amr->interpGhost(a_E, m_fluid_Realm);

}

void ito_plasma_stepper::computeElectricField(EBAMRCellData& a_E, const phase::which_phase a_phase){
  CH_TIME("ito_plasma_stepper::computeElectricField(ebamrcell, phase)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeElectricField(ebamrcell, phase" << endl;
  }

  this->computeElectricField(a_E, a_phase, m_fieldSolver->getPotential());
}

void ito_plasma_stepper::computeElectricField(EBAMRCellData& a_E, const phase::which_phase a_phase, const MFAMRCellData& a_potential){
  CH_TIME("ito_plasma_stepper::computeElectricField(ebamrcell, phase, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeElectricField(ebamrcell, phase mfamrcell" << endl;
  }
  
  EBAMRCellData pot_gas;
  m_amr->allocatePointer(pot_gas);
  m_amr->alias(pot_gas, a_phase, a_potential);

  m_amr->computeGradient(a_E, pot_gas, m_fluid_Realm, a_phase);
  DataOps::scale(a_E, -1.0);

  m_amr->averageDown(a_E, m_fluid_Realm, a_phase);
  m_amr->interpGhost(a_E, m_fluid_Realm, a_phase);
}

void ito_plasma_stepper::computeElectricField(EBAMRFluxData& a_E_face, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell){
  CH_TIME("ito_plasma_stepper::computeElectricField(ebamrflux, phase, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeElectricField(ebamrflux, phase mfamrcell" << endl;
  }

  CH_assert(a_E_face[0]->nComp() == SpaceDim);
  CH_assert(a_E_cell[0]->nComp() == SpaceDim);

  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++){

    const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluid_Realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_fluid_Realm, a_phase)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];

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

void ito_plasma_stepper::computeElectricField(EBAMRIVData& a_E_eb,  const phase::which_phase a_phase, const EBAMRCellData& a_E_cell){
  CH_TIME("ito_plasma_stepper::computeElectricField(ebamriv, phase, ebamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeElectricField(ebamriv, phase ebamrcell)" << endl;
  }

  CH_assert(a_E_eb[0]->nComp()   == SpaceDim);
  CH_assert(a_E_cell[0]->nComp() == SpaceDim);

  const IrregAmrStencil<EbCentroidInterpolationStencil>& interp_stencil = m_amr->getEbCentroidInterpolationStencilStencils(m_fluid_Realm, a_phase);
  interp_stencil.apply(a_E_eb, a_E_cell);
}

void ito_plasma_stepper::computeSpaceChargeDensity(){
  CH_TIME("ito_plasma_stepper::computeSpaceChargeDensity()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeSpaceChargeDensity()" << endl;
  }
  
  this->computeSpaceChargeDensity(m_fieldSolver->getRho(), m_ito->getDensities());
}

void ito_plasma_stepper::computeSpaceChargeDensity(MFAMRCellData& a_rho, const Vector<EBAMRCellData*>&  a_densities){
  CH_TIME("ito_plasma_stepper::computeSpaceChargeDensity(rho, densities)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeSpaceChargeDensity(rho, densities)" << endl;
  }

  // TLDR: a_densities is from the ito solvers so it is defined over the particle Realm. But a_rho is defined over
  //       the fluid Realm so we need scratch storage we can copy into. We use m_fluid_scratch1 for that. 

  // Reset
  DataOps::setValue(a_rho, 0.0);

  // Make alias
  EBAMRCellData rhoPhase;
  m_amr->allocatePointer(rhoPhase);
  m_amr->alias(rhoPhase, m_phase, a_rho);

  // Increment each solver
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<ItoSolver>& solver   = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();
    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    if(species->getChargeNumber() != 0){
      m_fluid_scratch1.copy(*a_densities[idx]);
      DataOps::incr(rhoPhase, m_fluid_scratch1, q);
    }
  }

  DataOps::scale(a_rho, Units::Qe);

  m_amr->averageDown(a_rho, m_fluid_Realm);
  m_amr->interpGhost(a_rho, m_fluid_Realm);


  // Add potential filters.
  if(false){//m_filter_rho){
    for (const auto& f : m_filters){
      const Real alpha  = std::get<0>(f);
      const int stride  = std::get<1>(f);
      const int num_app = std::get<2>(f);

      for (int iapp = 0; iapp < num_app; iapp++){
	DataOps::setValue(m_fluid_scratch1, 0.0);
	m_fluid_scratch1.copy(rhoPhase);
	DataOps::setCoveredValue(m_fluid_scratch1, 0.0, 0);
	DataOps::filterSmooth(rhoPhase, m_fluid_scratch1, stride, alpha);

	m_amr->averageDown(rhoPhase, m_fluid_Realm, m_phase);
	m_amr->interpGhost(rhoPhase, m_fluid_Realm, m_phase);
      }
    }
  }

  // Interpolate to centroids
  m_amr->interpToCentroids(rhoPhase, m_fluid_Realm, m_phase);
}

void ito_plasma_stepper::compute_conductivity(EBAMRCellData& a_conductivity){
  CH_TIME("ito_plasma_stepper::compute_conductivity(conductivity)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_conductivity(conductivity)" << endl;
  }

  this->compute_conductivity(a_conductivity, m_ito->getParticles(ItoSolver::WhichContainer::bulk));
  
}

void ito_plasma_stepper::compute_conductivity(EBAMRCellData& a_conductivity, const Vector<ParticleContainer<ItoParticle>* >& a_particles){
  CH_TIME("ito_plasma_stepper::compute_conductivity(conductivity)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_conductivity(conductivity)" << endl;
  }
  
  DataOps::setValue(a_conductivity, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>&   solver = solver_it();
    RefCountedPtr<ItoSpecies>& species = solver->getSpecies();
    
    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    if(Abs(q) > 0 && solver->isMobile()){
      DataOps::setValue(m_particle_scratch1, 0.0);

      solver->depositConductivity(m_particle_scratch1, *a_particles[idx]);

      // Copy to fluid Realm and add to fluid stuff
      m_fluid_scratch1.copy(m_particle_scratch1);
      DataOps::incr(a_conductivity, m_fluid_scratch1, Abs(q));
    }
  }

  DataOps::scale(a_conductivity, Units::Qe);

  m_amr->averageDown(a_conductivity, m_fluid_Realm, m_phase);
  m_amr->interpGhostPwl(a_conductivity, m_fluid_Realm, m_phase);

  // See if this helps....
  m_amr->interpToCentroids(a_conductivity, m_fluid_Realm, m_phase);
}

void ito_plasma_stepper::computeJ(EBAMRCellData& a_J, const Real a_dt){
  CH_TIME("ito_plasma_stepper::computeJ(J)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeJ(J)" << endl;
  }

  // TLDR: a_J is defined over the fluid Realm but the computation takes place on the particle Realm.
  //       If the Realms are different we compute on a scratch storage instead


  this->compute_conductivity(m_fluid_scratch1);
  DataOps::copy(a_J, m_fluid_E);

  DataOps::multiplyScalar(a_J, m_fluid_scratch1);
}

Real ito_plasma_stepper::computeRelaxationTime(){
  CH_TIME("ito_plasma_stepper::computeRelaxationTime()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeRelaxationTime()" << endl;
  }

  Real dt = 1.E99;
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const Real thisDt = this->computeRelaxationTime(lvl);

    dt = Min(dt, thisDt);
  }

  return dt;
}

Real ito_plasma_stepper::computeRelaxationTime(const int a_level){
  CH_TIME("ito_plasma_stepper::computeRelaxationTime(level)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeRelaxationTime(level)" << endl;
  }

  Real dt = 1.E99;

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluid_Realm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Real thisDt = this->computeRelaxationTime(a_level, dit());

    dt = Min(dt, thisDt);
  }

#ifdef CH_MPI
  Real tmp = dt;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("CdrSolver::compute_cfl_dt() - communication error on norm");
  }
  dt = tmp;
#endif

  return dt;
}

Real ito_plasma_stepper::computeRelaxationTime(const int a_level, const DataIndex a_dit){
  CH_TIME("ito_plasma_stepper::computeRelaxationTime(level, dit)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeRelaxationTime(level, dit)" << endl;
  }

  const int comp    = 0;
  const Real SAFETY = 1.E-10;

  const Box box = m_amr->getGrids(m_fluid_Realm)[a_level].get(a_dit);
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_fluid_Realm, m_phase)[a_level][a_dit];

  // Get a handle to the E-field
  EBAMRCellData amrE;
  m_amr->allocatePointer(amrE);
  m_amr->alias(amrE, m_phase, m_fieldSolver->getElectricField());
  
  const EBCellFAB& E = (*amrE[a_level])[a_dit];
  const EBCellFAB& J = (*m_J[a_level])[a_dit];

  EBCellFAB dt(ebisbox, box, 1);
  EBCellFAB e_magnitude(ebisbox, box, 1);
  EBCellFAB j_magnitude(ebisbox, box, 1);

  e_magnitude.setVal(0.0);
  j_magnitude.setVal(0.0);

  DataOps::vectorLength(e_magnitude, E, box);
  DataOps::vectorLength(j_magnitude, J, box);
  j_magnitude += SAFETY;

  dt.setVal(Units::eps0);
  dt *= e_magnitude;
  dt /= j_magnitude;
  
  return dt.min(comp);
}

bool ito_plasma_stepper::solvePoisson(){
  CH_TIME("ito_plasma_stepper::solvePoisson()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::solvePoisson()" << endl;
  }

  // Computes rho
  this->computeSpaceChargeDensity();

  // Solves the Poisson equation
  const bool converged = m_fieldSolver->solve(m_fieldSolver->getPotential(),
					  m_fieldSolver->getRho(),
					  m_sigma->getPhi(),
					  false);

  // Computes cell-centered E onto storage in the field solver
  m_fieldSolver->computeElectricField();

  // Code below here interpolates E to centroids on both Realms
  EBAMRCellData E;
  m_amr->allocatePointer(E);
  m_amr->alias(E, m_phase, m_fieldSolver->getElectricField());
  
  // Fluid Realm
  m_fluid_E.copy(E);
  m_amr->averageDown(m_fluid_E, m_fluid_Realm, m_phase);
  m_amr->interpGhostPwl(m_fluid_E, m_fluid_Realm, m_phase);
  m_amr->interpToCentroids(m_fluid_E, m_fluid_Realm, m_phase);

  // Particle Realm
  m_particle_E.copy(E);
  m_amr->averageDown(m_particle_E, m_particleRealm, m_phase);
  m_amr->interpGhostPwl(m_particle_E, m_particleRealm, m_phase);
  m_amr->interpToCentroids(m_particle_E, m_particleRealm, m_phase);
    
  return converged;
}

bool ito_plasma_stepper::solvePoisson(MFAMRCellData&                a_potential,
				       MFAMRCellData&                a_rho,
				       const Vector<EBAMRCellData*>& a_densities,
				       const EBAMRIVData&            a_sigma){
  CH_TIME("ito_plasma_stepper::solvePoisson(phi, rho, densities, sigma)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::solvePoisson(phi, rho, densities, sigma)" << endl;
  }

  this->computeSpaceChargeDensity(a_rho, a_densities);

  const bool converged = m_fieldSolver->solve(a_potential,
					  a_rho,
					  a_sigma,
					  false);
  m_fieldSolver->computeElectricField();
  
  return converged;
}

void ito_plasma_stepper::intersectParticles(const which_particles a_which_particles, const EbRepresentation a_representation, const bool a_delete){
  CH_TIME("ito_plasma_stepper::intersectParticles(which_particles, EbRepresentation)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::intersectParticles(which_particles, EbRepresentation)" << endl;
  }

  this->intersectParticles(a_which_particles,
			    ItoSolver::WhichContainer::bulk,
			    ItoSolver::WhichContainer::eb,
			    ItoSolver::WhichContainer::domain,
			    a_representation,
			    a_delete);
}

void ito_plasma_stepper::intersectParticles(const which_particles             a_which_particles,
					     const ItoSolver::WhichContainer a_particles,
					     const ItoSolver::WhichContainer a_eb_particles,
					     const ItoSolver::WhichContainer a_domain_particles,
					     const EbRepresentation           a_representation,
					     const bool                        a_delete) {
  CH_TIME("ito_plasma_stepper::intersectParticles(which_particles, string, string, string, EbRepresentation, bool)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::intersectParticles(which_particles, string, string, string, EbRepresentation, bool)" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>&   solver = solver_it();
    RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch(a_which_particles) {
    case which_particles::all:
      solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::all_mobile:
      if(mobile) solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::all_diffusive:
      if(diffusive) solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::charged_mobile:
      if(charged && mobile) solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::charged_diffusive:
      if(charged && diffusive) solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::all_mobile_or_diffusive:
      if(mobile || diffusive) solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::charged_and_mobile_or_diffusive:
      if(charged && (mobile || diffusive)) solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::stationary:
      if(!mobile && !diffusive) solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    default:
      MayDay::Abort("ito_plasma_stepper::intersectParticles_particles(which_particles, string, string, string, EbRepresentation, bool) - logic bust");
    }
  }  
}

void ito_plasma_stepper::removeCoveredParticles(const which_particles a_which_particles, const EbRepresentation a_representation, const Real a_tolerance){
  CH_TIME("ito_plasma_stepper::removeCoveredParticles(which_particles, representation, tolerance)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::removeCoveredParticles(which_particles, representation, tolerance)" << endl;
  }

  this->removeCoveredParticles(a_which_particles, ItoSolver::WhichContainer::bulk, a_representation, a_tolerance);
}

void ito_plasma_stepper::removeCoveredParticles(const which_particles             a_which,
						  const ItoSolver::WhichContainer a_container,
						  const EbRepresentation           a_representation,
						  const Real                        a_tolerance){
  CH_TIME("ito_plasma_stepper::removeCoveredParticles(which_particles, container, EbRepresentation, tolerance)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::removeCoveredParticles(which_particles, container, EbRepresentation, tolerance)" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>&   solver = solver_it();
    RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch(a_which) {
    case which_particles::all:
      solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::all_mobile:
      if(mobile) solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::all_diffusive:
      if(diffusive) solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::charged_mobile:
      if(charged && mobile) solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::charged_diffusive:
      if(charged && diffusive) solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::all_mobile_or_diffusive:
      if(mobile || diffusive) solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::charged_and_mobile_or_diffusive:
      if(charged && (mobile || diffusive)) solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::stationary:
      if(!mobile && !diffusive) solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    default:
      MayDay::Abort("ito_plasma_stepper::removeCoveredParticles_particles(which particles) - logic bust");
    }
  }  
}

void ito_plasma_stepper::transferCoveredParticles(const which_particles a_which, const EbRepresentation a_representation, const Real a_tolerance){
  CH_TIME("ito_plasma_stepper::transferCoveredParticles_particles(which_particles, EbRepresentation, tolerance)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::transferCoveredParticles_particles(which_particles, EbRepresentation, tolerance)" << endl;
  }

  this->transferCoveredParticles(a_which, ItoSolver::WhichContainer::bulk, ItoSolver::WhichContainer::covered, a_representation, a_tolerance);
}

void ito_plasma_stepper::transferCoveredParticles(const which_particles             a_which,
						    const ItoSolver::WhichContainer a_containerFrom,
						    const ItoSolver::WhichContainer a_containerTo,
						    const EbRepresentation           a_representation,
						    const Real                        a_tolerance){
  CH_TIME("ito_plasma_stepper::transferCoveredParticles_particles(which_particles, container, container, EbRepresentation, tolerance)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::transferCoveredParticles_particles(which_particles, container, container, EbRepresentation, tolerance)" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>&   solver = solver_it();
    RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch(a_which) {
    case which_particles::all:
      solver->transferCoveredParticles(a_containerFrom, a_containerTo,a_representation, a_tolerance);
      break;
    case which_particles::all_mobile:
      if(mobile) solver->transferCoveredParticles(a_containerFrom, a_containerTo,a_representation, a_tolerance);
      break;
    case which_particles::all_diffusive:
      if(diffusive) solver->transferCoveredParticles(a_containerFrom, a_containerTo,a_representation, a_tolerance);
      break;
    case which_particles::charged_mobile:
      if(charged && mobile) solver->transferCoveredParticles(a_containerFrom, a_containerTo,a_representation, a_tolerance);
      break;
    case which_particles::charged_diffusive:
      if(charged && diffusive) solver->transferCoveredParticles(a_containerFrom, a_containerTo,a_representation, a_tolerance);
      break;
    case which_particles::all_mobile_or_diffusive:
      if(mobile || diffusive) solver->transferCoveredParticles(a_containerFrom, a_containerTo,a_representation, a_tolerance);
      break;
    case which_particles::charged_and_mobile_or_diffusive:
      if(charged && (mobile || diffusive)) solver->transferCoveredParticles(a_containerFrom, a_containerTo,a_representation, a_tolerance);
      break;
    case which_particles::stationary:
      if(!mobile && !diffusive) solver->transferCoveredParticles(a_containerFrom, a_containerTo,a_representation, a_tolerance);
      break;
    default:
      MayDay::Abort("ito_plasma_stepper::transferCoveredParticles_particles(...) - logic bust");
    }    
  }
}

void ito_plasma_stepper::remap_particles(const which_particles a_which_particles){
  CH_TIME("ito_plasma_stepper::remap_particles(which_particles)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::remap_particles(which_particles)" << endl;
  }

  this->remap_particles(a_which_particles, ItoSolver::WhichContainer::bulk);
}

void ito_plasma_stepper::remap_particles(const which_particles a_which_particles, const ItoSolver::WhichContainer a_container){
  CH_TIME("ito_plasma_stepper::remap_particles(which_particles)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::remap_particles(which_particles)" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>&   solver = solver_it();
    RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch(a_which_particles) {
    case which_particles::all:
      solver->remap(a_container);
      break;
    case which_particles::all_mobile:
      if(mobile) solver->remap(a_container);
      break;
    case which_particles::all_diffusive:
      if(diffusive) solver->remap(a_container);
      break;
    case which_particles::charged_mobile:
      if(charged && mobile) solver->remap(a_container);
      break;
    case which_particles::charged_diffusive:
      if(charged && diffusive) solver->remap(a_container);
      break;
    case which_particles::all_mobile_or_diffusive:
      if(mobile || diffusive) solver->remap(a_container);
      break;
    case which_particles::charged_and_mobile_or_diffusive:
      if(charged && (mobile || diffusive)) solver->remap(a_container);
      break;
    case which_particles::stationary:
      if(!mobile && !diffusive) solver->remap(a_container);
      break;
    default:
      MayDay::Abort("ito_plasma_stepper::remap_particles(which particles) - logic bust");
    }
  }
}

void ito_plasma_stepper::depositParticles(const which_particles a_which_particles){
  CH_TIME("ito_plasma_stepper::depositParticles(which_particles)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::depositParticles(which_particles)" << endl;
  }

  this->depositParticles(a_which_particles, ItoSolver::WhichContainer::bulk);

}

void ito_plasma_stepper::depositParticles(const which_particles a_which_particles, const ItoSolver::WhichContainer a_container){
  CH_TIME("ito_plasma_stepper::depositParticles(which_particles)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::depositParticles(which_particles)" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>&   solver = solver_it();
    RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch(a_which_particles) {
    case which_particles::all:
      solver->depositParticles(a_container);
      break;
    case which_particles::all_mobile:
      if(mobile) solver->depositParticles(a_container);
      break;
    case which_particles::all_diffusive:
      if(diffusive) solver->depositParticles(a_container);
      break;
    case which_particles::charged_mobile:
      if(charged && mobile) solver->depositParticles(a_container);
      break;
    case which_particles::charged_diffusive:
      if(charged && diffusive) solver->depositParticles(a_container);
      break;
    case which_particles::all_mobile_or_diffusive:
      if(mobile || diffusive) solver->depositParticles(a_container);
      break;
    case which_particles::charged_and_mobile_or_diffusive:
      if(charged && (mobile || diffusive)) solver->depositParticles(a_container);
      break;
    case which_particles::stationary:
      if(!mobile && !diffusive) solver->depositParticles(a_container);
      break;
    default:
      MayDay::Abort("ito_plasma_stepper::depositParticles(which_particles) - logic bust");
    }
  }
}

void ito_plasma_stepper::set_ito_velocity_funcs(){
  CH_TIME("ito_plasma_stepper::set_ito_velocity_funcs");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_ito_velocity_funcs" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>& solver   = solver_it();
    RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    if(solver->isMobile()){
      EBAMRCellData& velo_func = solver->getVelocityFunction();
      velo_func.copy(m_particle_E);
      
      const int q = species->getChargeNumber();
      const int s = (q > 0) - (q < 0);
      
      DataOps::scale(velo_func, s);
    }
  }
}

void ito_plasma_stepper::compute_ito_velocities(){
  CH_TIME("ito_plasma_stepper::compute_ito_velocities()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_velocities()" << endl;
  }

  const ItoPlasmaPhysics::coupling which_coupling = m_physics->getCoupling();

  // Set velocity functions
  this->set_ito_velocity_funcs();

  // Compute mobilities based on appropriate coupling
  switch(which_coupling){
  case ItoPlasmaPhysics::coupling::LFA:
    this->computeItoMobilitiesLFA();
    break;
  case ItoPlasmaPhysics::coupling::LEA:
    this->compute_ito_mobilities_lea();
    break;
  }

  // Interpolate velocity function to particle position
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->interpolateVelocities(); // Interpolates v = +/- mu*E
  }
}

void ito_plasma_stepper::compute_ito_diffusion(){
  CH_TIME("ito_plasma_stepper::compute_ito_diffusion()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_diffusion()" << endl;
  }

  const ItoPlasmaPhysics::coupling which_coupling = m_physics->getCoupling();

  // Compute mobilities based on appropriate coupling
  switch(which_coupling){
  case ItoPlasmaPhysics::coupling::LFA:
    this->computeItoDiffusionLFA();
    break;
  case ItoPlasmaPhysics::coupling::LEA:
    this->compute_ito_diffusion_lea();
    break;
  }
}

void ito_plasma_stepper::computeItoMobilitiesLFA(){
  CH_TIME("ito_plasma_stepper::computeItoMobilitiesLFA()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeItoMobilitiesLFA()" << endl;
  }

  Vector<EBAMRCellData*> meshMobilities = m_ito->getMobilityFunctions();
  this->computeItoMobilitiesLFA(meshMobilities, m_fluid_E, m_time);
}

void ito_plasma_stepper::computeItoMobilitiesLFA(Vector<EBAMRCellData*>& a_meshMobilities, const EBAMRCellData& a_E, const Real a_time){
  CH_TIME("ito_plasma_stepper::computeItoMobilitiesLFA(mobilities, E, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeItoMobilitiesLFA(mobilities, E, time)" << endl;
  }


  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){

    // Computation is done on the fluid Realm. 
    Vector<LevelData<EBCellFAB>* > meshMobilities;
    for (int i = 0; i < a_meshMobilities.size(); i++){
      meshMobilities.push_back(&(*(m_fscratch1[i])[lvl]));
    }
    
    this->computeItoMobilitiesLFA(meshMobilities, *a_E[lvl], lvl, a_time);
  }

  // Average down and interpolate ghost cells. Then interpolate mobilities to particle positions. 
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    RefCountedPtr<ItoSolver>& solver = solver_it();

    if(solver->isMobile()){

      
#if 0 // In principle, we should be able to average down and interpolate on the fluid Realm and then copy directly to the particle Realm.
      // But we need to make sure that EBAMRData::copy also gets ghost cells 
      m_amr->averageDown(m_fscratch1[idx], m_fluid_Realm, m_phase);
      m_amr->interpGhost(m_fscratch1[idx], m_fluid_Realm, m_phase);

      a_meshMobilities[idx]->copy(m_fscratch1[idx]);
#else
      // Copy to particle Realm, build ghost cells and the interpolate the mobilities to particle positions. 
      a_meshMobilities[idx]->copy(m_fscratch1[idx]);
      
      m_amr->averageDown(*a_meshMobilities[idx], m_particleRealm, m_phase);
      m_amr->interpGhost(*a_meshMobilities[idx], m_particleRealm, m_phase);
#endif

      solver->interpolateMobilities();
    }
  }
}

void ito_plasma_stepper::computeItoMobilitiesLFA(Vector<LevelData<EBCellFAB>* >& a_meshMobilities,
						    const LevelData<EBCellFAB>&     a_E,
						    const int                       a_level,
						    const Real                      a_time){
  CH_TIME("ito_plasma_stepper::computeItoMobilitiesLFA(mobilities, E, level, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeItoMobilitiesLFA(mobilities, E, level, time)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluid_Realm)[a_level];
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const EBCellFAB& E = a_E[dit()];
    const Box bx       = dbl.get(dit());

    Vector<EBCellFAB*> meshMobilities;
    for (int i = 0; i < a_meshMobilities.size(); i++){
      meshMobilities.push_back(&((*a_meshMobilities[i])[dit()]));
    }

    this->computeItoMobilitiesLFA(meshMobilities, E, a_level, dit(), bx, a_time);
  }
}

void ito_plasma_stepper::computeItoMobilitiesLFA(Vector<EBCellFAB*>& a_meshMobilities,
						    const EBCellFAB&    a_E,
						    const int           a_level,
						    const DataIndex     a_dit,
						    const Box           a_box,
						    const Real          a_time) {
  CH_TIME("ito_plasma_stepper::computeItoMobilitiesLFA(meshMobilities, E, level, dit, box, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeItoMobilitiesLFA(meshMobilities, E, level, dit, box, time)" << endl;
  }

  const int comp         = 0;
  const Real dx          = m_amr->getDx()[a_level];
  const RealVect prob_lo = m_amr->getProbLo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();

  // Do regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv   = bit();
    const RealVect pos = m_amr->getProbLo() + dx*(RealVect(iv) + 0.5*RealVect::Unit);
    const RealVect e   = RealVect(D_DECL(E(iv,0), E(iv, 1), E(iv, 2)));

    // Call ito_physics and compute diffusion for each particle species
    const Vector<Real> mobilities = m_physics->computeItoMobilitiesLFA(a_time, pos, e);
    
    // Put mobilities in data holder
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      (*a_meshMobilities[idx]).getSingleValuedFAB()(iv, comp) = mobilities[idx];
    }
  }

  // Do irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_fluid_Realm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect e    = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, prob_lo);
    
    // Compute diffusion
    const Vector<Real> mobilities = m_physics->computeItoMobilitiesLFA(a_time, pos, e);

    // Put diffusion in the appropriate place.
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      (*a_meshMobilities[idx])(vof, comp) = mobilities[idx];
    }
  }
  

  // Covered is bogus.
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    a_meshMobilities[idx]->setCoveredCellVal(0.0, comp);
  }
}

void ito_plasma_stepper::compute_ito_mobilities_lea(){
  CH_TIME("ito_plasma_stepper");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_mobilities_lea()" << endl;
  }

  // This is really simple because the solvers do this directly...
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->updateMobilities();
  }
}

void ito_plasma_stepper::computeItoDiffusionLFA(){
  CH_TIME("ito_plasma_stepper::computeItoDiffusionLFA()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeItoDiffusionLFA()" << endl;
  }

  Vector<EBAMRCellData*> diffco_funcs = m_ito->getDiffusionFunctions();
  Vector<EBAMRCellData*> densities    = m_ito->getDensities();

  this->computeItoDiffusionLFA(diffco_funcs, densities, m_particle_E, m_time);
}

void ito_plasma_stepper::computeItoDiffusionLFA(Vector<EBAMRCellData*>&       a_diffusionCoefficient_funcs,
						   const Vector<EBAMRCellData*>& a_densities,
						   const EBAMRCellData&          a_E,
						   const Real                    a_time){
  CH_TIME("ito_plasma_stepper::computeItoDiffusionLFA(velo, E, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeItoDiffusionLFA(velo, E, time)" << endl;
  }

  // TLDR: In this routine we make m_fscratch1 hold the diffusion coefficients on the fluid Realm and m_fscratch2 hold the particle densities on the fluid Realm.
  //       This requires a couple of copies. 

  // 1. Copy particle Realm densities to fluid Realm scratch data. 
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();

    m_fscratch2[idx].copy(*a_densities[idx]);
  }

  // 2. Compute on each level. On the fluid Realm. 
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){

    const int num_ItoSpecies = m_physics->getNumItoSpecies();

    Vector<LevelData<EBCellFAB>* > diffco_funcs(num_ItoSpecies);
    Vector<LevelData<EBCellFAB>* > densities(num_ItoSpecies);
    
    for (int idx = 0; idx < a_diffusionCoefficient_funcs.size(); idx++){
      diffco_funcs[idx] = &(*(m_fscratch1[idx])[lvl]);
      densities[idx]    = &(*(m_fscratch2[idx])[lvl]);
    }

    this->computeItoDiffusionLFA(diffco_funcs, densities, *m_fluid_E[lvl], lvl, a_time);
  }

  // Average down, interpolate ghost cells, and then interpolate to particle positions
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    RefCountedPtr<ItoSolver>& solver = solver_it();

    if(solver->isDiffusive()){

#if 0 // In principle, we should be able to average down and interpolate ghost cells on the fluid Realm, and copy the entire result over to the particle Realm.
      m_amr->averageDown(m_fscratch1[idx], m_fluid_Realm, m_phase);
      m_amr->interpGhost(m_fscratch2[idx], m_fluid_Realm, m_phase);
      a_diffusionCoefficient_funcs[idx]->copy(m_fluid_scratch1[idx]);
#else // Instead, we copy to the particle Realm and average down there, then interpolate. 
      a_diffusionCoefficient_funcs[idx]->copy(m_fscratch1[idx]);
      
      m_amr->averageDown(*a_diffusionCoefficient_funcs[idx], m_particleRealm, m_phase);
      m_amr->interpGhost(*a_diffusionCoefficient_funcs[idx], m_particleRealm, m_phase);
#endif

      solver->interpolateDiffusion();
    }
  }
}

void ito_plasma_stepper::computeItoDiffusionLFA(Vector<LevelData<EBCellFAB>* >&       a_diffusionCoefficient,
						   const Vector<LevelData<EBCellFAB>* >& a_densities,
						   const LevelData<EBCellFAB>&           a_E,
						   const int                             a_level,
						   const Real                            a_time){
  CH_TIME("ito_plasma_stepper::computeItoDiffusionLFA(velo, E, level, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeItoDiffusionLFA(velo, E, level, time)" << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluid_Realm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());

    Vector<EBCellFAB*> diffusion(num_ItoSpecies);
    Vector<EBCellFAB*> densities(num_ItoSpecies);;
    
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();

      if(solver_it()->isDiffusive()){
	diffusion[idx] = &(*a_diffusionCoefficient[idx])[dit()];
      }
      densities[idx] = &(*a_densities[idx])[dit()];
    }

    this->computeItoDiffusionLFA(diffusion, densities, a_E[dit()], a_level, dit(), box, a_time);
  }
}

void ito_plasma_stepper::computeItoDiffusionLFA(Vector<EBCellFAB*>&       a_diffusionCoefficient,
						   const Vector<EBCellFAB*>& a_densities,
						   const EBCellFAB&          a_E,
						   const int                 a_level,
						   const DataIndex           a_dit,
						   const Box                 a_box,
						   const Real                a_time){
  CH_TIME("ito_plasma_stepper::computeItoDiffusionLFA(velo, E, level, dit, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeItoDiffusionLFA(velo, E, level, dit, time)" << endl;
  }

  const int comp         = 0;
  const Real dx          = m_amr->getDx()[a_level];
  const RealVect prob_lo = m_amr->getProbLo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();

  // Do regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv   = bit();
    const RealVect pos = m_amr->getProbLo() + dx*(RealVect(iv) + 0.5*RealVect::Unit);
    const RealVect e   = RealVect(D_DECL(E(iv,0), E(iv, 1), E(iv, 2)));

    // Make grid densities
    Vector<Real> densities;
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      densities.push_back((*a_densities[idx]).getSingleValuedFAB()(iv, comp));
    }

    // Call ito_physics and compute diffusion for each particle species
    const Vector<Real> diffusion = m_physics->computeItoDiffusionLFA(a_time, pos, e, densities);
    
    // Put diffusion where they belong
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<ItoSolver>& solver = solver_it();
      if(solver->isDiffusive()){
	const int idx = solver_it.index();
	(*a_diffusionCoefficient[idx]).getSingleValuedFAB()(iv, comp) = diffusion[idx];
      }
    }
  }

  // Do irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_fluid_Realm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect e    = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, prob_lo);
    
    // Get densities
    Vector<Real> densities;
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      densities.push_back((*a_densities[idx])(vof, comp));
    }
    
    // Compute diffusion
    const Vector<Real> diffusion = m_physics->computeItoDiffusionLFA(a_time, pos, e, densities);

    // Put diffusion in the appropriate place.
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      if(solver_it()->isDiffusive()){
	const int idx = solver_it.index();
	(*a_diffusionCoefficient[idx])(vof, comp) = diffusion[idx];
      }
    }
  }

  // Covered is bogus.
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    if(solver_it()->isDiffusive()){
      const int idx = solver_it.index();
      a_diffusionCoefficient[idx]->setCoveredCellVal(0.0, comp);
    }
  }
}

void ito_plasma_stepper::compute_ito_diffusion_lea(){
  CH_TIME("ito_plasma_stepper::compute_ito_diffusion_lea()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_diffusion_lea()" << endl;
  }
  
  // This is really simple because the solvers do this directly... No monkeying with interpolations or anything. 
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->updateDiffusion();
  }
}

void ito_plasma_stepper::compute_reactive_particles_per_cell(EBAMRCellData& a_ppc){
  CH_TIME("ito_plasma_stepper::compute_reactive_particles_per_cell(ppc)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_reactive_particles_per_cell(ppc)" << endl;
  }

  DataOps::setValue(a_ppc, 0.0);
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    this->compute_reactive_particles_per_cell(*a_ppc[lvl], lvl);
  }
}

void ito_plasma_stepper::compute_reactive_particles_per_cell(LevelData<EBCellFAB>& a_ppc, const int a_level){
  CH_TIME("ito_plasma_stepper::compute_reactive_particles_per_cell(ppc, lvl)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_reactive_particles_per_cell(ppc, lvl)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[a_level];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      
    const Box box = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];
    
    this->compute_reactive_particles_per_cell(a_ppc[dit()], a_level, dit(), box, ebisbox);
  }
}

void ito_plasma_stepper::compute_reactive_particles_per_cell(EBCellFAB& a_ppc, const int a_level, const DataIndex a_dit, const Box a_box, const EBISBox& a_ebisbox){
  CH_TIME("ito_plasma_stepper::compute_reactive_particles_per_cell(ppc, lvl, dit, box)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_reactive_particles_per_cell(ppc, lvl, dit, box)" << endl;
  }

  BaseFab<Real>& numFab = a_ppc.getSingleValuedFAB();

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>& solver = solver_it();
    const int idx                     = solver_it.index();

    const ParticleContainer<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::bulk);
    const BinFab<ItoParticle>& cellParticles         = particles.getCellParticles(a_level, a_dit);


    // Regular cells. 
    for (BoxIterator bit(a_box); bit.ok(); ++bit){
      const IntVect iv = bit();

      Real num = 0.0;
      
      if(a_ebisbox.isRegular(iv)){

	
	const List<ItoParticle>& listParticles = cellParticles(iv, 0);
	for (ListIterator<ItoParticle> lit(listParticles); lit.ok(); ++lit){
	  num += lit().mass();
	}
      }
      
      numFab(iv, idx) = num;
    }

    // Irregular cells.
    VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_phase)[a_level])[a_dit];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof       = vofit();
      const IntVect iv          = vof.gridIndex();
      const RealVect normal     = a_ebisbox.normal(vof);
      const RealVect ebCentroid = a_ebisbox.bndryCentroid(vof);

      Real num = 0.0;

      const List<ItoParticle>& listParticles = cellParticles(iv, 0);
      
      for (ListIterator<ItoParticle> lit(listParticles); lit.ok(); ++lit){
	const RealVect& pos = lit().position();

	if(PolyGeom::dot((pos-ebCentroid), normal) >= 0.0){
	  num += lit().mass();
	}
      }
      
      numFab(iv, idx) = num;
    }
  }
}

void ito_plasma_stepper::compute_reactive_mean_energies_per_cell(EBAMRCellData& a_mean_energies){
  CH_TIME("ito_plasma_stepper::compute_mean-energies_per_cell(EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_reactive_particles_per_cell(EBAMRCellData)" << endl;
  }

  DataOps::setValue(a_mean_energies, 0.0);
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    this->compute_reactive_mean_energies_per_cell(*a_mean_energies[lvl], lvl);
  }
}

void ito_plasma_stepper::compute_reactive_mean_energies_per_cell(LevelData<EBCellFAB>& a_mean_energies, const int a_level){
  CH_TIME("ito_plasma_stepper::compute_reactive_mean_energies_per_cell(ppc, lvl)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_reactive_mean_energies_per_cell(ppc, lvl)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[a_level];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      
    const Box box = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];
      
    this->compute_reactive_mean_energies_per_cell(a_mean_energies[dit()], a_level, dit(), box, ebisbox);
  }
}

void ito_plasma_stepper::compute_reactive_mean_energies_per_cell(EBCellFAB&      a_mean_energies,
								 const int       a_level,
								 const DataIndex a_dit,
								 const Box       a_box,
								 const EBISBox&  a_ebisbox){
  CH_TIME("ito_plasma_stepper::compute_reactive_mean_energies_per_cell(ppc, lvl, dit, box)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_reactive_mean_energies_per_cell(ppc, lvl, dit, box)" << endl;
  }

  BaseFab<Real>& numFab = a_mean_energies.getSingleValuedFAB();

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>& solver = solver_it();
    const int idx                     = solver_it.index();

    const ParticleContainer<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::bulk);
    const BinFab<ItoParticle>& cellParticles         = particles.getCellParticles(a_level, a_dit);

    // Regular cells. 
    for (BoxIterator bit(a_box); bit.ok(); ++bit){
      const IntVect iv = bit();

      if (a_ebisbox.isRegular(iv)){
	Real m = 0.0;
	Real E = 0.0;
	
	const List<ItoParticle>& listParticles = cellParticles(iv, 0);

	for (ListIterator<ItoParticle> lit(listParticles); lit.ok(); ++lit){
	  m += lit().mass();
	  E += lit().mass()*lit().energy();
	}

	numFab(iv, idx) = E/m;
      }
    }

    // Irregular cells.
    VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_phase)[a_level])[a_dit];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof       = vofit();
      const IntVect iv          = vof.gridIndex();
      const RealVect normal     = a_ebisbox.normal(vof);
      const RealVect ebCentroid = a_ebisbox.bndryCentroid(vof);

      Real m = 0.0;
      Real E = 0.0;

      const List<ItoParticle>& listParticles = cellParticles(iv, 0);
      
      for (ListIterator<ItoParticle> lit(listParticles); lit.ok(); ++lit){
	const RealVect& pos = lit().position();

	if(PolyGeom::dot((pos-ebCentroid), normal) >= 0.0){
	  m += lit().mass();
	  E += lit().mass()*lit().energy();
	}
      }

      numFab(iv, idx) = E/m;
    }
  }
}

void ito_plasma_stepper::advanceReactionNetwork_nwo(const Real a_dt){
  CH_TIME("ito_plasma_stepper::advanceReactionNetwork_nwo(dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advanceReactionNetwork_nwo(dt)" << endl;
  }

  this->advanceReactionNetwork_nwo(m_fluid_E, m_EdotJ, a_dt);
}

void ito_plasma_stepper::advanceReactionNetwork_nwo(const EBAMRCellData& a_E, const EBAMRCellData& a_EdotJ, const Real a_dt){
  CH_TIME("ito_plasma_stepper::advanceReactionNetwork(ppc, ypc, E, sources, dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advanceReactionNetwork(ppc, ypc, E, sources, dt)" << endl;
  }

  // 1. Compute the number of particles per cell. Set the number of Photons to be generated per cell to zero. 
  this->compute_reactive_particles_per_cell(m_particle_ppc);
  this->compute_reactive_mean_energies_per_cell(m_particle_eps);

  m_fluid_ppc.copy(m_particle_ppc);
  m_fluid_eps.copy(m_particle_eps);
  
  DataOps::setValue(m_fluid_ypc,    0.0);
  DataOps::setValue(m_particle_ypc, 0.0);
  DataOps::copy(m_particle_old, m_particle_ppc);

  // 2. Solve for the new number of particles per cell. This also obtains the number of Photons to be generated in each cell. 
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    this->advanceReactionNetwork_nwo(*m_fluid_ppc[lvl], *m_fluid_ypc[lvl], *m_fluid_eps[lvl], *a_E[lvl], *a_EdotJ[lvl], lvl, a_dt);
  }

  // 3. Copy the results to the particle Realm.
  m_particle_ppc.copy(m_fluid_ppc);
  m_particle_ypc.copy(m_fluid_ypc);
  m_particle_eps.copy(m_fluid_eps);

  // 4. Reconcile particles on the particle Realm. Not implemented (yet).
  this->reconcileParticles(m_particle_ppc, m_particle_old, m_particle_eps, m_particle_ypc);
}


void ito_plasma_stepper::advanceReactionNetwork_nwo(LevelData<EBCellFAB>&       a_particlesPerCell,
						      LevelData<EBCellFAB>&       a_newPhotonsPerCell,
						      LevelData<EBCellFAB>&       a_meanParticleEnergies,
						      const LevelData<EBCellFAB>& a_E,
						      const LevelData<EBCellFAB>& a_EdotJ,
						      const int                   a_level,
						      const Real                  a_dt){
  CH_TIME("ito_plasma_stepper::advanceReactionNetwork(ppc, ypc, energies, E, sources, level, dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advanceReactionNetwork(ppc, ypc, energies, E, sources, level, dt)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluid_Realm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    this->advanceReactionNetwork_nwo(a_particlesPerCell[dit()],
				       a_newPhotonsPerCell[dit()],
				       a_meanParticleEnergies[dit()],
				       a_E[dit()],
				       a_EdotJ[dit()],
				       a_level,
				       dit(),
				       dbl[dit()],
				       m_amr->getDx()[a_level],
				       a_dt);
  }
}

void ito_plasma_stepper::advanceReactionNetwork_nwo(EBCellFAB&       a_particlesPerCell,
						      EBCellFAB&       a_newPhotonsPerCell,
						      EBCellFAB&       a_meanParticleEnergies,
						      const EBCellFAB& a_E,
						      const EBCellFAB& a_EdotJ,
						      const int        a_level,
						      const DataIndex  a_dit,
						      const Box        a_box,
						      const Real       a_dx,
						      const Real       a_dt){
  CH_TIME("ito_plasma_stepper::advanceReactionNetwork_nwo(ppc, ypc, E, sources, level, dit, box, dx, dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advanceReactionNetwork_nwo(ppc, ypc, E, sources, level, dit, box, dx, dt)" << endl;
  }
  
  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_RtSpecies = m_physics->getNumRtSpecies();

  const RealVect prob_lo = m_amr->getProbLo();

  const EBISBox& ebisbox = m_amr->getEBISLayout(m_fluid_Realm, m_phase)[a_level][a_dit];
  const EBISBox& ebgraph = m_amr->getEBISLayout(m_fluid_Realm, m_phase)[a_level][a_dit];

  const BaseFab<Real>& Efab = a_E.getSingleValuedFAB();

  Vector<long long> particles(num_ItoSpecies);
  Vector<long long> newPhotons(num_RtSpecies);
  Vector<Real>      meanEnergies(num_ItoSpecies);
  Vector<Real>      energySources(num_ItoSpecies);

  // Regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();

    if(ebisbox.isRegular(iv)){
      const Real kappa = 1.0;
      const Real dV    = kappa*pow(a_dx, SpaceDim);

      const RealVect pos = prob_lo + a_dx*(RealVect(iv) + 0.5*RealVect::Unit);
      const RealVect E   = RealVect(D_DECL(Efab(iv, 0), Efab(iv, 1), Efab(iv, 2)));

      // Initialize for this cell. 
      for (int i = 0; i < num_ItoSpecies; i++){
	particles[i]     = llround(a_particlesPerCell.getSingleValuedFAB()(iv, i));
	meanEnergies[i]  = a_meanParticleEnergies.getSingleValuedFAB()(iv,i);
	energySources[i] = a_EdotJ.getSingleValuedFAB()(iv, i)*dV/Units::Qe;
      }

      for (int i = 0; i < num_RtSpecies; i++){
	newPhotons[i]= 0LL;
      }
	   
      // Do the physics advance
      m_physics->advanceParticles(particles, newPhotons, meanEnergies, energySources, a_dt, E, a_dx, kappa);

      // Set result
      for (int i = 0; i < num_ItoSpecies; i++){
	a_particlesPerCell.getSingleValuedFAB()(iv, i)     = 1.0*particles[i];
	a_meanParticleEnergies.getSingleValuedFAB()(iv, i) = 1.0*meanEnergies[i];
      }

      for (int i = 0; i < num_RtSpecies; i++){
	a_newPhotonsPerCell.getSingleValuedFAB()(iv, i) = 1.0*newPhotons[i];
      }
    }
  }

  // Irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_fluid_Realm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const Real kappa    = ebisbox.volFrac(vof);
    const Real dV       = kappa*pow(a_dx, SpaceDim);
    const RealVect pos  = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, prob_lo);
    const RealVect E    = RealVect(D_DECL(a_E(vof,0), a_E(vof,1), a_E(vof,2)));


    // Initialize for this cell. 
    for (int i = 0; i < num_ItoSpecies; i++){
      particles[i]     = a_particlesPerCell(vof, i);
      meanEnergies[i]  = a_meanParticleEnergies(vof, i);
      energySources[i] = a_EdotJ(vof, i)*dV/Units::Qe;
    }

    for (int i = 0; i < num_RtSpecies; i++){
      newPhotons[i]= 0LL;
    }

    m_physics->advanceParticles(particles, newPhotons, meanEnergies, energySources, a_dt, E, a_dx, kappa);


    // Set result
    for (int i = 0; i < num_ItoSpecies; i++){
      a_particlesPerCell(vof, i)     = 1.0*particles[i];
      a_meanParticleEnergies(vof, i) = 1.0*meanEnergies[i];
    }
    
    for (int i = 0; i < num_RtSpecies; i++){
      a_newPhotonsPerCell(vof, i) = 1.0*newPhotons[i];
    }
  }
}

void ito_plasma_stepper::reconcileParticles(const EBAMRCellData& a_newParticlesPerCell,
					     const EBAMRCellData& a_oldParticlesPerCell,
					     const EBAMRCellData& a_meanParticleEnergies,
					     const EBAMRCellData& a_newPhotonsPerCell){
  CH_TIME("ito_plasma_stepper::reconcileParticles(EBAMRCellData, EBAMRCellData, EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::reconcileParticles(EBAMRCellData, EBAMRCellData, EBAMRCellData)";
  }

  for(int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){

    this->reconcileParticles(*a_newParticlesPerCell[lvl], *a_oldParticlesPerCell[lvl], *a_meanParticleEnergies[lvl], *a_newPhotonsPerCell[lvl], lvl);
  }
}

void ito_plasma_stepper::reconcileParticles(const LevelData<EBCellFAB>& a_newParticlesPerCell,
					     const LevelData<EBCellFAB>& a_oldParticlesPerCell,
					     const LevelData<EBCellFAB>& a_meanParticleEnergies,
					     const LevelData<EBCellFAB>& a_newPhotonsPerCell,
					     const int                   a_level){
  CH_TIME("ito_plasma_stepper::reconcileParticles(LevelData<EBCellFAB>x3, int)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::reconcileParticles(LevelData<EBCellFAB>x3, int)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    this->reconcileParticles(a_newParticlesPerCell[dit()],
			      a_oldParticlesPerCell[dit()],
			      a_meanParticleEnergies[dit()],
			      a_newPhotonsPerCell[dit()],
			      a_level,
			      dit(),
			      dbl[dit()],
			      m_amr->getDx()[a_level]);
  }
}

void ito_plasma_stepper::reconcileParticles(const EBCellFAB& a_newParticlesPerCell,
					     const EBCellFAB& a_oldParticlesPerCell,
					     const EBCellFAB& a_meanParticleEnergies,
					     const EBCellFAB& a_newPhotonsPerCell,
					     const int        a_level,
					     const DataIndex  a_dit,
					     const Box        a_box,
					     const Real       a_dx){
  CH_TIME("ito_plasma_stepper::reconcileParticles(EBCellFABx3, int, DataIndex, Box, Real)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::reconcileParticles(EBCellFABx3, int, DataIndex, Box, Real)" << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_RtSpecies = m_physics->getNumRtSpecies();

  const RealVect prob_lo = m_amr->getProbLo();

  const EBISBox& ebisbox = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_level][a_dit];
  const EBISBox& ebgraph = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_level][a_dit];

  Vector<BinFab<ItoParticle>* > particlesFAB(num_ItoSpecies);
  Vector<BinFab<Photon>* >       sourcePhotonsFAB(num_RtSpecies);
  Vector<BinFab<Photon>* >       bulkPhotonsFAB(num_RtSpecies);


  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>& solver = solver_it();
    const int idx = solver_it.index();
    
    ParticleContainer<ItoParticle>& solverParticles = solver->getParticles(ItoSolver::WhichContainer::bulk);
    
    particlesFAB[idx] = &(solverParticles.getCellParticles(a_level, a_dit));
  }
  
  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<McPhoto>& solver = solver_it();
    const int idx = solver_it.index();
    
    ParticleContainer<Photon>& solverBulkPhotons = solver->getBulkPhotons();
    ParticleContainer<Photon>& solverSourPhotons = solver->getSourcePhotons();
    
    bulkPhotonsFAB[idx]   = &(solverBulkPhotons.getCellParticles(a_level, a_dit));
    sourcePhotonsFAB[idx] = &(solverSourPhotons.getCellParticles(a_level, a_dit));
  }

  // Regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    if(ebisbox.isRegular(iv)){
      const RealVect cellPos       = prob_lo + a_dx*(RealVect(iv) + 0.5*RealVect::Unit);
      const RealVect centroidPos   = cellPos;
      const RealVect lo            = -0.5*RealVect::Unit;
      const RealVect hi            = 0.5*RealVect::Unit;
      const RealVect bndryCentroid = RealVect::Zero;
      const RealVect bndryNormal   = RealVect::Zero;
      const Real     kappa         = 1.0;

      Vector<List<ItoParticle>* > particles(num_ItoSpecies);
      Vector<List<Photon>* >       bulkPhotons(num_RtSpecies);
      Vector<List<Photon>* >       sourcePhotons(num_RtSpecies);
      Vector<RefCountedPtr<RtSpecies> > photoSpecies(num_RtSpecies);

      Vector<Real>      particleMeanEnergies(num_ItoSpecies);
      Vector<long long> numNewParticles(num_ItoSpecies);
      Vector<long long> numOldParticles(num_ItoSpecies);
      Vector<long long> numNewPhotons(num_RtSpecies);

      for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();

	particles[idx]            = &((*particlesFAB[idx])(iv, 0));
	particleMeanEnergies[idx] = a_meanParticleEnergies.getSingleValuedFAB()(iv, idx);
	numNewParticles[idx]      = llround(a_newParticlesPerCell.getSingleValuedFAB()(iv, idx));
	numOldParticles[idx]      = llround(a_oldParticlesPerCell.getSingleValuedFAB()(iv, idx));
      }

      for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();

	bulkPhotons[idx]   = &((*bulkPhotonsFAB[idx])(iv, 0));
	sourcePhotons[idx] = &((*sourcePhotonsFAB[idx])(iv, 0));
	photoSpecies[idx]  = solver_it()->getSpecies();
	numNewPhotons[idx] = llround(a_newPhotonsPerCell.getSingleValuedFAB()(iv, idx));

	sourcePhotons[idx]->clear();
      }

      // Reconcile particles, Photons, and photoionization
      m_physics->reconcileParticles(particles, numNewParticles, numOldParticles, cellPos, centroidPos, lo, hi, bndryCentroid, bndryNormal, a_dx, kappa);
      m_physics->reconcilePhotons(sourcePhotons, numNewPhotons, cellPos, centroidPos, lo, hi, bndryCentroid, bndryNormal, a_dx, kappa);
      m_physics->reconcilePhotoionization(particles, particleMeanEnergies, numNewParticles, bulkPhotons);
      m_physics->setMeanParticleEnergy(particles, particleMeanEnergies);
      
      // Clear the bulk Photons - they have now been absorbed on the mesh. 
      for (int i = 0; i < num_RtSpecies; i++){
	//	bulkPhotons[i]->clear();
      }
    }
  }

  // Irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof          = vofit();
    const IntVect iv             = vof.gridIndex();
    const Real kappa             = ebisbox.volFrac(vof);
    const RealVect cellPos       = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, prob_lo);
    const RealVect centroidPos   = ebisbox.centroid(vof);
    const RealVect bndryCentroid = ebisbox.bndryCentroid(vof);
    const RealVect bndryNormal   = ebisbox.normal(vof);

    RealVect lo            = -0.5*RealVect::Unit;
    RealVect hi            = 0.5*RealVect::Unit;
    if(kappa < 1.0){
      DataOps::computeMinValidBox(lo, hi, bndryNormal, bndryCentroid);
    }

    Vector<List<ItoParticle>* > particles(num_ItoSpecies);
    Vector<List<Photon>* >       bulkPhotons(num_RtSpecies);
    Vector<List<Photon>* >       sourcePhotons(num_RtSpecies);
    Vector<RefCountedPtr<RtSpecies> > photoSpecies(num_RtSpecies);

    Vector<Real>      particleMeanEnergies(num_ItoSpecies);
    Vector<long long> numNewParticles(num_ItoSpecies);
    Vector<long long> numOldParticles(num_ItoSpecies);
    Vector<long long> numNewPhotons(num_RtSpecies);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();

      particles[idx]            = &((*particlesFAB[idx])(iv, 0));
      particleMeanEnergies[idx] = a_meanParticleEnergies(vof, idx);
      numNewParticles[idx]      = llround(a_newParticlesPerCell(vof, idx));
      numOldParticles[idx]      = llround(a_oldParticlesPerCell(vof, idx));
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();

      bulkPhotons[idx]   = &((*bulkPhotonsFAB[idx])(iv, 0));
      sourcePhotons[idx] = &((*sourcePhotonsFAB[idx])(iv, 0));
      photoSpecies[idx]  = solver_it()->getSpecies();
      numNewPhotons[idx] = llround(a_newPhotonsPerCell(vof, idx));

      sourcePhotons[idx]->clear();
    }

    // Reconcile particles, Photons, and photoionization
    m_physics->reconcileParticles(particles, numNewParticles, numOldParticles, cellPos, centroidPos, lo, hi, bndryCentroid, bndryNormal, a_dx, kappa);
    m_physics->reconcilePhotons(sourcePhotons, numNewPhotons, cellPos, centroidPos, lo, hi, bndryCentroid, bndryNormal, a_dx, kappa);
    m_physics->reconcilePhotoionization(particles, particleMeanEnergies, numNewParticles, bulkPhotons);
    m_physics->setMeanParticleEnergy(particles, particleMeanEnergies);
      
    // Clear the bulk Photons - they have now been absorbed on the mesh. 
    for (int i = 0; i < num_RtSpecies; i++){
      //      bulkPhotons[i]->clear();
    }
  }
}

void ito_plasma_stepper::advanceReactionNetwork(const Real a_dt){
  CH_TIME("ito_plasma_stepper::advanceReactionNetwork(a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advanceReactionNetwork(a_dt)" << endl;
  }

  if(m_nwo_reactions){
    this->advanceReactionNetwork_nwo(a_dt);
  }
  else{
    const int num_ItoSpecies = m_physics->getNumItoSpecies();
    const int num_RtSpecies = m_physics->getNumRtSpecies();
  
    Vector<ParticleContainer<ItoParticle>* > particles(num_ItoSpecies);  // Current particles. 
    Vector<ParticleContainer<Photon>* > bulk_Photons(num_RtSpecies);     // Photons absorbed on mesh
    Vector<ParticleContainer<Photon>* > new_Photons(num_RtSpecies);      // Produced Photons go here.

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      particles[solver_it.index()] = &(solver_it()->getParticles(ItoSolver::WhichContainer::bulk));
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      bulk_Photons[solver_it.index()] = &(solver_it()->getBulkPhotons());
      new_Photons[solver_it.index()] = &(solver_it()->getSourcePhotons());
    }

    this->advanceReactionNetwork(particles, bulk_Photons, new_Photons, m_energy_sources, m_particle_E, a_dt);
  }
}

void ito_plasma_stepper::advanceReactionNetwork(Vector<ParticleContainer<ItoParticle>* >& a_particles,
						  Vector<ParticleContainer<Photon>* >&       a_Photons,
						  Vector<ParticleContainer<Photon>* >&       a_newPhotons,
						  const Vector<EBAMRCellData>&                a_sources,
						  const EBAMRCellData&                        a_E,
						  const Real                                  a_dt){
  CH_TIME("ito_plasma_stepper::advanceReactionNetwork(Vector<ParticleContainer*> x3, Vector<EBAMRCellData*>, EBAMRCellData, Real)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advanceReactionNetwork(Vector<ParticleContainer*> x3, Vector<EBAMRCellData*>, EBAMRCellData, Real)" << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_RtSpecies = m_physics->getNumRtSpecies();

  Vector<AMRCellParticles<ItoParticle>* > particles(num_ItoSpecies);
  Vector<AMRCellParticles<Photon>* >       Photons(num_ItoSpecies);
  Vector<AMRCellParticles<Photon>* >       newPhotons(num_ItoSpecies);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    particles[idx] = &(a_particles[idx]->getCellParticles());
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    Photons[idx]    = &(a_Photons[idx]->getCellParticles());
    newPhotons[idx] = &(a_newPhotons[idx]->getCellParticles());
  }
				
  //Advance reaction network
  this->advanceReactionNetwork(particles, Photons, newPhotons, a_sources, m_particle_E, a_dt);

}

void ito_plasma_stepper::advanceReactionNetwork(Vector<AMRCellParticles<ItoParticle>* >& a_particles,
						  Vector<AMRCellParticles<Photon>* >&       a_Photons,
						  Vector<AMRCellParticles<Photon>* >&       a_newPhotons,
						  const Vector<EBAMRCellData>&              a_sources,
						  const EBAMRCellData&                      a_E,
						  const Real                                a_dt){
  CH_TIME("ito_plasma_stepper::advanceReactionNetwork(Vector<AMRCellParticles*> x 3, Vector<EBAMRCellData*>, EBAMRCellData, Real)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advanceReactionNetwork(Vector<AMRCellParticles*> x 3, Vector<EBAMRCellData*>, EBAMRCellData, Real)" << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_RtSpecies = m_physics->getNumRtSpecies();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    Vector<LayoutData<BinFab<ItoParticle> >* > particles(num_ItoSpecies);
    Vector<LayoutData<BinFab<Photon> >* >       Photons(num_RtSpecies);
    Vector<LayoutData<BinFab<Photon> >* >       newPhotons(num_RtSpecies);
    Vector<LevelData<EBCellFAB>* >              sources(num_ItoSpecies);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      particles[idx] = &(*(*a_particles[idx])[lvl]);
      sources[idx]   = &(*(a_sources[idx])[lvl]);
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      Photons[idx]    = &(*(*a_Photons[idx])[lvl]);
      newPhotons[idx] = &(*(*a_newPhotons[idx])[lvl]);
    }
      
    this->advanceReactionNetwork(particles, Photons, newPhotons, sources, *a_E[lvl], lvl, a_dt);
  }
}

void ito_plasma_stepper::advanceReactionNetwork(Vector<LayoutData<BinFab<ItoParticle> >* >& a_particles,
						  Vector<LayoutData<BinFab<Photon> >* >&       a_Photons,
						  Vector<LayoutData<BinFab<Photon> >* >&       a_newPhotons,
						  const Vector<LevelData<EBCellFAB>* >&        a_sources,
						  const LevelData<EBCellFAB>&                  a_E,
						  const int                                    a_lvl,
						  const Real                                   a_dt){
  CH_TIME("ito_plasma_stepper::advanceReactionNetwork(Vector<LD<BinFab>* > x 3, Vector<LD<EBCellFAB>*>, EBAMRCellData, level, dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advanceReactionNetwork(Vector<LD<BinFab>* > x 3, Vector<LD<EBCellFAB>*>, EBAMRCellData, level, dt)" << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_RtSpecies = m_physics->getNumRtSpecies();

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[a_lvl];
  const Real dx = m_amr->getDx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());

    Vector<BinFab<ItoParticle>* > particles(num_ItoSpecies);
    Vector<BinFab<Photon>* >       Photons(num_RtSpecies);;
    Vector<BinFab<Photon>* >       newPhotons(num_RtSpecies);
    Vector<EBCellFAB*>             sources(num_ItoSpecies);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      particles[idx] = &((*a_particles[idx])[dit()]);
      sources[idx]   = &((*a_sources[idx])[dit()]);
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      Photons[idx]    = &((*a_Photons[idx])[dit()]);
      newPhotons[idx] = &((*a_newPhotons[idx])[dit()]);
    }

    this->advanceReactionNetwork(particles, Photons, newPhotons, sources, a_E[dit()], a_lvl, dit(), box, dx, a_dt);
  }
}

void ito_plasma_stepper::advanceReactionNetwork(Vector<BinFab<ItoParticle>* >& a_particles,
						  Vector<BinFab<Photon>* >&       a_Photons,
						  Vector<BinFab<Photon>* >&       a_newPhotons,
						  const Vector<EBCellFAB*>&       a_sources,
						  const EBCellFAB&                a_E,
						  const int                       a_lvl,
						  const DataIndex                 a_dit,
						  const Box                       a_box,
						  const Real                      a_dx,
						  const Real                      a_dt){
  CH_TIME("ito_plasma_stepper::advanceReactionNetwork(Vector<BinFab*> x 3, Vector<EBCellFAB*>, EBCellFAB, level, dit, box, dx, dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advanceReactionNetwork(Vector<BinFab*> x 3, Vector<EBCellFAB*>, EBCellFAB, level, dit, box, dx, dt)" << endl;
  }

  const int comp = 0;

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_RtSpecies = m_physics->getNumRtSpecies();

  const RealVect prob_lo = m_amr->getProbLo();
  const RealVect dx      = a_dx*RealVect::Unit;

  const EBISBox& ebisbox = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_lvl][a_dit];
  const EBISBox& ebgraph = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_lvl][a_dit];

  const BaseFab<Real>& Efab = a_E.getSingleValuedFAB();

  // Regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    
    if(ebisbox.isRegular(iv)){
      const Real kappa   = 1.0;
      const RealVect pos = prob_lo + a_dx*(RealVect(iv) + 0.5*RealVect::Unit);
      const RealVect e   = RealVect(D_DECL(Efab(iv, 0), Efab(iv, 1), Efab(iv, 2)));
      
      Vector<List<ItoParticle>* > particles(num_ItoSpecies);
      Vector<List<Photon>* >       Photons(num_RtSpecies);
      Vector<List<Photon>* >       newPhotons(num_RtSpecies);
      Vector<Real>                 sources(num_ItoSpecies);

      for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
      
	List<ItoParticle>& bp = (*a_particles[idx])(iv, comp);
	particles[idx] = &bp;

	const BaseFab<Real>& sourcesFAB = a_sources[idx]->getSingleValuedFAB();
	sources[idx] = sourcesFAB(iv, comp);
      }

      for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
      
	List<Photon>& bp    = (*a_Photons[idx])(iv, comp);
	List<Photon>& bpNew = (*a_newPhotons[idx])(iv, comp);
      
	Photons[idx]    = &bp;
	newPhotons[idx] = &bpNew;
      }

      // Dummy stuff for regular cells
      const RealVect lo = -0.5*RealVect::Unit;
      const RealVect hi =  0.5*RealVect::Unit;
      const RealVect n  = RealVect::Zero;
      const RealVect c  = RealVect::Zero;

      // Advance reactions
      m_physics->advanceReactionNetwork(particles, Photons, newPhotons, sources, e, pos, c, c, n, lo, hi, a_dx, kappa, a_dt);
    }
  }

  // Now do the irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_phase)[a_lvl])[a_dit];
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
      DataOps::computeMinValidBox(lo, hi, n, ebc);
    }

    Vector<List<ItoParticle>* > particles(num_ItoSpecies);
    Vector<List<Photon>* >       Photons(num_RtSpecies);
    Vector<List<Photon>* >       newPhotons(num_RtSpecies);
    Vector<Real>                 sources(num_ItoSpecies);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      
      List<ItoParticle>& bp = (*a_particles[idx])(iv, comp);
      particles[idx] = &bp;

      const BaseFab<Real>& sourcesFAB = a_sources[idx]->getSingleValuedFAB();
      sources[idx] = sourcesFAB(iv, comp);
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      
      List<Photon>& bp    = (*a_Photons[idx])(iv, comp);
      List<Photon>& bpNew = (*a_newPhotons[idx])(iv, comp);
      
      Photons[idx]    = &bp;
      newPhotons[idx] = &bpNew;
    }

    // Advance reactions
    m_physics->advanceReactionNetwork(particles, Photons, newPhotons, sources, e, pos, cen, ebc, n, lo, hi, a_dx, kappa, a_dt);
  }
}

Real ito_plasma_stepper::compute_physics_dt() const{
  CH_TIME("ito_plasma_stepper::compute_physics_dt()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_physics_dt()" << endl;
  }

  // TLDR: This is done on the particle Realm because of the densities (which are defined on the particle Realm). 

  const Real dt = this->compute_physics_dt(m_particle_E, m_ito->getDensities());

  return dt;
}

Real ito_plasma_stepper::compute_physics_dt(const EBAMRCellData& a_E, const Vector<EBAMRCellData*> a_densities) const {
  CH_TIME("ito_plasma_stepper::compute_physics_dt(EBAMRCellFAB, Vector<EBAMRCellFAB*>)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_physics_dt(EBAMRCellFAB, Vector<EBAMRCellFAB*>)" << endl;
  }



  const int num_ItoSpecies = m_physics->getNumItoSpecies();

  Real minDt = 1.E99;
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    
    Vector<LevelData<EBCellFAB>*> densities(num_ItoSpecies);
    
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();

      densities[idx] = &(*(*a_densities[idx])[lvl]);
    }

    const Real levelDt = this->compute_physics_dt(*a_E[lvl], densities, lvl);

    minDt = Min(minDt, levelDt);
  }

  return minDt;
}

Real ito_plasma_stepper::compute_physics_dt(const LevelData<EBCellFAB>& a_E, const Vector<LevelData<EBCellFAB> *> a_densities, const int a_level) const {
  CH_TIME("ito_plasma_stepper::compute_physics_dt(LevelData<EBCellFAB>, Vector<LevelData<EBCellFAB> *>, int)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_physics_dt(LevelData<EBCellFAB>, Vector<LevelData<EBCellFAB> *>, int)" << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[a_level];

  Real minDt = 1.E99;
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

    Vector<EBCellFAB*> densities(num_ItoSpecies);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      
      densities[idx] = &((*a_densities[idx])[dit()]);
    }

    const Real boxDt = this->compute_physics_dt(a_E[dit()], densities, a_level, dit(), dbl.get(dit()));

    minDt = Min(minDt, boxDt);
  }

  // MPI reduction....
#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&minDt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("ito_plasma_stepper::compute_physics_dt(lvl) - communication error on norm");
  }
  minDt = tmp;
#endif  

  return minDt;
}

Real ito_plasma_stepper::compute_physics_dt(const EBCellFAB& a_E, const Vector<EBCellFAB*> a_densities, const int a_level, const DataIndex a_dit, const Box a_box) const {
  CH_TIME("ito_plasma_stepper::compute_physics_dt(EBCellFAB, Vector<EBCellFAB*>, lvl, dit, box)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_physics_dt(EBCellFAB, Vector<EBCellFAB*>, lvl, dit, box)" << endl;
  }

  Real minDt = 1.E99;

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  
  const int comp         = 0;
  const Real dx          = m_amr->getDx()[a_level];
  const RealVect prob_lo = m_amr->getProbLo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_level][a_dit];

  // Do regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv   = bit();
    const RealVect pos = m_amr->getProbLo() + dx*(RealVect(iv) + 0.5*RealVect::Unit);
    const RealVect e   = RealVect(D_DECL(E(iv,0), E(iv, 1), E(iv, 2)));

    if(ebisbox.isRegular(iv)){
      Vector<Real> densities(num_ItoSpecies);
      for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.index();
	
	const BaseFab<Real>& basefab = a_densities[idx]->getSingleValuedFAB();
	densities[idx] = basefab(iv, comp);
      }

      const Real cellDt = m_physics->computeDt(e, pos, densities);

      minDt = Min(minDt, cellDt);
    }
  }

  // Do irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect e    = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, prob_lo);

    Vector<Real> densities(num_ItoSpecies);
    
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.index();
      densities[idx] = (*a_densities[idx])(vof, comp);
    }

    const Real cellDt = m_physics->computeDt(e, pos, densities);

    minDt = Min(minDt, cellDt);
  }

  return minDt;
}

void ito_plasma_stepper::advance_Photons(const Real a_dt){
  CH_TIME("ito_plasma_stepper::advance_Photons(a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_advance_Photons(a_dt)" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<McPhoto>& solver = solver_it();
    
    // Add source Photons and move the Photons
    ParticleContainer<Photon>& Photons        = solver->getPhotons();
    ParticleContainer<Photon>& bulkPhotons    = solver->getBulkPhotons();
    ParticleContainer<Photon>& ebPhotons      = solver->getEbPhotons();
    ParticleContainer<Photon>& domainPhotons  = solver->getDomainPhotons();
    ParticleContainer<Photon>& sourcePhotons  = solver->getSourcePhotons();

    if(solver->isInstantaneous()){
      solver->clear(Photons);

      // Add source Photons
      Photons.addParticles(sourcePhotons);
      solver->clear(sourcePhotons);

      // Instantaneous advance
      solver->advancePhotonsStationary(bulkPhotons, ebPhotons, domainPhotons, Photons);
    }
    else{
      // Add source Photons
      Photons.addParticles(sourcePhotons);
      solver->clear(sourcePhotons);

      // Stationary advance
      solver->advancePhotonsTransient(bulkPhotons, ebPhotons, domainPhotons, Photons, a_dt);
    }
  }
}

void ito_plasma_stepper::sortPhotonsByCell(){
  CH_TIME("ito_plasma_stepper::sortPhotonsByCell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sortPhotonsByCell()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sortPhotonsByCell();
  }
}

void ito_plasma_stepper::sortPhotonsByPatch(){
  CH_TIME("ito_plasma_stepper::sortPhotonsByPatch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sortPhotonsByPatch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sortPhotonsByPatch();
  }
}

void ito_plasma_stepper::sortSourcePhotonsByCell(){
  CH_TIME("ito_plasma_stepper::sortSourcePhotonsByCell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sortSourcePhotonsByCell()" << endl;
  }
  
  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sortSourcePhotonsByCell();
  }
}

void ito_plasma_stepper::sortSourcePhotonsByPatch(){
  CH_TIME("ito_plasma_stepper::sortSourcePhotonsByPatch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sortSourcePhotonsByPatch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sortSourcePhotonsByPatch();
  }
}

void ito_plasma_stepper::sortBulkPhotonsByCell(){
  CH_TIME("ito_plasma_stepper::sortBulkPhotonsByCell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sortBulkPhotonsByCell()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sortBulkPhotonsByCell();
  }
}

void ito_plasma_stepper::sortBulkPhotonsByPatch(){
  CH_TIME("ito_plasma_stepper::sortBulkPhotonsByPatch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sortBulkPhotonsByPatch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sortBulkPhotonsByPatch();
  }
}

void ito_plasma_stepper::sortEbPhotonsByCell(){
  CH_TIME("ito_plasma_stepper::sortEbPhotonsByCell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sortEbPhotonsByCell()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sortEbPhotonsByCell();
  }
}

void ito_plasma_stepper::sortEbPhotonsByPatch(){
  CH_TIME("ito_plasma_stepper::sortEbPhotonsByPatch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sortEbPhotonsByPatch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sortEbPhotonsByPatch();
  }
}

void ito_plasma_stepper::sortDomainPhotonsByCell(){
  CH_TIME("ito_plasma_stepper::sortDomainPhotonsByCell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sortDomainPhotonsByCell()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sortDomainPhotonsByCell();
  }
}

void ito_plasma_stepper::sortDomainPhotonsByPatch(){
  CH_TIME("ito_plasma_stepper::sortDomainPhotonsByPatch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sortDomainPhotonsByPatch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sortDomainPhotonsByPatch();
  }

}

bool ito_plasma_stepper::loadBalanceThisRealm(const std::string a_realm) const{
  CH_TIME("TimeStepper::loadBalanceThisRealm");
  if(m_verbosity > 5){
    pout() << "TimeStepper::loadBalanceThisRealm" << endl;
  }

  bool ret = false;

  if(a_realm == m_particleRealm && m_LoadBalancing){
    ret = true;
  }
  
  return ret;
}

void ito_plasma_stepper::loadBalanceBoxes(Vector<Vector<int> >&            a_procs,
					    Vector<Vector<Box> >&            a_boxes,
					    const std::string                a_realm,
					    const Vector<DisjointBoxLayout>& a_grids,
					    const int                        a_lmin,
					    const int                        a_finestLevel){
  CH_TIME("ito_plasma_stepper::loadBalanceBoxes");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper_stepper::loadBalanceBoxes" << endl;
  }

  if(m_LoadBalancing && a_realm == m_particleRealm){
    this->LoadBalancing_particle_Realm(a_procs, a_boxes, a_realm, a_grids, a_lmin, a_finestLevel);
  }
  else{
    MayDay::Abort("ito_plasma_stepper::loadBalanceBoxes - shouldn't happen, how did you get here..?");
  }
}

Vector<long int> ito_plasma_stepper::getCheckpointLoads(const std::string a_realm, const int a_level) const {
  CH_TIME("ito_plasma_stepper::getCheckpointLoads(...)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper_stepper::getCheckpointLoads(...)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(a_realm)[a_level];
  const int nbox = dbl.size();

  Vector<long int> loads(nbox, 0L);
  if(m_LoadBalancing && a_realm == m_particleRealm){
    Vector<RefCountedPtr<ItoSolver> > lb_solvers = this->get_lb_solvers();
    
    for (int isolver = 0; isolver < lb_solvers.size(); isolver++){
      Vector<long int> solver_loads(nbox, 0L);
      lb_solvers[isolver]->computeLoads(solver_loads, dbl, a_level);

      for (int ibox = 0; ibox < nbox; ibox++){
	loads[ibox] += solver_loads[ibox];
      }
    }

    // Now add the "constant" loads
    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit){
      const Box box  = dbl[lit()];
      const int ibox = lit().intCode();

      loads[ibox] += lround(m_load_ppc*box.numPts());
    }
  }
  else{
    loads = TimeStepper::getCheckpointLoads(a_realm, a_level);
  }

  return loads;
}

void ito_plasma_stepper::LoadBalancing_particle_Realm(Vector<Vector<int> >&            a_procs,
						     Vector<Vector<Box> >&            a_boxes,
						     const std::string                a_realm,
						     const Vector<DisjointBoxLayout>& a_grids,
						     const int                        a_lmin,
						     const int                        a_finestLevel){
  CH_TIME("ito_plasma_stepper::LoadBalancing_particle_Realm(...)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper_stepper::LoadBalancing_particle_Realm(...)" << endl;
  }
  
  // Decompose the DisjointBoxLayout
  a_procs.resize(1 + a_finestLevel);
  a_boxes.resize(1 + a_finestLevel);
  
  for (int lvl = a_lmin; lvl <= a_finestLevel; lvl++){
    a_procs[lvl] = a_grids[lvl].procIDs();
    a_boxes[lvl] = a_grids[lvl].boxArray();
  }

  // Get the particles that we will use for load balancing. 
  Vector<RefCountedPtr<ItoSolver> > lb_solvers = this->get_lb_solvers();

  // Regrid particles onto the "dummy grids" a_grids
  for (int i = 0; i < lb_solvers.size(); i++){
    ParticleContainer<ItoParticle>& particles = lb_solvers[i]->getParticles(ItoSolver::WhichContainer::bulk);
    
    particles.regrid(a_grids, m_amr->getDomains(), m_amr->getDx(), m_amr->getRefinementRatios(), a_lmin, a_finestLevel);

    // If we make superparticles during regrids, do it here so we can better estimate the computational loads for each patch. This way, if a grid is removed the realistic
    // load estimate of the underlying grid(s) is improved.
    if(m_regrid_superparticles){
      particles.sortParticlesByCell();
      lb_solvers[i]->makeSuperparticles(ItoSolver::WhichContainer::bulk, m_ppc);
      particles.sortParticlesByPatch();
    }
  }

  // Get loads on each level
  Vector<Vector<long int> > loads(1 + a_finestLevel);
  for (int lvl = 0; lvl <= a_finestLevel; lvl++){
    loads[lvl] = this->getCheckpointLoads(a_realm, lvl);
  }


  // Do the actual load balancing
  LoadBalancing::sort(a_boxes, loads, m_boxSort);
  LoadBalancing::balanceLevelByLevel(a_procs, loads, a_boxes);
  //  LoadBalancing::hierarchy(a_procs, loads, a_boxes); If you want to try something crazy...

  // Go back to "pre-regrid" mode so we can get particles to the correct patches after load balancing. 
  for (int i = 0; i < lb_solvers.size(); i++){
    ParticleContainer<ItoParticle>& particles = lb_solvers[i]->getParticles(ItoSolver::WhichContainer::bulk);
    particles.preRegrid(a_lmin);
  }
}

Vector<RefCountedPtr<ItoSolver> > ito_plasma_stepper::get_lb_solvers() const {
  CH_TIME("ito_plasma_stepper::get_lb_solvers()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::get_lb_solvers()" << endl;
  }

  Vector<RefCountedPtr<ItoSolver> > lb_solvers;

  if(m_LoadBalancing_idx < 0){
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<ItoSolver>& solver = solver_it();
      
      lb_solvers.push_back(solver);
    }
  }
  else {
    RefCountedPtr<ItoSolver>& solver = m_ito->getSolvers()[m_LoadBalancing_idx];
    lb_solvers.push_back(solver);
  }

  return lb_solvers;
}

void ito_plasma_stepper::computeElectricFielddotJ_source(const Real a_dt){
  CH_TIME("ito_plasma_stepper::computeElectricFielddotJ_source(a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeElectricFielddotJ_source(a_dt)" << endl;
  }

  // Swap between these two. 
  if(m_nwo_reactions){
    //    this->computeElectricFielddotJ_source_nwo();
    this->computeElectricFielddotJ_source_nwo2(a_dt);
  }
  else{
    this->computeElectricFielddotJ_source();
  }
}

void ito_plasma_stepper::computeElectricFielddotJ_source(){
  CH_TIME("ito_plasma_stepper::computeElectricFielddotJ_source()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeElectricFielddotJ_source()" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>& solver   = solver_it();
    RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    DataOps::setValue(m_energy_sources[idx], 0.0);

    // Do mobile contribution. 
    if(q != 0 && solver->isMobile()){

      // Drift contribution
      solver->depositConductivity(m_particle_scratch1, solver->getParticles(ItoSolver::WhichContainer::bulk)); // Deposit mu*n
      DataOps::copy(m_particle_scratchD, m_particle_E); // Could use m_particle_E or solver's m_velo_func here, but m_velo_func = +/- E (depends on q)
      
      DataOps::multiplyScalar(m_particle_scratchD, m_particle_scratch1);        // m_particle_scratchD = mu*n*E
      DataOps::dotProduct(m_particle_scratch1, m_particle_E, m_particle_scratchD); // m_particle_scratch1 = mu*n*E*E
      DataOps::incr(m_energy_sources[idx], m_particle_scratch1, 1.0);            // a_source[idx] += mu*n*E*E
    }

    // Diffusive contribution
    if(q != 0 && solver->isDiffusive()){

      // Compute the negative gradient of the diffusion term
      solver->depositDiffusivity(m_particle_scratch1, solver->getParticles(ItoSolver::WhichContainer::bulk));
      m_amr->computeGradient(m_particle_scratchD, m_particle_scratch1, m_particleRealm, m_phase);
      DataOps::scale(m_particle_scratchD, -1.0); // scratchD = -grad(D*n)
      
      DataOps::dotProduct(m_particle_scratch1, m_particle_scratchD, m_particle_E); // m_particle_scratch1 = -E*grad(D*n)
      DataOps::incr(m_energy_sources[idx], m_particle_scratch1, 1.0);            // a_source[idx]
    }
    
    if (q != 0 && (solver->isMobile() || solver->isDiffusive())){
      DataOps::scale(m_energy_sources[idx], Abs(q)*Units::Qe);
    }
  }
}

void ito_plasma_stepper::computeElectricFielddotJ_source_nwo(){
  CH_TIME("ito_plasma_stepper::computeElectricFielddotJ_source_nwo()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeElectricFielddotJ_source_nwo()" << endl;
  }

  DataOps::setValue(m_EdotJ, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>& solver   = solver_it();
    RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    // Do mobile contribution. Computes Z*e*E*mu*n*E*E
    if(q != 0 && solver->isMobile()){
      solver->depositConductivity(m_particle_scratch1, solver->getParticles(ItoSolver::WhichContainer::bulk)); // Deposit mu*n
      m_fluid_scratch1.copy(m_particle_scratch1);                                 // Copy mu*n to fluid Realm
      DataOps::copy(m_fluid_scratchD, m_fluid_E);                                // m_fluid_scratchD = E
      DataOps::multiplyScalar(m_fluid_scratchD, m_fluid_scratch1);              // m_fluid_scratchD = E*mu*n
      DataOps::dotProduct(m_fluid_scratch1, m_fluid_E, m_fluid_scratchD);          // m_particle_scratch1 = E.dot.(E*mu*n)
      DataOps::scale(m_fluid_scratch1, Abs(q)*Units::Qe);                      // m_particle_scratch1 = Z*e*mu*n*E*E

      m_amr->averageDown(m_fluid_scratch1, m_fluid_Realm, m_phase);
      m_amr->interpGhost(m_fluid_scratch1, m_fluid_Realm, m_phase);
      DataOps::plus(m_EdotJ, m_fluid_scratch1, 0, idx, 1);                       // a_source[idx] += Z*e*mu*n*E*E
    }

    // Diffusive contribution. Computes -Z*e*E*grad(D*n)
    if(q != 0 && solver->isDiffusive()){
      solver->depositDiffusivity(m_particle_scratch1, solver->getParticles(ItoSolver::WhichContainer::bulk));            // Deposit D*n
      m_fluid_scratch1.copy(m_particle_scratch1);                                           // Copy D*n to fluid Realm
      m_amr->computeGradient(m_fluid_scratchD, m_fluid_scratch1, m_fluid_Realm, m_phase);  // scratchD = grad(D*n)
      DataOps::scale(m_fluid_scratchD, -1.0);                                              // scratchD = -grad(D*n)
      DataOps::dotProduct(m_fluid_scratch1,  m_fluid_scratchD, m_fluid_E);                   // scratch1 = -E.dot.grad(D*n)
      DataOps::scale(m_fluid_scratch1, Abs(q)*Units::Qe);                                // scratch1 = -Z*e*E*grad(D*n)

      m_amr->averageDown(m_fluid_scratch1, m_fluid_Realm, m_phase);
      m_amr->interpGhost(m_fluid_scratch1, m_fluid_Realm, m_phase);
      
      DataOps::plus(m_EdotJ, m_fluid_scratch1, 0, idx, 1);                                 // source  += -Z*e*E*grad(D*n)
    }
  }
}

void ito_plasma_stepper::computeElectricFielddotJ_source_nwo2(const Real a_dt){
  CH_TIME("ito_plasma_stepper::computeElectricFielddotJ_source_nwo2(a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::computeElectricFielddotJ_source_nwo2(a_dt)" << endl;
  }

  DataOps::setValue(m_EdotJ, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ItoSolver>& solver   = solver_it();
    RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    ParticleContainer<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::bulk);

    const DepositionType::Which deposition = solver->getDeposition();

    if((mobile || diffusive) && q != 0){

      // We will interpolate m_particle_E onto particle velocity vectors.
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
	const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	  const EBCellFAB& E     = (*m_particle_E[lvl])[dit()];
	  const EBISBox& ebisbox = E.getEBISBox();
	  const FArrayBox& Efab  = E.getFArrayBox();
	  const RealVect dx      = m_amr->getDx()[lvl]*RealVect::Unit;
	  const RealVect origin  = m_amr->getProbLo();
	  const Box box          = dbl[dit()];

	  List<ItoParticle>& particleList = particles[lvl][dit()].listItems();

	  // This interpolates the velocity function on to the particle velocities
	  EbParticleInterp meshInterp(box, ebisbox, dx, origin, true);
	  meshInterp.interpolateVelocity(particleList, Efab, deposition);

	  // Go through the particles and set their mass to E.dot(X^new - X^old)
	  for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit){
	    ItoParticle& p = lit();

	    const Real m         = p.mass();
	    const RealVect& v    = p.velocity(); // Actually = E(X^new)
	    const RealVect& Xnew = p.position();
	    const RealVect& Xold = p.oldPosition();

	    p.tmp()  = m;
	    p.mass() = m*PolyGeom::dot(v, Xnew-Xold);
	  } 
	}
      }

      // Deposit the result
      solver->depositParticles(m_particle_scratch1, particles, deposition);
      m_fluid_scratch1.copy(m_particle_scratch1);

      // Scale by Qe/dt to make it Joule/dt. Then add to correct index
      DataOps::scale(m_fluid_scratch1, q*Units::Qe/a_dt);
      DataOps::plus(m_EdotJ, m_fluid_scratch1, 0, idx, 1);

      // Set p.mass() back to the original value
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
	const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];
	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	  List<ItoParticle>& particleList = particles[lvl][dit()].listItems();

	  for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit){
	    ItoParticle& p = lit();
	    p.mass() = p.tmp();
	  }
	}
      }
    }
  }
}
#include "CD_NamespaceFooter.H"
