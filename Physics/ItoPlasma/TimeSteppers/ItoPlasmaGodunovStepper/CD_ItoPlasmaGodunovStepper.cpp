/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaGodunovStepper.cpp
  @brief  Implementation of CD_ItoPlasmaGodunovStepper.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_ItoPlasmaGodunovStepper.H>
#include <CD_Timer.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

ItoPlasmaGodunovStepper::ItoPlasmaGodunovStepper(RefCountedPtr<ItoPlasmaPhysics>& a_physics)
{
  m_name    = "ItoPlasmaGodunovStepper";
  m_physics = a_physics;

  m_dt_relax = 1.E99;

  ParmParse pp("ItoPlasmaGodunovStepper");
  pp.get("particle_realm", m_particleRealm);
  pp.get("profile", m_profile);
  pp.get("load_ppc", m_load_ppc);
  pp.get("nwo_reactions", m_nwo_reactions);

  m_avg_cfl = 0.0;
}

ItoPlasmaGodunovStepper::~ItoPlasmaGodunovStepper() {}

int
ItoPlasmaGodunovStepper::getNumberOfPlotVariables() const
{
  CH_TIME("ItoPlasmaGodunovStepper::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::getNumberOfPlotVariables" << endl;
  }

  int ncomp = ItoPlasmaStepper::getNumberOfPlotVariables();

  ncomp++; // Add conductivity

  return ncomp;
}

void
ItoPlasmaGodunovStepper::writePlotData(EBAMRCellData&       a_output,
                                       Vector<std::string>& a_plotVariableNames,
                                       int&                 a_icomp) const
{
  CH_TIME("ItoPlasmaGodunovStepper::writeConductivity");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::writeConductivity" << endl;
  }

  ItoPlasmaStepper::writePlotData(a_output, a_plotVariableNames, a_icomp);

  // Do conductivity
  this->writeConductivity(a_output, a_icomp);
  a_plotVariableNames.push_back("conductivity");
}

void
ItoPlasmaGodunovStepper::writeConductivity(EBAMRCellData& a_output, int& a_icomp) const
{
  CH_TIME("ItoPlasmaStepper::writeConductivity");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::writeConductivity" << endl;
  }

  const Interval src_interv(0, 0);
  const Interval dst_interv(a_icomp, a_icomp);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    if (m_conduct_cell.getRealm() == a_output.getRealm()) {
      m_conduct_cell[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else {
      m_conduct_cell[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }
  a_icomp += 1;
}

void
ItoPlasmaGodunovStepper::allocate()
{
  CH_TIME("ItoPlasmaGodunovStepper::allocate");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::allocate" << endl;
  }

  m_ito->allocateInternals();
  m_rte->allocateInternals();
  m_fieldSolver->allocateInternals();
  m_sigma->allocate();

  // Now allocate for the conductivity particles and rho^dagger particles
  const int num_ItoSpecies = m_physics->getNumItoSpecies();

  m_conductivity_particles.resize(num_ItoSpecies);
  m_rho_dagger_particles.resize(num_ItoSpecies);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const RefCountedPtr<ItoSolver>& solver = solver_it();

    const int idx = solver_it.index();

    m_conductivity_particles[idx] = new ParticleContainer<PointParticle>();
    m_rho_dagger_particles[idx]   = new ParticleContainer<PointParticle>();

    m_amr->allocate(*m_conductivity_particles[idx], m_particleRealm);
    m_amr->allocate(*m_rho_dagger_particles[idx], m_particleRealm);
  }
}

void
ItoPlasmaGodunovStepper::parseOptions()
{
  CH_TIME("ItoPlasmaGodunovStepper::parseOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseOptions" << endl;
  }

  ParmParse   pp(m_name.c_str());
  std::string str;

  pp.get("verbosity", m_verbosity);
  pp.get("ppc", m_ppc);
  pp.get("max_cells_hop", m_max_cells_hop);
  pp.get("merge_interval", m_merge_interval);
  pp.get("relax_factor", m_relax_factor);
  pp.get("regrid_super", m_regridSuperparticles);
  pp.get("algorithm", str);
  pp.get("load_balance", m_LoadBalancing);
  pp.get("load_index", m_LoadBalancing_idx);
  pp.get("min_dt", m_min_dt);
  pp.get("max_dt", m_max_dt);
  pp.get("filter_rho", m_filter_rho);
  pp.get("filter_cond", m_filter_cond);
  pp.get("eb_tolerance", m_eb_tolerance);

  // Get algorithm
  if (str == "euler_maruyama") {
    m_algorithm = which_algorithm::euler_maruyama;
  }
  else if (str == "trapezoidal") {
    m_algorithm = which_algorithm::trapezoidal;
  }
  else {
    MayDay::Abort("ItoPlasmaGodunovStepper::parseOptions - unknown algorithm requested");
  }

  // Dt limitation
  pp.get("which_dt", str);
  if (str == "advection") {
    m_whichDt = which_dt::advection;
  }
  else if (str == "diffusion") {
    m_whichDt = which_dt::diffusion;
  }
  else if (str == "AdvectionDiffusion") {
    m_whichDt = which_dt::AdvectionDiffusion;
  }
  else {
    MayDay::Abort("ItoPlasmaGodunovStepper::parseOptions - unknown 'which_dt' requested");
  }

  // Box sorting for load balancing
  pp.get("box_sorting", str);
  if (str == "none") {
    m_boxSort = BoxSorting::None;
  }
  else if (str == "std") {
    m_boxSort = BoxSorting::Std;
  }
  else if (str == "shuffle") {
    m_boxSort = BoxSorting::Shuffle;
  }
  else if (str == "morton") {
    m_boxSort = BoxSorting::Morton;
  }
  else {
    MayDay::Abort(
      "ItoPlasmaGodunovStepper::parseOptions - unknown box sorting method requested for argument 'BoxSorting'");
  }

  // Parse filterse
  this->parseFilters();

  // Setup runtime storage (requirements change with algorithm)
  this->setupRuntimeStorage();
}

void
ItoPlasmaGodunovStepper::parseRuntimeOptions()
{
  CH_TIME("ItoPlasmaGodunovStepper::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }

  ParmParse   pp(m_name.c_str());
  std::string str;

  pp.get("verbosity", m_verbosity);
  pp.get("ppc", m_ppc);
  pp.get("max_cells_hop", m_max_cells_hop);
  pp.get("merge_interval", m_merge_interval);
  pp.get("relax_factor", m_relax_factor);
  pp.get("regrid_super", m_regridSuperparticles);
  pp.get("algorithm", str);
  pp.get("load_balance", m_LoadBalancing);
  pp.get("load_index", m_LoadBalancing_idx);
  pp.get("min_dt", m_min_dt);
  pp.get("max_dt", m_max_dt);
  pp.get("filter_rho", m_filter_rho);
  pp.get("filter_cond", m_filter_cond);
  pp.get("eb_tolerance", m_eb_tolerance);

  // Get algorithm
  if (str == "euler_maruyama") {
    m_algorithm = which_algorithm::euler_maruyama;
  }
  else if (str == "trapezoidal") {
    m_algorithm = which_algorithm::trapezoidal;
  }
  else {
    MayDay::Abort("ItoPlasmaGodunovStepper::parseOptions - unknown algorithm requested");
  }

  // Dt limitation
  pp.get("which_dt", str);
  if (str == "advection") {
    m_whichDt = which_dt::advection;
  }
  else if (str == "diffusion") {
    m_whichDt = which_dt::diffusion;
  }
  else if (str == "AdvectionDiffusion") {
    m_whichDt = which_dt::AdvectionDiffusion;
  }
  else {
    MayDay::Abort("ItoPlasmaGodunovStepper::parseOptions - unknown 'which_dt' requested");
  }

  // Box sorting for load balancing
  pp.get("box_sorting", str);
  if (str == "none") {
    m_boxSort = BoxSorting::None;
  }
  else if (str == "std") {
    m_boxSort = BoxSorting::Std;
  }
  else if (str == "shuffle") {
    m_boxSort = BoxSorting::Shuffle;
  }
  else if (str == "morton") {
    m_boxSort = BoxSorting::Morton;
  }
  else {
    MayDay::Abort(
      "ItoPlasmaGodunovStepper::parseOptions - unknown box sorting method requested for argument 'BoxSorting'");
  }

  // Parse filterse
  this->parseFilters();

  // Setup runtime storage (requirements change with algorithm)
  this->setupRuntimeStorage();

  //
  m_ito->parseRuntimeOptions();
  m_fieldSolver->parseRuntimeOptions();
  m_rte->parseRuntimeOptions();
}

void
ItoPlasmaGodunovStepper::allocateInternals()
{
  CH_TIME("ItoPlasmaGodunovStepper::allocateInternals");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocateInternals" << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_rtSpecies  = m_physics->getNumRtSpecies();

  m_amr->allocate(m_fluid_scratch1, m_fluidRealm, m_phase, 1);
  m_amr->allocate(m_fluid_scratchD, m_fluidRealm, m_phase, SpaceDim);

  m_amr->allocate(m_particle_scratch1, m_particleRealm, m_phase, 1);
  m_amr->allocate(m_particle_scratchD, m_particleRealm, m_phase, SpaceDim);
  m_amr->allocate(m_particle_E, m_particleRealm, m_phase, SpaceDim);

  m_amr->allocate(m_J, m_fluidRealm, m_phase, SpaceDim);
  m_amr->allocate(m_scratch1, m_fluidRealm, m_phase, 1);
  m_amr->allocate(m_scratch2, m_fluidRealm, m_phase, 1);
  m_amr->allocate(m_conduct_cell, m_fluidRealm, m_phase, 1);
  m_amr->allocate(m_conduct_face, m_fluidRealm, m_phase, 1);
  m_amr->allocate(m_conduct_eb, m_fluidRealm, m_phase, 1);
  m_amr->allocate(m_fluid_E, m_fluidRealm, m_phase, SpaceDim);

  // Allocate for energy sources
  m_energy_sources.resize(num_ItoSpecies);
  for (int i = 0; i < m_energy_sources.size(); i++) {
    m_amr->allocate(m_energy_sources[i], m_particleRealm, m_phase, 1);
  }

  // Allocate fluid scratch storage
  m_fscratch1.resize(num_ItoSpecies);
  m_fscratch2.resize(num_ItoSpecies);
  for (int i = 0; i < num_ItoSpecies; i++) {
    m_amr->allocate(m_fscratch1[i], m_fluidRealm, m_phase, 1);
    m_amr->allocate(m_fscratch2[i], m_fluidRealm, m_phase, 1);
  }

  // Allocate for PPC and YPC on both Realm. Also do EdotJ.
  m_amr->allocate(m_particle_ppc, m_particleRealm, m_phase, num_ItoSpecies);
  m_amr->allocate(m_particle_old, m_particleRealm, m_phase, num_ItoSpecies);
  m_amr->allocate(m_particle_eps, m_particleRealm, m_phase, num_ItoSpecies);
  m_amr->allocate(m_particle_ypc, m_particleRealm, m_phase, num_rtSpecies);

  m_amr->allocate(m_fluid_ppc, m_fluidRealm, m_phase, num_ItoSpecies);
  m_amr->allocate(m_fluid_eps, m_fluidRealm, m_phase, num_ItoSpecies);
  m_amr->allocate(m_fluid_ypc, m_fluidRealm, m_phase, num_rtSpecies);

  m_amr->allocate(m_EdotJ, m_fluidRealm, m_phase, num_ItoSpecies);
}

Real
ItoPlasmaGodunovStepper::advance(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::advance");
  if (m_verbosity > 5) {
    pout() << m_name + "::advance" << endl;
  }

  Real particle_time = 0.0;
  Real relax_time    = 0.0;
  Real Photon_time   = 0.0;
  Real sort_time     = 0.0;
  Real super_time    = 0.0;
  Real reaction_time = 0.0;
  Real clear_time    = 0.0;
  Real deposit_time  = 0.0;
  Real velo_time     = 0.0;
  Real diff_time     = 0.0;
  Real total_time    = 0.0;

  // Particle algorithms
  MPI_Barrier(Chombo_MPI::comm);
  total_time = -Timer::wallClock();
  particle_time -= Timer::wallClock();
  switch (m_algorithm) {
  case which_algorithm::euler_maruyama:
    this->advanceParticlesEulerMaruyama(a_dt);
    break;
  case which_algorithm::trapezoidal:
    this->advanceParticlesTrapezoidal(a_dt);
    break;
  default:
    MayDay::Abort("ItoPlasmaGodunovStepper::advance - logic bust");
  }
  particle_time += Timer::wallClock();

  // Compute current and relaxation time.
  MPI_Barrier(Chombo_MPI::comm);
  relax_time = -Timer::wallClock();
  this->computeJ(m_J, a_dt);
  m_dt_relax = this->computeRelaxationTime(); // This is for the restricting the next step.
  relax_time += Timer::wallClock();

  // Move Photons
  MPI_Barrier(Chombo_MPI::comm);
  Photon_time = -Timer::wallClock();
  this->advancePhotons(a_dt);
  Photon_time += Timer::wallClock();

  // If we are using the LEA, we must compute the Ohmic heating term. This must be done
  // BEFORE sorting the particles per cell.
  if (m_physics->getCoupling() == ItoPlasmaPhysics::coupling::LEA) {
    this->computeEdotJSource(a_dt);
  }

  // Sort the particles and Photons per cell so we can call reaction algorithms
  MPI_Barrier(Chombo_MPI::comm);
  sort_time = -Timer::wallClock();
  m_ito->sortParticlesByCell(ItoSolver::WhichContainer::Bulk);
  this->sortPhotonsByCell(McPhoto::WhichContainer::Bulk);
  this->sortPhotonsByCell(McPhoto::WhichContainer::Source);
  sort_time += Timer::wallClock();

  // Chemistry kernel.
  MPI_Barrier(Chombo_MPI::comm);
  reaction_time = -Timer::wallClock();
  this->advanceReactionNetwork(a_dt);
  reaction_time += Timer::wallClock();

  // Make superparticles.
  MPI_Barrier(Chombo_MPI::comm);
  super_time = -Timer::wallClock();
  if ((m_timeStep + 1) % m_merge_interval == 0 && m_merge_interval > 0) {
    m_ito->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_ppc);
  }
  super_time += Timer::wallClock();

  // Sort particles per patch.
  MPI_Barrier(Chombo_MPI::comm);
  sort_time -= Timer::wallClock();
  m_ito->sortParticlesByPatch(ItoSolver::WhichContainer::Bulk);
  this->sortPhotonsByPatch(McPhoto::WhichContainer::Bulk);
  this->sortPhotonsByPatch(McPhoto::WhichContainer::Source);
  sort_time += Timer::wallClock();

  // Clear other data holders for now. BC comes later...
  MPI_Barrier(Chombo_MPI::comm);
  clear_time = -Timer::wallClock();
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->clear(ItoSolver::WhichContainer::EB);
    solver_it()->clear(ItoSolver::WhichContainer::Domain);
  }
  clear_time += Timer::wallClock();

  //
  MPI_Barrier(Chombo_MPI::comm);
  deposit_time -= Timer::wallClock();
  m_ito->depositParticles();
  deposit_time += Timer::wallClock();

  // Prepare next step
  MPI_Barrier(Chombo_MPI::comm);
  velo_time -= Timer::wallClock();
  this->computeItoVelocities();
  velo_time += Timer::wallClock();

  MPI_Barrier(Chombo_MPI::comm);
  diff_time -= Timer::wallClock();
  this->computeItoDiffusion();
  diff_time += Timer::wallClock();

  total_time += Timer::wallClock();

  if (m_profile) {

    // Convert to %
    particle_time *= 100. / total_time;
    relax_time *= 100. / total_time;
    Photon_time *= 100. / total_time;
    sort_time *= 100. / total_time;
    super_time *= 100. / total_time;
    reaction_time *= 100. / total_time;
    clear_time *= 100. / total_time;
    deposit_time *= 100. / total_time;
    velo_time *= 100. / total_time;
    diff_time *= 100. / total_time;

    // Total percentage/imbalance
    Real imbalance = 0.0;
    imbalance += particle_time;
    imbalance += relax_time;
    imbalance += Photon_time;
    imbalance += sort_time;
    imbalance += super_time;
    imbalance += reaction_time;
    imbalance += clear_time;
    imbalance += deposit_time;
    imbalance += velo_time;
    imbalance += diff_time;
    imbalance = 100. - imbalance;

    pout() << "\n";
    pout() << "ItoPlasmaGodunovStepper::advance breakdown:" << endl << "======================================" << endl;
    printTimerHead();
    printTimerDiagnostics(particle_time, "Transport & Poisson (%)");
    printTimerDiagnostics(relax_time, "Relax time (%)");
    printTimerDiagnostics(Photon_time, "Photons (%)");
    printTimerDiagnostics(sort_time, "Sort (%)");
    printTimerDiagnostics(super_time, "Superparticles (%)");
    printTimerDiagnostics(reaction_time, "Reaction network (%)");
    printTimerDiagnostics(clear_time, "EB removal (%)");
    printTimerDiagnostics(deposit_time, "Deposition (%)");
    printTimerDiagnostics(velo_time, "Velo comp (%)");
    printTimerDiagnostics(diff_time, "Diff comp (%)");
    printTimerDiagnostics(imbalance, "Imbalance (%)");
    printTimerDiagnostics(total_time, "Total time (s)");
    printTimerTail();
    pout() << "\n";
  }

  return a_dt;
}

Real
ItoPlasmaGodunovStepper::computeDt()
{
  CH_TIME("ItoPlasmaGodunovStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::computeDt" << endl;
  }

  Real a_dt = std::numeric_limits<Real>::max();

  if (m_whichDt == which_dt::advection) {
    a_dt = m_ito->computeAdvectiveDt();
  }
  else if (m_whichDt == which_dt::diffusion) {
    a_dt = m_ito->computeDiffusiveDt();
  }
  else if (m_whichDt == which_dt::AdvectionDiffusion) {
    a_dt = m_ito->computeDt();
  }

  a_dt                = a_dt * m_max_cells_hop;
  TimeCode a_timeCode = TimeCode::Advection;

  // Physics-based restriction
  const Real physicsDt = this->computePhysicsDt();
  if (physicsDt < a_dt) {
    a_dt       = physicsDt;
    a_timeCode = TimeCode::Physics;
  }

  if (a_dt < m_min_dt) {
    a_dt       = m_min_dt;
    a_timeCode = TimeCode::Hardcap;
  }

  if (a_dt > m_max_dt) {
    a_dt       = m_max_dt;
    a_timeCode = TimeCode::Hardcap;
  }

  m_timeCode = a_timeCode;

#if 0 // Debug code
  const Real dtCFL = m_ito->computeDt();
  m_avg_cfl += a_dt/dtCFL;
  if(procID() == 0) std::cout << "dt = " << a_dt
			      << "\t relax dt = " << m_dt_relax
			      << "\t factor = " << a_dt/m_dt_relax
			      << "\t CFL = " << a_dt/dtCFL
			      << "\t avgCFL = " << m_avg_cfl/(1+m_timeStep)
			      << std::endl;
#endif

  return a_dt;
}

void
ItoPlasmaGodunovStepper::preRegrid(const int a_lmin, const int a_oldFinestLevel)
{
  CH_TIME("ItoPlasmaGodunovStepper::preRegrid");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::preRegrid" << endl;
  }

  ItoPlasmaStepper::preRegrid(a_lmin, a_oldFinestLevel);

  // Copy conductivity to scratch storage
  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();
  m_amr->allocate(m_cache, m_fluidRealm, m_phase, ncomp);
  for (int lvl = 0; lvl <= a_oldFinestLevel; lvl++) {
    m_conduct_cell[lvl]->localCopyTo(*m_cache[lvl]);
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const int idx = solver_it.index();

    m_conductivity_particles[idx]->preRegrid(a_lmin);
    m_rho_dagger_particles[idx]->preRegrid(a_lmin);
  }
}

void
ItoPlasmaGodunovStepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("ItoPlasmaGodunovStepper::regrid");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::regrid" << endl;
  }

  Real ito_time      = 0.0;
  Real poisson_time  = 0.0;
  Real rte_time      = 0.0;
  Real sigma_time    = 0.0;
  Real internal_time = 0.0;

  Real gdnv_time    = 0.0;
  Real setup_time   = 0.0;
  Real solve_time   = 0.0;
  Real super_time   = 0.0;
  Real cleanup_time = 0.0;
  Real total_time   = 0.0;

  // Regrid solvers
  total_time -= Timer::wallClock();
  ito_time -= Timer::wallClock();
  m_ito->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  ito_time += Timer::wallClock();

  MPI_Barrier(Chombo_MPI::comm);
  poisson_time -= Timer::wallClock();
  m_fieldSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  poisson_time += Timer::wallClock();

  MPI_Barrier(Chombo_MPI::comm);
  rte_time -= Timer::wallClock();
  m_rte->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  rte_time += Timer::wallClock();

  MPI_Barrier(Chombo_MPI::comm);
  sigma_time -= Timer::wallClock();
  m_sigma->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  sigma_time += Timer::wallClock();

  // Allocate internal memory for ItoPlasmaGodunovStepper now....
  MPI_Barrier(Chombo_MPI::comm);
  internal_time -= Timer::wallClock();
  this->allocateInternals();
  internal_time += Timer::wallClock();

  // We need to remap/regrid the stored particles as well.
  MPI_Barrier(Chombo_MPI::comm);
  gdnv_time -= Timer::wallClock();
  const Vector<DisjointBoxLayout>& grids   = m_amr->getGrids(m_particleRealm);
  const Vector<ProblemDomain>&     domains = m_amr->getDomains();
  const Vector<Real>&              dx      = m_amr->getDx();
  const Vector<int>&               ref_rat = m_amr->getRefinementRatios();

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const int idx = solver_it.index();
    //    m_amr->regrid(*m_rho_dagger_particles[idx], a_lmin, a_newFinestLevel);
    //    m_amr->regrid(*m_conductivity_particles[idx], a_lmin, a_newFinestLevel);
  }
  gdnv_time += Timer::wallClock();

  // Recompute the conductivity and space charge densities.
  MPI_Barrier(Chombo_MPI::comm);
  setup_time -= Timer::wallClock();
  this->computeRegridConductivity();
  this->computeRegridRho();
  this->setupSemiImplicitPoisson(m_prevDt);
  setup_time += Timer::wallClock();

  MPI_Barrier(Chombo_MPI::comm);
  solve_time -= Timer::wallClock();
  const bool converged = this->solvePoisson();
  if (!converged) {
    MayDay::Abort("ItoPlasmaGodunovStepper::regrid - Poisson solve did not converge after regrid!!!");
  }
  solve_time += Timer::wallClock();

  // Regrid superparticles.
  MPI_Barrier(Chombo_MPI::comm);
  super_time -= Timer::wallClock();
  if (m_regridSuperparticles) {
    m_ito->sortParticlesByCell(ItoSolver::WhichContainer::Bulk);
    m_ito->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_ppc);
    m_ito->sortParticlesByPatch(ItoSolver::WhichContainer::Bulk);
  }
  super_time += Timer::wallClock();

  // Now let the ito solver deposit its actual particles... In the above it deposit m_rho_dagger_particles.
  MPI_Barrier(Chombo_MPI::comm);
  cleanup_time -= Timer::wallClock();
  m_ito->depositParticles();

  // Recompute new velocities and diffusion coefficients
  this->computeItoVelocities();
  this->computeItoDiffusion();
  cleanup_time += Timer::wallClock();

  MPI_Barrier(Chombo_MPI::comm);
  total_time += Timer::wallClock();

  if (m_profile) {

    // Convert to %
    ito_time *= 100. / total_time;
    poisson_time *= 100. / total_time;
    rte_time *= 100. / total_time;
    sigma_time *= 100. / total_time;
    internal_time *= 100. / total_time;
    gdnv_time *= 100. / total_time;
    setup_time *= 100. / total_time;
    solve_time *= 100. / total_time;
    super_time *= 100. / total_time;
    cleanup_time *= 100. / total_time;

    // Total percentage/imbalance
    Real imbalance = 0.0;
    imbalance += ito_time;
    imbalance += poisson_time;
    imbalance += rte_time;
    imbalance += sigma_time;
    imbalance += internal_time;
    imbalance += solve_time;
    imbalance += super_time;
    imbalance += cleanup_time;
    imbalance = 100. - imbalance;

    pout() << "\n";
    pout() << "ItoPlasmaGodunovStepper::regrid breakdown:" << endl << "======================================" << endl;
    printTimerHead();
    printTimerDiagnostics(ito_time, "Ito regrid (%)");
    printTimerDiagnostics(poisson_time, "Poisson regrid (%)");
    printTimerDiagnostics(rte_time, "RTE regrid (%)");
    printTimerDiagnostics(sigma_time, "Sigma regrid (%)");
    printTimerDiagnostics(internal_time, "Internal regrid (%)");
    printTimerDiagnostics(gdnv_time, "Gdnv particles (%)");
    printTimerDiagnostics(setup_time, "Poisson setup (%)");
    printTimerDiagnostics(solve_time, "Poisson solve (%)");
    printTimerDiagnostics(super_time, "Super time (%)");
    printTimerDiagnostics(cleanup_time, "Cleanup (%)");
    printTimerDiagnostics(total_time, "Total time (s)");
    printTimerTail();
    pout() << "\n";
  }
}

void
ItoPlasmaGodunovStepper::setupRuntimeStorage()
{
  CH_TIME("ItoPlasmaGodunovStepper::setupRuntimeStorage");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupRuntimeStorage" << endl;
  }
#if 0
  switch (m_algorithm) {
  case which_algorithm::euler_maruyama:
    ItoParticle::setNumRuntimeVectors(1);
    break;
  case which_algorithm::trapezoidal:
    ItoParticle::setNumRuntimeVectors(2); // For V^k and the diffusion hop.
    break;
  default:
    MayDay::Abort("ItoPlasmaGodunovStepper::setupRuntimeStorage - logic bust");
  }
#else
  MayDay::Error(
    "ItoPlasmaGodunovStepper::setupRuntimeStorage -- need to figure out how to add more fields to ItoParticle");
#endif
}

void
ItoPlasmaGodunovStepper::setOldPositions()
{
  CH_TIME("ItoPlasmaGodunovStepper::setOldPositions()");
  if (m_verbosity > 5) {
    pout() << m_name + "::setOldPositions()" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
      ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

        List<ItoParticle>& particleList = particles[dit()].listItems();

        for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
          ItoParticle& p  = particleList[lit];
          p.oldPosition() = p.position();
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::remapGodunovParticles(Vector<ParticleContainer<PointParticle>*>& a_particles,
                                               const WhichParticles                       a_WhichParticles)
{
  CH_TIME("ItoPlasmaGodunovStepper::remapGodunovParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::remapGodunovParticles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_WhichParticles) {
    case WhichParticles::All:
      a_particles[idx]->remap();
      break;
    case WhichParticles::AllMobile:
      if (mobile)
        a_particles[idx]->remap();
      break;
    case WhichParticles::AllDiffusive:
      if (diffusive)
        a_particles[idx]->remap();
      break;
    case WhichParticles::ChargedMobile:
      if (charged && mobile)
        a_particles[idx]->remap();
      break;
    case WhichParticles::ChargedDiffusive:
      if (charged && diffusive)
        a_particles[idx]->remap();
      break;
    case WhichParticles::AllMobileOrDiffusive:
      if (mobile || diffusive)
        a_particles[idx]->remap();
      break;
    case WhichParticles::ChargedAndMobileOrDiffusive:
      if (charged && (mobile || diffusive))
        a_particles[idx]->remap();
      break;
    case WhichParticles::Stationary:
      if (!mobile && !diffusive)
        a_particles[idx]->remap();
      break;
    default:
      MayDay::Abort("ItoPlasmaGodunovStepper::remapGodunovParticles - logic bust");
    }
  }
}

void
ItoPlasmaGodunovStepper::deposit_PointParticles(const Vector<ParticleContainer<PointParticle>*>& a_particles,
                                                const WhichParticles                             a_WhichParticles)
{
  CH_TIME("ItoPlasmaGodunovStepper::deposit_PointParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::deposit_PointParticles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_WhichParticles) {
    case WhichParticles::All:
      solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case WhichParticles::AllMobile:
      if (mobile)
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case WhichParticles::AllDiffusive:
      if (diffusive)
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case WhichParticles::ChargedMobile:
      if (charged && mobile)
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case WhichParticles::ChargedDiffusive:
      if (charged && diffusive)
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case WhichParticles::AllMobileOrDiffusive:
      if (mobile || diffusive)
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case WhichParticles::ChargedAndMobileOrDiffusive:
      if (charged && (mobile || diffusive))
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case WhichParticles::Stationary:
      if (!mobile && !diffusive)
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    default:
      MayDay::Abort("ItoPlasmaGodunovStepper::deposit_PointParticles - logic bust");
    }
  }
}

void
ItoPlasmaGodunovStepper::clearGodunovParticles(const Vector<ParticleContainer<PointParticle>*>& a_particles,
                                               const WhichParticles                             a_WhichParticles)
{
  CH_TIME("ItoPlasmaGodunovStepper::clearGodunovParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::deposit_clearParticles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_WhichParticles) {
    case WhichParticles::All:
      a_particles[idx]->clearParticles();
      break;
    case WhichParticles::AllMobile:
      if (mobile)
        a_particles[idx]->clearParticles();
      break;
    case WhichParticles::AllDiffusive:
      if (diffusive)
        a_particles[idx]->clearParticles();
      break;
    case WhichParticles::ChargedMobile:
      if (charged && mobile)
        a_particles[idx]->clearParticles();
      break;
    case WhichParticles::ChargedDiffusive:
      if (charged && diffusive)
        a_particles[idx]->clearParticles();
      break;
    case WhichParticles::AllMobileOrDiffusive:
      if (mobile || diffusive)
        a_particles[idx]->clearParticles();
      break;
    case WhichParticles::ChargedAndMobileOrDiffusive:
      if (charged && (mobile || diffusive))
        a_particles[idx]->clearParticles();
      break;
    case WhichParticles::Stationary:
      if (!mobile && !diffusive)
        a_particles[idx]->clearParticles();
      break;
    default:
      MayDay::Abort("ItoPlasmaGodunovStepper::clearGodunovParticles - logic bust");
    }
  }
}

void
ItoPlasmaGodunovStepper::computeAllConductivities(const Vector<ParticleContainer<PointParticle>*>& a_particles)
{
  CH_TIME("ItoPlasmaGodunovStepper::computeAllConductivities");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAllConductivities" << endl;
  }

  this->compute_cell_conductivity(m_conduct_cell, a_particles);

  // Now do the faces
  this->compute_face_conductivity();
}

void
ItoPlasmaGodunovStepper::compute_cell_conductivity(EBAMRCellData&                                   a_conductivity,
                                                   const Vector<ParticleContainer<PointParticle>*>& a_particles)
{
  CH_TIME("ItoPlasmaGodunovStepper::compute_cell_conductivity(conductivity, PointParticle");
  if (m_verbosity > 5) {
    pout() << m_name + "::compute_cell_conductivity(conductivity, PointParticle)" << endl;
  }

  DataOps::setValue(a_conductivity, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    if (q != 0 && solver->isMobile()) {
      DataOps::setValue(m_particle_scratch1, 0.0);
#if 1 // Original code
      solver->depositParticles<PointParticle,
                               &PointParticle::weight>(m_particle_scratch1,
                                                       *a_particles[idx]); // The particles should have "masses" = m*mu
#else
      const EBAMRCellData& mu  = solver->getMobilityFunction();
      const EBAMRCellData& phi = solver->getPhi();
      DataOps::copy(m_particle_scratch1, mu);
      DataOps::multiply(m_particle_scratch1, phi);
#endif

      // Copy to fluid Realm and add to fluid stuff
      m_fluid_scratch1.copy(m_particle_scratch1);
      DataOps::incr(a_conductivity, m_fluid_scratch1, Abs(q));
    }
  }

  // Test code
  if (m_filter_cond) {
    for (const auto& f : m_filters) {
      const Real alpha   = std::get<0>(f);
      const int  stride  = std::get<1>(f);
      const int  num_app = std::get<2>(f);

      for (int iapp = 0; iapp < num_app; iapp++) {
        DataOps::setValue(m_fluid_scratch1, 0.0);
        m_fluid_scratch1.copy(a_conductivity);
        DataOps::setCoveredValue(m_fluid_scratch1, 0.0, 0);
        //        DataOps::filterSmooth(a_conductivity, m_fluid_scratch1, stride, alpha);

        MayDay::Error("ItoPlasmaGodunovStepper::compute_cell_conductivity - don't use filterSmooth");

        m_amr->conservativeAverage(a_conductivity, m_fluidRealm, m_phase);
        m_amr->interpGhost(a_conductivity, m_fluidRealm, m_phase);
      }
    }
  }

  DataOps::scale(a_conductivity, Units::Qe);

  m_amr->conservativeAverage(a_conductivity, m_fluidRealm, m_phase);
  m_amr->interpGhostPwl(a_conductivity, m_fluidRealm, m_phase);

  // See if this helps....
  m_amr->interpToCentroids(a_conductivity, m_fluidRealm, m_phase);
}

void
ItoPlasmaGodunovStepper::compute_face_conductivity()
{
  CH_TIME("ItoPlasmaGodunovStepper::compute_face_conductivity");
  if (m_verbosity > 5) {
    pout() << m_name + "::compute_face_conductivity" << endl;
  }

  DataOps::setValue(m_conduct_face, 0.0);
  DataOps::setValue(m_conduct_eb, 0.0);

  // This code does averaging from cell to face.
  DataOps::averageCellToFace(m_conduct_face, m_conduct_cell, m_amr->getDomains());

  // This code extrapolates the conductivity to the EB. This should actually be the EB centroid but since the stencils
  // for EB extrapolation can be a bit nasty (e.g. Negative weights), we do the centroid instead and take that as an approximation.
#if 0
  const IrregAmrStencil<CentroidInterpolationStencil>& ebsten = m_amr->getCentroidInterpolationStencils(m_fluidRealm, m_phase);
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    ebsten.apply(m_conduct_eb, m_conduct_cell, lvl);
  }
#else
  DataOps::incr(m_conduct_eb, m_conduct_cell, 1.0);
#endif
}

void
ItoPlasmaGodunovStepper::setupSemiImplicitPoisson(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::setupSemiImplicitPoisson");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupSemiImplicitPoisson" << endl;
  }

  // Set coefficients as usual
  m_fieldSolver->setPermittivities();

  // Get bco and increment with mobilities
  MFAMRFluxData& bco     = m_fieldSolver->getPermittivityFace();
  MFAMRIVData&   bco_irr = m_fieldSolver->getPermittivityEB();

  EBAMRFluxData bco_gas;
  EBAMRIVData   bco_irr_gas;

  m_amr->allocatePointer(bco_gas);
  m_amr->allocatePointer(bco_irr_gas);

  m_amr->alias(bco_gas, phase::gas, bco);
  m_amr->alias(bco_irr_gas, phase::gas, bco_irr);

  DataOps::scale(m_conduct_face, a_dt / Units::eps0);
  DataOps::scale(m_conduct_eb, a_dt / Units::eps0);

  DataOps::multiply(m_conduct_face, bco_gas);
  DataOps::multiply(m_conduct_eb, bco_irr_gas);

  DataOps::incr(bco_gas, m_conduct_face, 1.0);
  DataOps::incr(bco_irr_gas, m_conduct_eb, 1.0);

  m_amr->conservativeAverage(bco_gas, m_fluidRealm, phase::gas);
  m_amr->conservativeAverage(bco_irr_gas, m_fluidRealm, phase::gas);

  // Set up the solver
  m_fieldSolver->setupSolver();
}

void
ItoPlasmaGodunovStepper::setupStandardPoisson()
{
  CH_TIME("ItoPlasmaGodunovStepper::setupStandardPoisson");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupStandardPoisson" << endl;
  }

  // Set coefficients as usual
  m_fieldSolver->setPermittivities();
  m_fieldSolver->setupSolver();
}

void
ItoPlasmaGodunovStepper::copyConductivityParticles(Vector<ParticleContainer<PointParticle>*>& a_conductivity_particles)
{
  CH_TIME("ItoPlasmaGodunovStepper::copyConductivityParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::copyConductivityParticles" << endl;
  }

  this->clearGodunovParticles(a_conductivity_particles, WhichParticles::All);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const RefCountedPtr<ItoSolver>&  solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
        const List<ItoParticle>& ito_parts =
          solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl][dit()].listItems();
        List<PointParticle>& gdnv_parts = (*a_conductivity_particles[idx])[lvl][dit()].listItems();

        if (q != 0 && solver->isMobile()) {
          for (ListIterator<ItoParticle> lit(ito_parts); lit.ok(); ++lit) {
            const ItoParticle& p        = lit();
            const RealVect&    pos      = p.position();
            const Real&        weight   = p.weight();
            const Real&        mobility = p.mobility();

            gdnv_parts.add(PointParticle(pos, weight * mobility));
          }
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::copyRhoDaggerParticles(Vector<ParticleContainer<PointParticle>*>& a_rho_dagger_particles)
{
  CH_TIME("ItoPlasmaGodunovStepper::copyRhoDaggerParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::copyRhoDaggerParticles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const RefCountedPtr<ItoSolver>&  solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
        const List<ItoParticle>& ito_parts =
          solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl][dit()].listItems();
        List<PointParticle>& gdnv_parts = (*a_rho_dagger_particles[idx])[lvl][dit()].listItems();

        gdnv_parts.clear();

        if (q != 0) {
          for (ListIterator<ItoParticle> lit(ito_parts); lit.ok(); ++lit) {
            const ItoParticle& p      = lit();
            const RealVect&    pos    = p.position();
            const Real&        weight = p.weight();

            gdnv_parts.add(PointParticle(pos, weight));
          }
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::computeRegridConductivity()
{
  CH_TIME("ItoPlasmaGodunovStepper::computeRegridConductivity");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeRegridConductivity" << endl;
  }

  this->computeAllConductivities(m_conductivity_particles);
}

void
ItoPlasmaGodunovStepper::computeRegridRho()
{
  CH_TIME("ItoPlasmaGodunovStepper::computeRegridRho");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeRegridRho" << endl;
  }

  this->deposit_PointParticles(m_rho_dagger_particles, WhichParticles::All);
}

void
ItoPlasmaGodunovStepper::advanceParticlesEulerMaruyama(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::advanceParticlesEulerMaruyama");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceParticlesEulerMaruyama" << endl;
  }

  Real posTime         = 0.0;
  Real diffuseTime     = 0.0;
  Real remapGdnvTime   = 0.0;
  Real depositGdnvTime = 0.0;
  Real copyCondTime    = 0.0;
  Real condTime        = 0.0;
  Real setupTime       = 0.0;
  Real poissonTime     = 0.0;
  Real velocityTime    = 0.0;
  Real particleTime    = 0.0;
  Real remapTime       = 0.0;
  Real isectTime       = 0.0;
  Real depositTime     = 0.0;
  Real totalTime       = 0.0;

  totalTime -= Timer::wallClock();

  m_prevDt = a_dt; // Needed for regrids.

  // 1. Store X^k positions.
  MPI_Barrier(Chombo_MPI::comm);
  posTime -= Timer::wallClock();
  this->setOldPositions();
  posTime += Timer::wallClock();

  // 2. Diffuse the particles. This copies onto m_rho_dagger_particles and stores the hop on the full particles.
  MPI_Barrier(Chombo_MPI::comm);
  diffuseTime -= Timer::wallClock();
  this->diffuseParticlesEulerMaruyama(m_rho_dagger_particles, a_dt);
  diffuseTime += Timer::wallClock();

  MPI_Barrier(Chombo_MPI::comm);
  remapGdnvTime -= Timer::wallClock();
  this->remapGodunovParticles(m_rho_dagger_particles, WhichParticles::AllDiffusive);
  remapGdnvTime += Timer::wallClock();

  // 3. Solve the semi-implicit Poisson equation. Also, copy the particles used for computing the conductivity to scratch.
  MPI_Barrier(Chombo_MPI::comm);
  copyCondTime -= Timer::wallClock();
  this->copyConductivityParticles(m_conductivity_particles); // Sets particle "weights" = w*mu
  copyCondTime += Timer::wallClock();

  // Compute conductivity on mesh
  MPI_Barrier(Chombo_MPI::comm);
  condTime -= Timer::wallClock();
  this->computeAllConductivities(m_conductivity_particles); // Deposits q_e*Z*w*mu on the mesh
  condTime += Timer::wallClock();

  // Setup Poisson solver
  MPI_Barrier(Chombo_MPI::comm);
  setupTime -= Timer::wallClock();
  this->setupSemiImplicitPoisson(a_dt); // Multigrid setup
  setupTime += Timer::wallClock();

  // Compute space charge density
  MPI_Barrier(Chombo_MPI::comm);
  depositGdnvTime -= Timer::wallClock();
  this->deposit_PointParticles(m_rho_dagger_particles,
                               WhichParticles::
                                 AllDiffusive); // Diffusive should be enough because state is not changed for others.
  depositGdnvTime += Timer::wallClock();

  MPI_Barrier(Chombo_MPI::comm);
  poissonTime -= Timer::wallClock();
  this->solvePoisson(); // Solve the stinking equation.
  poissonTime += Timer::wallClock();

  // 4. Recompute velocities with the new electric field, then do the actual semi-implicit Euler-Maruyama update.
  MPI_Barrier(Chombo_MPI::comm);
  velocityTime -= Timer::wallClock();
#if 1 // This is what the algorithm says.
  this->setItoVelocityFunctions();
  m_ito->interpolateVelocities();
#else // Have to use this for LEA - need to debug.
  this->computeItoVelocities();
#endif
  velocityTime += Timer::wallClock();

  MPI_Barrier(Chombo_MPI::comm);
  particleTime -= Timer::wallClock();
  this->stepEulerMaruyama(a_dt);
  particleTime += Timer::wallClock();

  MPI_Barrier(Chombo_MPI::comm);
  remapTime -= Timer::wallClock();
  this->remapParticles(WhichParticles::AllMobileOrDiffusive);
  remapTime += Timer::wallClock();

  // 5. Do intersection test and remove EB particles. These particles are NOT allowed to react later.
  MPI_Barrier(Chombo_MPI::comm);
  isectTime -= Timer::wallClock();
  const bool delete_eb_particles = true;
  this->intersectParticles(WhichParticles::AllMobileOrDiffusive,
                           EbRepresentation::ImplicitFunction,
                           delete_eb_particles);
  this->removeCoveredParticles(WhichParticles::AllMobileOrDiffusive,
                               EbRepresentation::ImplicitFunction,
                               m_eb_tolerance);
  isectTime += Timer::wallClock();

  // 6. Deposit particles. This shouldn't be necessary unless we want to compute (E,J)
  MPI_Barrier(Chombo_MPI::comm);
  depositTime -= Timer::wallClock();
  this->depositParticles(WhichParticles::AllMobileOrDiffusive);
  depositTime += Timer::wallClock();

  totalTime += Timer::wallClock();

  if (m_profile) {

    posTime *= 100. / totalTime;
    diffuseTime *= 100. / totalTime;
    remapGdnvTime *= 100. / totalTime;
    depositGdnvTime *= 100. / totalTime;
    copyCondTime *= 100. / totalTime;
    condTime *= 100. / totalTime;
    setupTime *= 100. / totalTime;
    poissonTime *= 100. / totalTime;
    velocityTime *= 100. / totalTime;
    particleTime *= 100. / totalTime;
    remapTime *= 100. / totalTime;
    isectTime *= 100. / totalTime;
    depositTime *= 100. / totalTime;

    Real imbalance = 0.0;
    imbalance += posTime;
    imbalance += diffuseTime;
    imbalance += remapGdnvTime;
    imbalance += copyCondTime;
    imbalance += condTime;
    imbalance += setupTime;
    imbalance += poissonTime;
    imbalance += velocityTime;
    imbalance += particleTime;
    imbalance += remapTime;
    imbalance += isectTime;
    imbalance += depositTime;
    imbalance = 100. - imbalance;

    pout() << "\n";
    pout() << "ItoPlasmaGodunovStepper::euler_maruyama breakdown:" << endl
           << "======================================" << endl;
    printTimerHead();
    printTimerDiagnostics(posTime, "Old position (%)");
    printTimerDiagnostics(diffuseTime, "Diffusion (%)");
    printTimerDiagnostics(remapGdnvTime, "Remap gdnv (%)");
    printTimerDiagnostics(depositGdnvTime, "Deposit gdnv (%)");
    printTimerDiagnostics(copyCondTime, "Copy conductivity (%)");
    printTimerDiagnostics(condTime, "Conductity comp (%)");
    printTimerDiagnostics(setupTime, "Poisson setup (%)");
    printTimerDiagnostics(poissonTime, "Poisson solve (%)");
    printTimerDiagnostics(velocityTime, "Velo comp (%)");
    printTimerDiagnostics(particleTime, "Advect particles (%)");
    printTimerDiagnostics(remapTime, "Remap particles (%)");
    printTimerDiagnostics(isectTime, "Particle intersection (%)");
    printTimerDiagnostics(depositTime, "Deposition (%)");
    printTimerDiagnostics(imbalance, "Imbalance (%)");
    printTimerDiagnostics(totalTime, "Total time (s)");
    printTimerTail();
    pout() << "\n";
  }
}

void
ItoPlasmaGodunovStepper::diffuseParticlesEulerMaruyama(Vector<ParticleContainer<PointParticle>*>& a_rho_dagger,
                                                       const Real                                 a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::diffuseParticlesEulerMaruyama");
  if (m_verbosity > 5) {
    pout() << m_name + "::diffuseParticlesEulerMaruyama" << endl;
  }

  this->clearGodunovParticles(a_rho_dagger, WhichParticles::All);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;    // Multiplication factor for mobility
    const Real g = diffusive ? 1.0 : 0.0; // Multiplication factor for diffusion

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
      ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

        List<ItoParticle>&   ItoParticles = particles[dit()].listItems();
        List<PointParticle>& gdnv_parts   = (*a_rho_dagger[idx])[lvl][dit()].listItems();

        if (diffusive) {
          for (ListIterator<ItoParticle> lit(ItoParticles); lit.ok(); ++lit) {
            ItoParticle&    p      = lit();
            const Real      factor = g * sqrt(2.0 * p.diffusion() * a_dt);
            const Real&     weight = p.weight();
            const RealVect& pos    = p.position();
#if 0
	    RealVect&       hop    = p.runtimeVector(0);
            hop                    = factor * solver->randomGaussian();


            // Add simpler particle
            gdnv_parts.add(PointParticle(pos + hop, weight));
#else
            MayDay::Error("ItoPlasmaGodunovStepper::diffuseParticlesEulerMaruayma -- runtime stuff");
#endif
          }
        }
        else { // Splitting up diffusion and non-diffusion because I dont want to generate random numbers where they're not required...
          for (ListIterator<ItoParticle> lit(ItoParticles); lit.ok(); ++lit) {
            ItoParticle&    p      = lit();
            const Real&     weight = p.weight();
            const RealVect& pos    = p.position();
#if 0
	    RealVect&       hop  = p.runtimeVector(0);
            hop                  = RealVect::Zero;

            // Add simpler particle
            gdnv_parts.add(PointParticle(pos, weight));
#else
            MayDay::Error("ItoPlasmaGodunovStepper::diffuseParticlesEulerMaruayma -- runtime stuff");
#endif
          }
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::stepEulerMaruyama(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::stepEulerMaruyama");
  if (m_verbosity > 5) {
    pout() << m_name + "::stepEulerMaruyama" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;
    const Real g = diffusive ? 1.0 : 0.0;

    if (mobile || diffusive) {

      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
        ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

        for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

          List<ItoParticle>& particleList = particles[dit()].listItems();

          for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
            ItoParticle& p = lit();

            // Add diffusion hop again. The position after the diffusion hop is oldPosition() and X^k is in position()
#if 0
	    const RealVect& hop = p.runtimeVector(0);
            p.position()        = p.oldPosition() + f * p.velocity() * a_dt + g * hop;
#else
            MayDay::Error("ItoPlasmaGodunovStepper::stopEulerMaruyama -- runtime stuff");
#endif
          }
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::advanceParticlesTrapezoidal(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::advanceParticlesTrapezoidal");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceParticlesTrapezoidal" << endl;
  }

  m_prevDt = 0.5 * a_dt; // Needed for regrids.

  this->setOldPositions();

  // ====== PREDICTOR BEGIN ======
  this->preTrapezoidalPredictor(m_rho_dagger_particles, a_dt);
  this->remapGodunovParticles(m_rho_dagger_particles,
                              WhichParticles::
                                AllDiffusive); // Particles that were copied but not moved are in the right box.
  this->deposit_PointParticles(m_rho_dagger_particles, WhichParticles::All); // All copies need to deposit.

  this->copyConductivityParticles(m_conductivity_particles);
  this->computeAllConductivities(m_conductivity_particles);
  this->setupSemiImplicitPoisson(a_dt);
  this->solvePoisson();

  this->setItoVelocityFunctions();
  m_ito->interpolateVelocities();
  this->trapezoidalPredictor(a_dt);
  this->remapParticles(WhichParticles::AllMobileOrDiffusive);
  // ====== PREDICTOR END ======

  // ====== CORRECTOR BEGIN =====
  this->preTrapezoidalCorrector(m_rho_dagger_particles,
                                a_dt); // Mobile or diffusive moves to X^dagger = X^k + 0.5*dt*V^k + hop
  this->remapGodunovParticles(m_rho_dagger_particles,
                              WhichParticles::
                                AllMobileOrDiffusive); // Only need to remap particles that were mobile or diffusive
  this->deposit_PointParticles(m_rho_dagger_particles,
                               WhichParticles::All); // Everything needs to deposit...

  this->copyConductivityParticles(m_conductivity_particles);
  this->computeAllConductivities(m_conductivity_particles);
  this->setupSemiImplicitPoisson(0.5 * a_dt);
  this->solvePoisson();

  this->setItoVelocityFunctions();
  m_ito->interpolateVelocities();
  this->trapezoidalCorrector(a_dt);
  this->remapParticles(WhichParticles::AllMobileOrDiffusive);
  // ====== CORRECTOR END =====

  // Do particle-boundary intersection.
  this->intersectParticles(WhichParticles::AllMobileOrDiffusive, EbRepresentation::ImplicitFunction, true);
  this->removeCoveredParticles(WhichParticles::AllMobileOrDiffusive,
                               EbRepresentation::ImplicitFunction,
                               m_eb_tolerance);

  // Finally, deposit particles.
  this->depositParticles(WhichParticles::AllMobileOrDiffusive);
}

void
ItoPlasmaGodunovStepper::preTrapezoidalPredictor(Vector<ParticleContainer<PointParticle>*>& a_rho_dagger,
                                                 const Real                                 a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::preTrapezoidalPredictor");
  if (m_verbosity > 5) {
    pout() << m_name + "::preTrapezoidalPredictor" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;    // Multiplication factor for mobility
    const Real g = diffusive ? 1.0 : 0.0; // Multiplication factor for diffusion

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
      ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

        List<ItoParticle>&   ItoParticles = particles[dit()].listItems();
        List<PointParticle>& gdnv_parts   = (*a_rho_dagger[idx])[lvl][dit()].listItems();

        gdnv_parts.clear();

        // Store the diffusion hop, and add the godunov particles
        for (ListIterator<ItoParticle> lit(ItoParticles); lit.ok(); ++lit) {
          ItoParticle&    p      = lit();
          const Real      factor = sqrt(2.0 * p.diffusion() * a_dt);
          const RealVect  hop    = factor * solver->randomGaussian();
          const RealVect& Xk     = p.oldPosition();
          const Real&     weight = p.weight();

          // Store the diffusion hop and the current velocity.
#if 0
	  p.runtimeVector(0) = g * hop;
	  p.runtimeVector(1) = f * p.velocity();

          // Add simpler particle
          gdnv_parts.add(PointParticle(Xk + g * hop, weight));
#else
          MayDay::Error("ItoPlasmaGodunovStepper::preTrapezoidalPredictor -- runtime stuff");
#endif
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::trapezoidalPredictor(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::trapezoidalPredictor");
  if (m_verbosity > 5) {
    pout() << m_name + "::trapezoidalPredictor" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;
    const Real g = diffusive ? 1.0 : 0.0;

    if (mobile || diffusive) {

      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
        ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

        for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

          List<ItoParticle>& particleList = particles[dit()].listItems();

          for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
            ItoParticle& p = lit();

            // Add diffusion hop again. The position after the diffusion hop is oldPosition() and X^k is in position()
#if 0
	    const RealVect& hop = p.runtimeVector(0);
	    const RealVect& Vk  = p.runtimeVector(1);

            p.position() = p.oldPosition() + f * p.velocity() * a_dt + g * hop;
#else
            MayDay::Error("ItoPlasmaGodunovStepper::trapezoidalPredictor -- runtime stuff");
#endif
          }
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::preTrapezoidalCorrector(Vector<ParticleContainer<PointParticle>*>& a_rho_dagger,
                                                 const Real                                 a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::preTrapezoidalCorrector");
  if (m_verbosity > 5) {
    pout() << m_name + "::preTrapezoidalCorrector" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;    // Multiplication factor for mobility
    const Real g = diffusive ? 1.0 : 0.0; // Multiplication factor for diffusion

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
      ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

        List<ItoParticle>&   ItoParticles = particles[dit()].listItems();
        List<PointParticle>& gdnv_parts   = (*a_rho_dagger[idx])[lvl][dit()].listItems();

        gdnv_parts.clear();

        // Store the diffusion hop, and add the godunov particles
        for (ListIterator<ItoParticle> lit(ItoParticles); lit.ok(); ++lit) {
          ItoParticle& p = lit();

          const Real&     weight = p.weight();
          const RealVect& Xk     = p.oldPosition();
#if 0
	  const RealVect& hop  = p.runtimeVector(0);
	  const RealVect& Vk   = p.runtimeVector(1);

          // Move particle.
          const RealVect pos = Xk + 0.5 * a_dt * f * Vk + g * hop;
          gdnv_parts.add(PointParticle(pos, weight));
#else
          MayDay::Error("ItoPlasmaGodunovStepper::trapezoidalPredictor -- preTrapezoidalCorrector");
#endif
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::trapezoidalCorrector(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::trapezoidalCorrector");
  if (m_verbosity > 5) {
    pout() << m_name + "::trapezoidalCorrector" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;
    const Real g = diffusive ? 1.0 : 0.0;

    if (mobile || diffusive) {

      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
        ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

        for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

          List<ItoParticle>& particleList = particles[dit()].listItems();

          for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
            ItoParticle& p = lit();

            const RealVect& Xk = p.oldPosition();
#if 0
	    const RealVect& hop = p.runtimeVector(0);
	    const RealVect& Vk  = p.runtimeVector(1);
            const RealVect& Vk1 = p.velocity();

            p.position() = Xk + 0.5 * f * a_dt * (Vk + Vk1) + g * hop;
#else
            MayDay::Error("ItoPlasmaGodunovStepper::trapezoidalCorrector -- runtime stuff");
#endif
          }
        }
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
