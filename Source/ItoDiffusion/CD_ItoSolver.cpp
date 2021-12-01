/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoSolver.cpp
  @brief  Implementation of CD_ItoSolver.H
  @author Robert Marskar
*/

// Std includes
#include <chrono>

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>
#include <ParticleIO.H>

// Our includes
#include <CD_SimpleItoParticle.H>
#include <CD_ItoSolver.H>
#include <CD_DataOps.H>
#include <CD_ParticleOps.H>
#include <CD_NamespaceHeader.H>

#define ITO_DEBUG 0

constexpr int ItoSolver::m_comp;
constexpr int ItoSolver::m_nComp;

ItoSolver::ItoSolver() {
  CH_TIME("ItoSolver::ItoSolver");

  // Default settings
  m_verbosity      = -1           ;
  m_name          = "ItoSolver"  ;
  m_className     = "ItoSolver"  ;
  m_realm         = Realm::primal;
  m_phase         = phase::gas   ;
  m_checkpointing = WhichCheckpoint::Particles;

  m_mobilityInterpolation = WhichMobilityInterpolation::Direct;
}

ItoSolver::~ItoSolver() {

}

std::string ItoSolver::getName() const {
  CH_TIME("ItoSolver::getName");
  
  return m_name;
}

const std::string ItoSolver::getRealm() const {
  CH_TIME("ItoSolver::getRealm");
  
  return m_realm;
}

void ItoSolver::setRealm(const std::string a_realm) {
  CH_TIME("ItoSolver::setRealm");
  
  m_realm = a_realm;
}

const RefCountedPtr<ItoSpecies>& ItoSolver::getSpecies() const {
  CH_TIME("ItoSolver::getSpecies");
  if(m_verbosity > 5){
    pout() << m_name + "::getSpecies" << endl;
  }

  CH_assert(!m_species.isNull());
  
  return m_species;
}

void ItoSolver::parseOptions() {
  CH_TIME("ItoSolver::parseOptions");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseOptions" << endl;
  }

  this->parseSuperParticles();        
  this->parseRng();
  this->parsePlotVariables();
  this->parseDeposition();
  this->parseBisectStep();
  this->parsePvrBuffer();
  this->parseDiffusionHop();
  this->parseRedistribution();
  this->parseDivergenceComputation();
  this->parseCheckpointing();
}

void ItoSolver::parseRuntimeOptions() {
  CH_TIME("ItoSolver::parseRuntimeOptions");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }

  this->parseSuperParticles();
  this->parseRng();
  this->parsePlotVariables();
  this->parseDeposition();
  this->parseBisectStep();
  this->parseDiffusionHop();
  this->parseRedistribution();
  this->parseDivergenceComputation();
  this->parseCheckpointing();
}

void ItoSolver::parseSuperParticles() {
  CH_TIME("ItoSolver::parseSuperParticles");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseSuperParticles" << endl;
  }

  // Seed the RNG
  ParmParse pp(m_className.c_str());
  pp.get("kd_direction", m_directionKD);

  m_directionKD = min(m_directionKD, SpaceDim-1);
}

void ItoSolver::parseRng() {
  CH_TIME("ItoSolver::parseRng");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseRng" << endl;
  }

  // Seed the RNG
  ParmParse pp(m_className.c_str());
  pp.get("seed",       m_rngSeed);
  pp.get("normal_max", m_normalDistributionTruncation);

  // Use a random seed if the input was < 0.
  if(m_rngSeed < 0) { 
    m_rngSeed = std::chrono::system_clock::now().time_since_epoch().count();
  }

  // Initialize random number generator and distributions
  constexpr Real zero = 0.0;
  constexpr Real one  = 1.0;
  
  m_rng                   = std::mt19937_64(m_rngSeed);
  m_uniformDistribution01 = std::uniform_real_distribution<Real>( zero, one       ); // <- Uniform real distribution from [ 0   1]
  m_uniformDistribution11 = std::uniform_real_distribution<Real>(-one , one       ); // <- Uniform real distribution from [-1, -1]
  m_uniformDistribution0d = std::uniform_int_distribution<int>  ( 0   , SpaceDim-1); // <- Uniform integer distribution from [0, SpaceDim-1]
  m_normalDistribution01  = std::normal_distribution<Real>      ( zero, one       ); // <- Normal distribution with mean value of zero and standard deviation of one. 
}

void ItoSolver::parsePlotVariables() {
  CH_TIME("McPhoto::parsePlotVariables");
  if(m_verbosity > 5) {
    pout() << m_name + "::parsePlotVariables" << endl;
  }

  m_plotPhi             = false;
  m_plotVelocity        = false;
  m_plotDiffCo          = false;
  m_plotParticles       = false;
  m_plotParticlesEB     = false;
  m_plotParticlesDomain = false;
  m_plotParticlesSource = false;
  m_plotEnergyDensity   = false;
  m_plotAverageEnergy   = false;

  ParmParse pp(m_className.c_str());
  
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++) {
    if(     str[i] == "phi")            m_plotPhi             = true;
    else if(str[i] == "vel")            m_plotVelocity        = true;
    else if(str[i] == "dco")            m_plotDiffCo          = true;
    else if(str[i] == "part")           m_plotParticles       = true;
    else if(str[i] == "eb_part")        m_plotParticlesEB     = true;
    else if(str[i] == "dom_part")       m_plotParticlesDomain = true;
    else if(str[i] == "src_part")       m_plotParticlesSource = true;
    else if(str[i] == "energy_density") m_plotEnergyDensity   = true;
    else if(str[i] == "average_energy") m_plotAverageEnergy   = true;
  }
}

void ItoSolver::parseDeposition() {
  CH_TIME("ItoSolver::parseRng");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseRng" << endl;
  }

  ParmParse pp(m_className.c_str());
  std::string str;

  // Deposition for particle-mesh operations
  pp.get("deposition", str);
  if(str == "ngp") {
    m_deposition = DepositionType::NGP;
  }
  else if(str == "cic") {
    m_deposition = DepositionType::CIC;
  }
  else if(str == "tsc") {
    m_deposition = DepositionType::TSC;
  }
  else if(str == "w4") {
    m_deposition = DepositionType::W4;
  }
  else{
    MayDay::Abort("ItoSolver::parseDeposition - unknown interpolant requested");
  }

  // Deposition for plotting only
  pp.get("plot_deposition", str);

  if(str == "ngp") {
    m_plotDeposition = DepositionType::NGP;
  }
  else if(str == "cic") {
    m_plotDeposition = DepositionType::CIC;
  }
  else if(str == "tsc") {
    m_plotDeposition = DepositionType::TSC;
  }
  else if(str == "w4") {
    m_plotDeposition = DepositionType::W4;
  }
  else{
    MayDay::Abort("ItoSolver::parseDeposition - unknown interpolant requested");
  }

  // Mobility interpolation.
  pp.get("WhichMobilityInterpolation",str);

  if(str == "mobility") {
    m_mobilityInterpolation = WhichMobilityInterpolation::Direct;
  }
  else if(str == "velocity") {
    m_mobilityInterpolation = WhichMobilityInterpolation::Velocity;
  }
  else{
    MayDay::Abort("ItoSolver::parseDeposition - unknown interpolation method for mobility");
  }

  pp.get("irr_ngp_deposition", m_forceIrregDepositionNGP);
  pp.get("irr_ngp_interp",     m_forceIrregInterpolationNGP);
}

void ItoSolver::parseBisectStep() {
  CH_TIME("ItoSolver::parseBisectStep");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseBisectStep" << endl;
  }

  ParmParse pp(m_className.c_str());
  pp.get("bisect_step", m_bisectionStep);
}

void ItoSolver::parsePvrBuffer() {
  CH_TIME("ItoSolver::parsePvrBuffer");
  if(m_verbosity > 5) {
    pout() << m_name + "::parsePvrBuffer" << endl;
  }

  ParmParse pp(m_className.c_str());
  pp.get("pvr_buffer",  m_pvrBuffer);
  pp.get("halo_buffer", m_haloBuffer);

  std::string str;
  pp.get("halo_deposition", str);
  if(str == "ngp") {
    m_forceHaloNGP = true;
  }
  else if(str == "native") {
    m_forceHaloNGP = false;
  }
  else{
    MayDay::Abort("ItoSolver::parsePvrBuffer - unknown argument to 'halo_deposition'");
  }
}

void ItoSolver::parseDiffusionHop() {
  CH_TIME("ItoSolver::parseDiffusionHop");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseDiffusionHop" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("max_diffusion_hop", m_maxDiffusionHop);
}

void ItoSolver::parseRedistribution() {
  CH_TIME("ItoSolver::parseRedistribution");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseRedistribution" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("redistribute", m_useRedistribution);
}

void ItoSolver::parseDivergenceComputation() {
  CH_TIME("ItoSolver::parseDivergenceComputation");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseDivergenceComputation" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("blend_conservation", m_blendConservation);
}

void ItoSolver::parseCheckpointing() {
  CH_TIME("ItoSolver::parseCheckpointing");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseCheckpointing" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("checkpointing", str);
  pp.get("ppc_restart", m_restartPPC);
  if(str == "particles") {
    m_checkpointing = WhichCheckpoint::Particles;
  }
  else if(str == "numbers") {
    m_checkpointing = WhichCheckpoint::Numbers;
  }
  else{
    MayDay::Abort("ItoSolver::parseCheckpointing - unknown checkpointing method requested");
  }
}

Vector<std::string> ItoSolver::getPlotVariableNames() const {
  CH_TIME("ItoSolver::getPlotVariableNames");
  if(m_verbosity > 5) {
    pout() << m_name + "::getPlotVariableNames" << endl;
  }

  Vector<std::string> names(0);
  
  if(m_plotPhi) {
    names.push_back(m_name + " phi");
  }
  if(m_plotDiffCo && m_isDiffusive) {
    names.push_back(m_name + " diffusion_coefficient");
  }
  if(m_plotVelocity && m_isMobile) {
    names.push_back("x-Velocity " + m_name);
    names.push_back("y-Velocity " + m_name);
    if(SpaceDim == 3) {
      names.push_back("z-Velocity " + m_name);
    }
  }
  if(m_plotParticles)         names.push_back(m_name + " particles");
  if(m_plotParticlesEB)      names.push_back(m_name + " eb_particles");
  if(m_plotParticlesDomain)  names.push_back(m_name + " domain_particles");
  if(m_plotParticlesSource)  names.push_back(m_name + " source_particles");
  if(m_plotParticlesCovered) names.push_back(m_name + " covered_particles");
  if(m_plotEnergyDensity)    names.push_back(m_name + " energy * phi");
  if(m_plotAverageEnergy)    names.push_back(m_name + " average_energy");

  return names;
}

int ItoSolver::getNumberOfPlotVariables() const {
  CH_TIME("ItoSolver::getNumberOfPlotVariables");
  if(m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }

  int numPlotVars = 0;
  
  if(m_plotPhi                      ) numPlotVars += 1       ;
  if(m_plotDiffCo   && m_isDiffusive) numPlotVars += 1       ;
  if(m_plotVelocity && m_isMobile   ) numPlotVars += SpaceDim;
  if(m_plotParticles                ) numPlotVars += 1       ;
  if(m_plotParticlesEB              ) numPlotVars += 1       ;
  if(m_plotParticlesDomain          ) numPlotVars += 1       ;
  if(m_plotParticlesSource          ) numPlotVars += 1       ;
  if(m_plotParticlesCovered         ) numPlotVars += 1       ;
  if(m_plotEnergyDensity            ) numPlotVars += 1       ;
  if(m_plotAverageEnergy            ) numPlotVars += 1       ;

  return numPlotVars;
}

int ItoSolver::getPVRBuffer() const {
  CH_TIME("ItoSolver::getPVRBuffer");
  if(m_verbosity > 5) {
    pout() << m_name + "::getPVRBuffer" << endl;
  }

  return m_pvrBuffer;
}

int ItoSolver::getHaloBuffer() const {
  CH_TIME("ItoSolver::getHaloBuffer");
  if(m_verbosity > 5) {
    pout() << m_name + "::getHaloBuffer" << endl;
  }

  return m_haloBuffer;
}

void ItoSolver::setPVRBuffer(const int a_pvrBuffer) {
  CH_TIME("ItoSolver::setPVRBuffer");
  if(m_verbosity > 5) {
    pout() << m_name + "::setPVRBuffer" << endl;
  }

  CH_assert(a_pvrBuffer >= 0);

  m_pvrBuffer  = a_pvrBuffer;
  m_haloBuffer = 0;
}

void ItoSolver::setHaloBuffer(const int a_haloBuffer)  {
  CH_TIME("ItoSolver::setHaloBuffer");
  if(m_verbosity > 5) {
    pout() << m_name + "::setHaloBuffer" << endl;
  }

  CH_assert(a_haloBuffer >= 0);

  m_haloBuffer = a_haloBuffer;
  m_pvrBuffer  = 0;
}


size_t ItoSolver::getNumParticles(const WhichContainer a_whichContainer, const bool a_localOnly) const {
  CH_TIME("ItoSolver::getNumParticles(WhichContainer, bool)");
  if(m_verbosity > 5) {
    pout() << m_name + "::getNumParticles(WhichContainer, bool)" << endl;
  }

  const ParticleContainer<ItoParticle>& particles = m_particleContainers.at(a_whichContainer);

  size_t N;  
  if(a_localOnly) {
    N = particles.getNumberOfValidParticesLocal();
  }
  else{
    N = particles.getNumberOfValidParticesGlobal();
  }

  return N;
}

void ItoSolver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry) {
  CH_TIME("ItoSolver::setComputationalGeometry");
  if(m_verbosity > 5) {
    pout() << m_name + "::setComputationalGeometry" << endl;
  }

  CH_assert(!a_computationalGeometry.isNull());
  
  m_computationalGeometry = a_computationalGeometry;
}

void ItoSolver::setAmr(const RefCountedPtr<AmrMesh>& a_amr) {
  CH_TIME("ItoSolver::setAmr");
  if(m_verbosity > 5) {
    pout() << m_name + "::setAmr" << endl;
  }

  CH_assert(!a_amr.isNull());

  m_amr = a_amr;
}

void ItoSolver::registerOperators() const {
  CH_TIME("ItoSolver::registerOperators");
  if(m_verbosity > 5) {
    pout() << m_name + "::registerOperators" << endl;
  }

  if(m_amr.isNull()) {
    MayDay::Abort("CdrSolver::registerOperators - need to set AmrMesh!");
  }
  else{
    m_amr->registerOperator(s_eb_coar_ave,     m_realm, m_phase);
    m_amr->registerOperator(s_eb_fill_patch,   m_realm, m_phase);
    m_amr->registerOperator(s_noncons_div,     m_realm, m_phase);
    m_amr->registerOperator(s_particle_mesh,   m_realm, m_phase);
    if(m_useRedistribution) {
      m_amr->registerOperator(s_eb_redist,  m_realm, m_phase);
    }

    // Register mask if using halo deposition
    if(m_haloBuffer > 0) {
      CH_assert(m_pvrBuffer == 0);
      
      m_amr->registerMask(s_particle_halo, m_haloBuffer, m_realm);
    }
  }
}

void ItoSolver::setPhase(const phase::which_phase a_phase) {
  CH_TIME("ItoSolver::setPhase");
  if(m_verbosity > 5) {
    pout() << m_name + "::setPhase" << endl;
  }

  m_phase = a_phase;
}

void ItoSolver::setVerbosity(const int a_verbosity) {
  CH_TIME("ItoSolver::setVerbosity");
  
  m_verbosity = a_verbosity;
  
  if(m_verbosity > 5) {
    pout() << m_name + "::setVerbosity" << endl;
  }
}

void ItoSolver::setTime(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("ItoSolver::setTime");
  if(m_verbosity > 5) {
    pout() << m_name + "::setTime" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt  ;
}

void ItoSolver::initialData() {
  CH_TIME("ItoSolver::initialData");
  if(m_verbosity > 5) {
    pout() << m_name + "::initialData" << endl;
  }

  CH_assert(!m_species.isNull());

  // TLDR: This function will fetch the initial particles from the species and deposit them on the mesh. In most cases the various MPI ranks
  //       will have drawn a different set of initial particles (the only sane way to do it) and so those particles are put directly in
  //       the 'bulk' particle container. After that we remove the particles that fell inside the EB and deposit the particles on the mesh. 

  ParticleContainer<ItoParticle>& bulkParticles = m_particleContainers.at(WhichContainer::Bulk);
  bulkParticles.clearParticles();
  bulkParticles.addParticles(m_species->getInitialParticles());

  constexpr Real tolerance = 0.0;

  // Add particles, remove the ones that are inside the EB, and then depsit
  this->removeCoveredParticles(bulkParticles, EbRepresentation::ImplicitFunction, tolerance);  // Remove particles that are less than tolerance away from the EB
  this->depositParticles<ItoParticle, &ItoParticle::mass>(m_phi, bulkParticles, m_deposition); // Deposit particles on the mesh. 
}

void ItoSolver::computeLoads(Vector<long int>& a_loads, const DisjointBoxLayout& a_dbl, const int a_level) {
  CH_TIME("ItoSolver::computeLoads");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeLoads" << endl;
  }

  const ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  a_loads.resize(a_dbl.size(), 0L);
  for (DataIterator dit(a_dbl); dit.ok(); ++dit) {
    a_loads[dit().intCode()] = particles[a_level][dit()].numItems();
  }

  // If using MPI, we must gather the loads. 
#ifdef CH_MPI
  int count = a_loads.size();
  Vector<long int> tmp(count);
  MPI_Allreduce(&(a_loads[0]),&(tmp[0]), count, MPI_LONG, MPI_SUM, Chombo_MPI::comm);
  a_loads = tmp;
#endif
}

void ItoSolver::removeCoveredParticles(const EbRepresentation a_representation, const Real a_tol) {
  CH_TIME("ItoSolver::removeCoveredParticles(EbRepresentation, tolerance)");
  if(m_verbosity > 5) {
    pout() << m_name + "::removeCoveredParticles(EbRepresentation, tolerance)" << endl;
  }

  this->removeCoveredParticles(WhichContainer::Bulk, a_representation, a_tol);
}

void ItoSolver::removeCoveredParticles(const WhichContainer a_container, const EbRepresentation a_representation, const Real a_tol) {
  CH_TIME("ItoSolver::removeCoveredParticles(WhichContainer, EbRepresentation, Real)");
  if(m_verbosity > 5) {
    pout() << m_name + "::removeCoveredParticles(WhichContainer, EbRepresentation, Real)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);
  
  this->removeCoveredParticles(particles, a_representation, a_tol);
}

void ItoSolver::removeCoveredParticles(ParticleContainer<ItoParticle>& a_particles, const EbRepresentation a_representation, const Real a_tol) const {
  CH_TIME("ItoSolver::removeCoveredParticles(particles, EbRepresentation)");
  if(m_verbosity > 5) {
    pout() << m_name + "::removeCoveredParticles(particles, EbRepresentation)" << endl;
  }

  switch(a_representation) {
  case EbRepresentation::ImplicitFunction:
    m_amr->removeCoveredParticlesIF(a_particles, m_phase, a_tol);
    break;
  case EbRepresentation::Discrete:
    m_amr->removeCoveredParticlesDiscrete(a_particles, m_phase, a_tol);
    break;
  case EbRepresentation::Voxel:
    m_amr->removeCoveredParticlesVoxels(a_particles, m_phase);
    break;
  default:
    MayDay::Error("ItoSolver::removeCoveredParticles - unsupported EB representation requested");
  }
}

void ItoSolver::transferCoveredParticles(const EbRepresentation a_representation, const Real a_tol) {
  CH_TIME("ItoSolver::transferCoveredParticles(EbRepresentation, Real)");
  if(m_verbosity > 5) {
    pout() << m_name + "::transferCoveredParticles(EbRepresentation, Real)" << endl;
  }

  this->transferCoveredParticles(WhichContainer::Bulk, WhichContainer::Covered, a_representation, a_tol);
}

void ItoSolver::transferCoveredParticles(const WhichContainer a_containerFrom, const WhichContainer a_containerTo, const EbRepresentation a_representation, const Real a_tol) {
  CH_TIME("ItoSolver::transferCoveredParticles(WhichContainer, WhichContainer, EbRepresentation, Real)");
  if(m_verbosity > 5) {
    pout() << m_name + "::transferCoveredParticles(WhichContainer, WhichContainer, EbRepresentation, Real)" << endl;
  }

  ParticleContainer<ItoParticle>& particlesFrom = this->getParticles(a_containerFrom);
  ParticleContainer<ItoParticle>& particlesTo   = this->getParticles(a_containerTo);

  this->transferCoveredParticles(particlesFrom, particlesTo, a_representation, a_tol);
}

void ItoSolver::transferCoveredParticles(ParticleContainer<ItoParticle>& a_particlesFrom,
					 ParticleContainer<ItoParticle>& a_particlesTo,
					 const EbRepresentation          a_representation,
					 const Real                      a_tol) const {
  CH_TIME("ItoSolver::transferCoveredParticles(ParticleContainer, ParticleContainer, EbRepresentation, Real)");
  if(m_verbosity > 5) {
    pout() << m_name + "::transferCoveredParticles(ParticleContainer, ParticleContainer, EbRepresentation, Real)" << endl;
  }

  switch(a_representation) {
  case EbRepresentation::ImplicitFunction:
    m_amr->transferCoveredParticlesIF(a_particlesFrom, a_particlesTo, m_phase, a_tol);
    break;
  case EbRepresentation::Discrete:
    m_amr->transferCoveredParticlesDiscrete(a_particlesFrom, a_particlesTo, m_phase, a_tol);    
    break;
  case EbRepresentation::Voxel:
    m_amr->transferCoveredParticlesVoxels(a_particlesFrom, a_particlesTo, m_phase);
    break;        
  default:
    MayDay::Error("ItoSolver::transferCoveredParticles -- logic bust");
  }
}

void ItoSolver::intersectParticles(const EbRepresentation a_ebRepresentation, const bool a_deleteParticles) {
  CH_TIME("ItoSolver::intersectParticles(EbRepresentation, bool)");
  if(m_verbosity > 5) {
    pout() << m_name + "::intersectParticles(EbRepresentation, bool)" << endl;
  }

  this->intersectParticles(WhichContainer::Bulk, WhichContainer::EB, WhichContainer::Domain, a_ebRepresentation, a_deleteParticles);
}

void ItoSolver::intersectParticles(const WhichContainer   a_particles,
				   const WhichContainer   a_ebParticles,
				   const WhichContainer   a_domainParticles,
				   const EbRepresentation a_ebRepresentation,				     
				   const bool             a_deleteParticles) {
  CH_TIME("ItoSolver::intersectParticles(WhichContainerx3, EbRepresentation, bool)");
  if(m_verbosity > 5) {
    pout() << m_name + "::intersectParticles(WhichContainerx3, EbRepresentation, bool)" << endl;
  }

  ParticleContainer<ItoParticle>& particles       = this->getParticles(a_particles);
  ParticleContainer<ItoParticle>& ebParticles     = this->getParticles(a_ebParticles);
  ParticleContainer<ItoParticle>& domainParticles = this->getParticles(a_domainParticles);

  this->intersectParticles(particles, ebParticles, domainParticles, a_ebRepresentation, a_deleteParticles);
}


void ItoSolver::intersectParticles(ParticleContainer<ItoParticle>& a_particles,
				   ParticleContainer<ItoParticle>& a_ebParticles,
				   ParticleContainer<ItoParticle>& a_domainParticles,
				   const EbRepresentation          a_ebRepresentation,
				   const bool                      a_deleteParticles) {
  CH_TIME("ItoSolver::intersectParticles(ParticleContainerx3, EbRepresentation, bool)");
  if(m_verbosity > 5) {
    pout() << m_name + "::intersectParticles(ParticleContainerx3, EbRepresentation, bool)" << endl;
  }

  CH_assert(!a_particles.      isCellSorted());
  CH_assert(!a_ebParticles.    isCellSorted());
  CH_assert(!a_domainParticles.isCellSorted());  

  switch(a_ebRepresentation) {
  case EbRepresentation::ImplicitFunction:
    this->intersectParticlesIF(a_particles, a_ebParticles, a_domainParticles, a_deleteParticles);
    break;
  default:
    MayDay::Error("ItoSolver::intersectParticles - unsupported EB representation requested");
  }
}

void ItoSolver::intersectParticlesIF(ParticleContainer<ItoParticle>& a_particles,
				     ParticleContainer<ItoParticle>& a_eb_particles,
				     ParticleContainer<ItoParticle>& a_domain_particles,
				     const bool                        a_delete) {
  CH_TIME("ItoSolver::intersectParticlesIF(container, container, container, bool)");
  if(m_verbosity > 5) {
    pout() << m_name + "::intersectParticlesIF(container, container, container, bool)" << endl;
  }

  const RealVect prob_lo = m_amr->getProbLo();
  const RealVect prob_hi = m_amr->getProbHi();
  const Real     SAFETY  = 1.E-6;

  // This is the implicit function used for intersection tests
  const RefCountedPtr<BaseIF>& impfunc = m_amr->getBaseImplicitFunction(m_phase);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    for (DataIterator dit = m_amr->getGrids(m_realm)[lvl]; dit.ok(); ++dit) {

      const Real     dx      = m_amr->getDx()[lvl];
      const EBISBox& ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[lvl][dit()];

      List<ItoParticle>& particles    = a_particles[lvl][dit()].listItems();
      List<ItoParticle>& ebParticles  = a_eb_particles[lvl][dit()].listItems();
      List<ItoParticle>& domParticles = a_domain_particles[lvl][dit()].listItems();

      ebParticles.clear();
      domParticles.clear();

      for (ListIterator<ItoParticle> lit(particles); lit.ok(); ++lit) {
	ItoParticle& p = lit();

	const RealVect newPos  = p.position();
	const RealVect oldPos  = p.oldPosition();
	const RealVect path    = newPos - oldPos;
	const Real     pathLen = path.vectorLength();

	// Check if we should check of different types of boundary intersections. These are checp initial tests that allow
	// us to skip intersection tests for some Photons.
	bool checkEB  = false;
	bool checkDom = false;

	if(!impfunc.isNull()) {
	  checkEB = true;
	}
	for (int dir = 0; dir < SpaceDim; dir++) {
	  if(newPos[dir] < prob_lo[dir] || newPos[dir] > prob_hi[dir]) { // Checks if we crossed a domain boundary. 
	    checkDom = true; 
	  }
	}

	// Must do intersection test on at least one of these. 
	if(checkEB || checkDom) { 
	  Real dom_s = 1.E99;
	  Real eb_s  = 1.E99;

	  bool contact_domain = false;
	  bool contact_eb     = false;
	      
	  if(checkDom) contact_domain = ParticleOps::domainIntersection(oldPos, newPos, prob_lo, prob_hi, dom_s);
#if 0
	  if(checkEB)  contact_eb     = ParticleOps::ebIntersectionBisect(impfunc, oldPos, newPos, dx, eb_s);
#else
	  if(checkEB)  contact_eb     = ParticleOps::ebIntersectionRaycast(impfunc, oldPos, newPos, 1.E-10*dx, eb_s);
#endif
	  
	  if(contact_eb || contact_domain) { // Particle trajectory crossed something. 
	    if(eb_s < dom_s) { // It was the EB first. 
	      p.position() = oldPos + eb_s*path;
	      if(a_delete) {
		ebParticles.transfer(lit);
	      }
	      else{
		ebParticles.add(lit());
	      }
	    }
	    else{ // It was the domain boundary. 
	      p.position() = oldPos + Max(0.0,dom_s-SAFETY)*path;
	      if(a_delete) {
		domParticles.transfer(lit);
	      }
	      else{
		domParticles.add(lit());
	      }
	    }
	  }
	}
      }
    }
  }
}

void ItoSolver::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) {
  CH_TIME("ItoSolver::regrid");
  if(m_verbosity > 5) {
    pout() << m_name + "::regrid" << endl;
  }

  CH_assert(a_lmin           >= 0);
  CH_assert(a_oldFinestLevel >= 0);
  CH_assert(a_newFinestLevel >= 0);    

  // Reallocate mesh data. 
  m_amr->reallocate(m_phi,              m_phase, a_lmin);
  m_amr->reallocate(m_scratch,          m_phase, a_lmin);
  m_amr->reallocate(m_depositionNC,     m_phase, a_lmin);
  m_amr->reallocate(m_massDiff,         m_phase, a_lmin);
  
  // Only allocate memory if we have advection. 
  if(m_isMobile) {
    m_amr->reallocate(m_mobilityFunction, m_phase, a_lmin);    
    m_amr->reallocate(m_velocityFunction, m_phase, a_lmin);
  }
  else{
    m_amr->allocatePointer(m_mobilityFunction);    
    m_amr->allocatePointer(m_velocityFunction);
  }

  // Only allocate memory if we have diffusion.
  if(m_isDiffusive) {
    m_amr->reallocate(m_diffusionFunction, m_phase, a_lmin);
  }
  else{
    m_amr->allocatePointer(m_diffusionFunction);
  }

  // Regrid particle containers. 
  const Vector<DisjointBoxLayout>& grids   = m_amr->getGrids(m_realm);
  const Vector<ProblemDomain>&     domains = m_amr->getDomains();
  const Vector<Real>&              dx      = m_amr->getDx();
  const Vector<int>&               refRat  = m_amr->getRefinementRatios();

  for (auto& container : m_particleContainers) {
    ParticleContainer<ItoParticle>& particles = container.second;
    
    particles.regrid(grids, domains, dx, refRat, a_lmin, a_newFinestLevel);
  }
}

void ItoSolver::setSpecies(const RefCountedPtr<ItoSpecies>& a_species) {
  CH_TIME("ItoSolver::setSpecies");
  if(m_verbosity > 5) {
    pout() << m_name + "::setSpecies" << endl;
  }
  
  m_species     = a_species;
  m_name        = a_species->getName();
  m_isDiffusive = m_species->isDiffusive();
  m_isMobile    = m_species->isMobile();
}

void ItoSolver::allocateInternals() {
  CH_TIME("ItoSolver::allocateInternals");
  if(m_verbosity > 5) {
    pout() << m_name + "::allocateInternals" << endl;
  }

  CH_assert(!m_species.isNull());
  
  const int ncomp = 1;

  m_amr->allocate(m_phi,              m_realm, m_phase, ncomp); // Mesh data -- always allocate it. 
  m_amr->allocate(m_scratch,          m_realm, m_phase, ncomp); // Scratch data -- needed because of IO when we deposit particles (should always allocate it). 
  m_amr->allocate(m_depositionNC,     m_realm, m_phase, ncomp); // This is BaseIVFAB data so it has a low memory overhead. It doesn't hurt to allocate it. 
  m_amr->allocate(m_massDiff,         m_realm, m_phase, ncomp); // This is BaseIVFAB data so it has a low memory overhead. It doesn't hurt to allocate it. 

  // Only allocate memory for velocity if we actually have a mobile solver
  if(m_isMobile) {
    m_amr->allocate(m_mobilityFunction, m_realm, m_phase, ncomp); //     
    m_amr->allocate(m_velocityFunction, m_realm, m_phase, SpaceDim);
  }
  else{
    m_amr->allocatePointer(m_mobilityFunction);
    m_amr->allocatePointer(m_velocityFunction);
  }

  // Only allocate memory if we actually have a diffusion solver
  if(m_isDiffusive) {
    m_amr->allocate(m_diffusionFunction, m_realm, m_phase, 1);
  }
  else{
    m_amr->allocatePointer(m_diffusionFunction);
  }

  m_particleContainers.emplace(WhichContainer::Bulk,    ParticleContainer<ItoParticle>());
  m_particleContainers.emplace(WhichContainer::EB,      ParticleContainer<ItoParticle>());
  m_particleContainers.emplace(WhichContainer::Domain,  ParticleContainer<ItoParticle>());
  m_particleContainers.emplace(WhichContainer::Source,  ParticleContainer<ItoParticle>());
  m_particleContainers.emplace(WhichContainer::Covered, ParticleContainer<ItoParticle>());
  m_particleContainers.emplace(WhichContainer::Scratch, ParticleContainer<ItoParticle>());

  for (auto& container : m_particleContainers) {
    m_amr->allocate(container.second, m_pvrBuffer, m_realm);
  }
}

#ifdef CH_USE_HDF5  
void ItoSolver::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("ItoSolver::writeCheckpointLevel");
  if(m_verbosity > 5) {
    pout() << m_name + "::writeCheckpointLevel" << endl;
  }

  // Write state. 
  write(a_handle, *m_phi[a_level], m_name);

  // Write particles.
  switch(m_checkpointing){
  case WhichCheckpoint::Particles:
    this->writeCheckPointLevelParticles(a_handle, a_level);
    break;
  case WhichCheckpoint::Numbers:
    this->writeCheckPointLevelFluid(a_handle, a_level);
    break;
  default:
    MayDay::Error("ItoSolver::writeCheckpointLevel -- logic bust");
  }
}
#endif

#ifdef CH_USE_HDF5  
void ItoSolver::writeCheckPointLevelParticles(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("ItoSolver::writeCheckPointLevelParticles");
  if(m_verbosity > 5) {
    pout() << m_name + "::writeCheckPointLevelParticles" << endl;
  }

  // TLDR: This routine writes particles directly to the HDF5 file -- since the ItoParticle is quite large (112 bytes, ish), we put the most
  //       important ItoParticle fields into the SimpleItoParticle and checkpoint that particle type instead. This is what happens below.

  const std::string str = m_name + "_particlesP"; // I call this _particlesP to distinguish it from the "fluid" checkpointing method.

  // Set up a particle container with SimpleItoParticle -- this is ItoParticle's low-memory cousin. 
  ParticleContainer<SimpleItoParticle> lowMemoryParticles;
  m_amr->allocate(lowMemoryParticles,  m_pvrBuffer, m_realm);

  // Handle to the particles that will be checkpointed. 
  const ParticleContainer<ItoParticle>& myParticles = this->getParticles(WhichContainer::Bulk);

  // Make ItoParticle into SimpleItoParticle. This saves a ton of disk space.
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  for (DataIterator dit(dbl); dit.ok(); ++dit) {

    // Handle to list of particles low-memory particles in the current grid patch. 
    List<SimpleItoParticle>& simpleParticles = (lowMemoryParticles[a_level])[dit()].listItems();

    // Make the ItoParticle into a SimpleItoParticle for checkpointing purposes -- we only need mass, position, and energy.
    simpleParticles.clear();
    for (ListIterator<ItoParticle> lit(myParticles[a_level][dit()].listItems()); lit.ok(); ++lit) {
      const ItoParticle& p = lit();
      simpleParticles.append(SimpleItoParticle(p.mass(), p.position(), p.energy()));
    }
  }

  // Finally, write the particles to HDF5. 
  writeParticlesToHDF(a_handle, lowMemoryParticles[a_level], str);
}
#endif

#ifdef CH_USE_HDF5  
void ItoSolver::writeCheckPointLevelFluid(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("ItoSolver::writeCheckPointLevelFluid");
  if(m_verbosity > 5) {
    pout() << m_name + "::writeCheckPointLevelFluid" << endl;
  }

  // TLDR: This routine checkpoints the particle data using the number of particles in a grid cell. When the simulation is restarted, we read the
  //       number of particles per cell from the HDF5 file and re-initialize the particles. However, this function does NOT currently store the energy
  //       (or other parameters of interest) in the HDF5 file. Only the number of particles is available. I don't expect this function to be widely used
  //       by anyone.

  // Checkpointing name. 
  const std::string str = m_name + "_particlesF"; // I call this _particlesF to distinguish it from the "particle" checkpointing method.

  // Handles to relevant grid information
  const DisjointBoxLayout& dbl    = m_amr->getGrids     (m_realm         )[a_level]; // Boxes.
  const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[a_level]; // EB-information
  const RealVect           dx     = m_amr->getDx()[a_level]*RealVect::Unit;          // Grid-resolution.
  const RealVect           probLo = m_amr->getProbLo();                              // Lower-left corner of physical domain. 

  // Handle to the particles that will be checkpointed. 
  const ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  // Create transient storage that holds the particle numbers. 
  LevelData<EBCellFAB> particleNumbers(dbl, m_nComp, IntVect::Zero, EBCellFactory(ebisl));
  EBLevelDataOps::setVal(particleNumbers, 0.0);

  // Now go through the grid and add the number of particles in each cell
  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box cellBox = dbl[dit()];
    
    BaseFab<Real>& particleNumbersFAB = particleNumbers[dit()].getSingleValuedFAB(); // No multivalued cells please. 

    // Add the patch particles into cell-sorted particles. 
    BinFab<ItoParticle> cellSortedParticles(cellBox, dx, probLo);
    cellSortedParticles.addItems(particles[a_level][dit()].listItems());

    // Go through the patch and get the number of particles per cell. 
    for (BoxIterator bit(cellBox); bit.ok(); ++bit) {
      const IntVect iv = bit();

      const List<ItoParticle>& cellParticles = cellSortedParticles(iv, m_comp);      

      // Go through the particles in the current grid cell and set the number of particles. 
      particleNumbersFAB(iv, m_comp) = 0.0;
      for (ListIterator<ItoParticle> lit(cellParticles); lit.ok(); ++lit) {
	const ItoParticle& p = lit();
	
	particleNumbersFAB(iv, m_comp) += p.mass();
      }
    }
  }

  // Finally, write the particle numbers to HDF5. 
  write(a_handle, particleNumbers, str);
}
#endif

#ifdef CH_USE_HDF5  
void ItoSolver::readCheckpointLevel(HDF5Handle& a_handle, const int a_level) {
  CH_TIME("ItoSolver::readCheckpointLevel");
  if(m_verbosity > 5) {
    pout() << m_name + "::readCheckpointLevel" << endl;
  }

  // Read state vector
  read<EBCellFAB>(a_handle, *m_phi[a_level], m_name, m_amr->getGrids(m_realm)[a_level], Interval(0,0), false);

  // Instantiate the particles
  switch(m_checkpointing){
  case WhichCheckpoint::Particles:
    this->readCheckpointLevelParticles(a_handle, a_level);
    break;
  case WhichCheckpoint::Numbers:
    this->readCheckpointLevelFluid(a_handle, a_level);
    break;
  default:
    MayDay::Error("ItoSolver::readCheckpointLevel -- logic bust");
  }
}
#endif

#ifdef CH_USE_HDF5
void ItoSolver::readCheckpointLevelParticles(HDF5Handle& a_handle, const int a_level) {
  CH_TIME("ItoSolver::readCheckpointLevelParticles");
  if(m_verbosity > 5) {
    pout() << m_name + "::readCheckpointLevelParticles" << endl;
  }

  // TLDR: This function is the one that reads SimpleItoParticles from the checkpoint file and instantiates full ItoParticle's from that. Recalling
  //       writeCheckpointLevelParticles we only stored the mass, position, and energy of the particles. Here we read that information back in.   

  // This is the particle container that we will fill. 
  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);  

  CH_assert(m_checkpointing == WhichCheckpoint::Particles);
  CH_assert(!particles.isCellSorted()                    );

  const std::string str = m_name + "_particlesP";
    
  // Allocate storage for holding simple particle data and then read the particles into simpleParticles[a_level]. The other
  // grid levels are not touched. 
  Vector<RefCountedPtr<ParticleData<SimpleItoParticle> > > simpleParticles;
  m_amr->allocate(simpleParticles, m_realm);

  readParticlesFromHDF(a_handle, *simpleParticles[a_level], str);


  // Go through the particles we read from the file and make them into true ItoParticles. 
  for (DataIterator dit(m_amr->getGrids(m_realm)[a_level]); dit.ok(); ++dit) {
    List<ItoParticle>&             itoParticles       =   particles      [a_level] [dit()].listItems();
    const List<SimpleItoParticle>& simpleItoParticles = (*simpleParticles[a_level])[dit()].listItems();

    for (ListIterator<SimpleItoParticle> lit(simpleItoParticles); lit.ok(); ++lit) {
      const SimpleItoParticle& simpleParticle = lit();
      const ItoParticle        itoParticle    = ItoParticle(simpleParticle.mass(),
							    simpleParticle.position(),
							    RealVect::Zero,
							    0.0,
							    0.0,
							    simpleParticle.energy());
      
      itoParticles.append(itoParticle);
    }
  }  
}
#endif

#ifdef CH_USE_HDF5
void ItoSolver::readCheckpointLevelFluid(HDF5Handle& a_handle, const int a_level) {
  CH_TIME("ItoSolver::readCheckpointLevelFluid");
  if(m_verbosity > 5) {
    pout() << m_name + "::readCheckpointLevelFluid" << endl;
  }

  CH_assert(m_checkpointing == WhichCheckpoint::Numbers);

  const std::string str = m_name + "_particlesF";

  constexpr int      comp     = 0;
  constexpr int      numComp  = 1;
  constexpr int      numGhost = 0;
  const     Interval interv   = Interval(comp, comp);
  
  // Allocate some storage that we can read into. 
  EBAMRCellData particlesPerCell;
  m_amr->allocate(particlesPerCell, m_realm, m_phase, numComp, numGhost);

  read<EBCellFAB>(a_handle, *particlesPerCell[a_level], str, m_amr->getGrids(m_realm)[a_level], interv, false);
    
  // particlesPerCell holds the number of particles per cell -- call the other version which instantiates new particles from that. 
  this->drawNewParticles(*particlesPerCell[a_level], a_level, m_restartPPC);  
}
#endif

void ItoSolver::drawNewParticles(const LevelData<EBCellFAB>& a_particlesPerCell, const int a_level, const int a_newPPC) {
  CH_TIME("ItoSolver::drawNewParticles");
  if(m_verbosity > 5) {
    pout() << m_name + "::drawNewParticles" << endl;
  }

  // Handle to grid information. 
  const RealVect           probLo = m_amr->getProbLo();
  const Real               dx     = m_amr->getDx()[a_level];    
  const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[a_level];
  const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[a_level];

  // Particle container that we will fill. 
  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  // Go through each patch and instantiate new particles. 
  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box&           cellBox = dbl  [dit()];
    const EBISBox&       ebisbox = ebisl[dit()];
    const BaseFab<Real>& ppc     = a_particlesPerCell[dit()].getSingleValuedFAB();

    // This should draw new particles rather than append -- so clear out any old particles. 
    List<ItoParticle>& myParticles = particles[a_level][dit()].listItems();
    myParticles.clear();

    // Do regular cells
    for (BoxIterator bit(cellBox); bit.ok(); ++bit) {
      const IntVect iv = bit();

      // Do regular cells -- in these cells we only need to draw a random position somewhere inside the cubic cell. Easy.
      if(ebisbox.isRegular(iv)) {

	// Compute weights and remainder. This bit of code will take the number of physical particles and divide them into a_newPPC particles with
	// approximately equal weights. It is possible that one of the particles will have a larger particle weight than the others. 
	const unsigned long long numPhysicalParticles = (unsigned long long) llround(ppc(iv));

	unsigned long long computationalParticleWeight;
	unsigned long long computationalParticleNum;
	unsigned long long computationalParticleRemainder;
	DataOps::computeParticleWeights(computationalParticleWeight, computationalParticleNum, computationalParticleRemainder, numPhysicalParticles, a_newPPC);

	// Settings for drawing new particles. 
	const RealVect minLo = -0.5*RealVect::Unit;
	const RealVect minHi =  0.5*RealVect::Unit;
	const RealVect norma =      RealVect::Zero;
	const RealVect centr =      RealVect::Zero;
	const RealVect pos   = probLo + (RealVect(iv) + 0.5*RealVect::Unit)*dx;
	const Real kappa     = 1.0;

	// Now add the partices. If the remainder was > 0 we add another one with weight w + r
	for (unsigned long long i = 0; i < computationalParticleNum; i++) {
	  const RealVect particlePosition = this->randomPosition(pos, minLo, minHi, centr, norma, dx, kappa);
	  const Real     particleWeight   = (Real) computationalParticleWeight;
	  
	  myParticles.add(ItoParticle(particleWeight, particlePosition));
	}

	// Add the "remainder" particle. 
	if(computationalParticleRemainder > 0ULL) {
	  const RealVect particlePosition = this->randomPosition(pos, minLo, minHi, centr, norma, dx, kappa);
	  const Real     particleWeight   = (Real) (computationalParticleWeight + computationalParticleRemainder);

	  myParticles.add(ItoParticle(particleWeight, particlePosition));
	}
      }
    }

    // Do the same for irregular cells. This differs from the regular-cell case only in that the positions 
    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_level])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit) {
      const VolIndex& vof   = vofit();
      const IntVect   iv    = vof.gridIndex();
      const RealVect  cent  = ebisbox.bndryCentroid(vof);
      const RealVect  norm  = ebisbox.normal(vof);
      const RealVect  pos   = probLo + dx*(RealVect(iv) + 0.5*RealVect::Unit);
      const Real      kappa = ebisbox.volFrac(vof);

      const unsigned long long numPhysicalParticles = (unsigned long long) llround(ppc(iv));

      if(numPhysicalParticles > 0ULL){

	// No multi-valued cells please -- I don't know how to handle them. 
	CH_assert(!ebisbox.isMultiValued(iv));	
	
	// Compute a small box that encloses the cut-cell volume
	RealVect minLo = -0.5*RealVect::Unit;
	RealVect minHi =  0.5*RealVect::Unit;
	if(kappa < 1.0) {
	  DataOps::computeMinValidBox(minLo, minHi, norm, cent);
	}

	// Compute weights and remainder
	unsigned long long computationalParticleWeight;
	unsigned long long computationalParticleNum;
	unsigned long long computationalParticleRemainder;
	DataOps::computeParticleWeights(computationalParticleWeight, computationalParticleNum, computationalParticleRemainder, numPhysicalParticles, a_newPPC);

	// Now add the partices. If r > 0 we add another one with weight w + r
	for (unsigned long long i = 0; i < computationalParticleNum; i++) {
	  const RealVect particlePosition = this->randomPosition(pos, minLo, minHi, cent, norm, dx, kappa);
	  const Real particleWeight       = (Real) computationalParticleWeight;
	  
	  myParticles.add(ItoParticle(particleWeight, particlePosition));
	}

	if(computationalParticleRemainder > 0ULL) {
	  const RealVect particlePosition = this->randomPosition(pos, minLo, minHi, cent, norm, dx, kappa);
	  const Real particleWeight       = (Real) (computationalParticleWeight + computationalParticleRemainder);

	  myParticles.add(ItoParticle(particleWeight, particlePosition));
	}
      }
    }
  }
}

void ItoSolver::writePlotData(EBAMRCellData& a_output, int& a_comp) const {
  CH_TIME("ItoSolver::writePlotData");
  if(m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }

  CH_assert(a_comp >= 0);

  // Write phi
  if(m_plotPhi) {
    const Interval srcInterval(m_comp, m_comp);
    const Interval dstInterval(a_comp, a_comp);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      if(m_realm == a_output.getRealm()) {
	m_phi[lvl]->localCopyTo(srcInterval, *a_output[lvl], dstInterval);
      }
      else{
	m_phi[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
      }
    }
    //DataOps::setCoveredValue(a_output, a_comp, 0.0);
    a_comp++;
  }

  // Plot diffusion coefficient
  if(m_plotDiffCo && m_isDiffusive) {
    const Interval srcInterval(m_comp, m_comp);
    const Interval dstInterval(a_comp, a_comp);
    
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      if(m_realm == a_output.getRealm()) {
	m_diffusionFunction[lvl]->localCopyTo(srcInterval, *a_output[lvl], dstInterval);
      }
      else{
	m_diffusionFunction[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
      }
    }

    // Set covered diffusion coefficient to zero.
    DataOps::setCoveredValue(a_output, a_comp, 0.0);
    
    a_comp++;
  }

  // Write velocities
  if(m_plotVelocity && m_isMobile) {
    const Interval srcInterval(m_comp,          SpaceDim-1);
    const Interval dstInterval(a_comp, a_comp + SpaceDim-1);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      if(m_realm == a_output.getRealm()) {
	m_velocityFunction[lvl]->localCopyTo(srcInterval, *a_output[lvl], dstInterval);
      }
      else{
	m_velocityFunction[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
      }
    }

    // Set covered velocities to zero.
    for (int dir = 0; dir < SpaceDim; dir++) {
      DataOps::setCoveredValue(a_output, a_comp + dir, 0.0);
    }

    a_comp += SpaceDim;
  }

  // Plot various particle data holders.
  constexpr bool interpolateToCentroids = false;
  
  if(m_plotParticles) {
    this->depositParticles<ItoParticle, &ItoParticle::mass>(m_scratch, m_particleContainers.at(WhichContainer::Bulk), m_plotDeposition);
    this->writeData(a_output, a_comp, m_scratch, interpolateToCentroids);
  }
  if(m_plotParticlesEB) {
    this->depositParticles<ItoParticle, &ItoParticle::mass>(m_scratch, m_particleContainers.at(WhichContainer::EB), m_plotDeposition);
    this->writeData(a_output, a_comp, m_scratch, interpolateToCentroids);
  }
  if(m_plotParticlesDomain) {
    this->depositParticles<ItoParticle, &ItoParticle::mass>(m_scratch, m_particleContainers.at(WhichContainer::Domain), m_plotDeposition);
    this->writeData(a_output, a_comp, m_scratch, interpolateToCentroids);
  }
  if(m_plotParticlesSource) {
    this->depositParticles<ItoParticle, &ItoParticle::mass>(m_scratch, m_particleContainers.at(WhichContainer::Source), m_plotDeposition);
    this->writeData(a_output, a_comp, m_scratch, interpolateToCentroids);
  }
  if(m_plotEnergyDensity) {
    this->depositEnergyDensity(m_scratch, m_particleContainers.at(WhichContainer::Bulk), m_plotDeposition);
    this->writeData(a_output, a_comp, m_scratch, interpolateToCentroids);
  }
  if(m_plotAverageEnergy) {
    this->computeAverageEnergy(m_scratch, m_particleContainers.at(WhichContainer::Bulk));
    this->writeData(a_output, a_comp, m_scratch, interpolateToCentroids);
  }
}

void ItoSolver::writeData(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp) const {
  CH_TIME("ItoSolver::writeData");
  if(m_verbosity > 5) {
    pout() << m_name + "::writeData" << endl;
  }

  CH_assert(a_output.getRealm() == m_realm);

  const int numComp = a_data[0]->nComp();

  const Interval srcInterval(0, numComp-1);
  const Interval dstInterval(a_comp, a_comp + numComp - 1);

  // Copy data onto scratch
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, numComp);
  DataOps::copy(scratch, a_data);

  // Interp if we should
  if(a_interp) {
    m_amr->interpToCentroids(scratch, m_realm, phase::gas);
  }

  m_amr->averageDown(scratch, m_realm, m_phase);
  m_amr->interpGhost(scratch, m_realm, m_phase);

  // Copy into source data holder.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    if(m_realm == a_output.getRealm()) {
      scratch[lvl]->localCopyTo(srcInterval, *a_output[lvl], dstInterval);
    }
    else{
      scratch[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
    }
  }

  // Set covered value to zero.
  for (int comp = 0; comp < numComp; comp++){
    DataOps::setCoveredValue(a_output, a_comp + comp, 0.0);
  }

  a_comp += numComp;
}

void ItoSolver::depositConductivity(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles) const {
  CH_TIME("ItoSolver::depositConductivity(EBAMRCellData, ParticleContainer)");
  if(m_verbosity > 5) {
    pout() << m_name + "::depositConductivity(EBAMRCellData, ParticleContainer)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1     );
  CH_assert(!a_particles.isCellSorted());

  this->depositConductivity(a_phi, a_particles, m_deposition);
}

void ItoSolver::depositConductivity(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles, const DepositionType a_deposition) const {
  CH_TIME("ItoSolver::depositConductivity(EBAMRCellData, ParticleContainer, DepositionType)");
  if(m_verbosity > 5) {
    pout() << m_name + "::depositConductivity(EBAMRCellData, ParticleContainer, DepositionType)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1     );
  CH_assert(!a_particles.isCellSorted());  

  if(m_isMobile) {
    this->depositParticles<ItoParticle, &ItoParticle::conductivity>(a_phi, a_particles, a_deposition);
  }
  else{
    DataOps::setValue(a_phi, 0.0);
  }
}

void ItoSolver::depositDiffusivity(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles) const {
  CH_TIME("ItoSolver::depositDiffusivity(EBAMRCellData, ParticleContainer)");
  if(m_verbosity > 5) {
    pout() << m_name + "::depositDiffusivity(EBAMRCellData, ParticleContainer)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1     );
  CH_assert(!a_particles.isCellSorted());  

  this->depositDiffusivity(a_phi, a_particles, m_deposition);
}

void ItoSolver::depositDiffusivity(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles, const DepositionType a_deposition) const {
  CH_TIME("ItoSolver::depositDiffusivity(EBAMRCellData, ParticleContainer, DepositionType)");
  if(m_verbosity > 5) {
    pout() << m_name + "::depositDiffusivity(EBAMRCellData, ParticleContainer, DepositionType)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1     );
  CH_assert(!a_particles.isCellSorted());  

  this->depositParticles<ItoParticle, &ItoParticle::diffusivity>(a_phi, a_particles, a_deposition);
}

void ItoSolver::depositEnergyDensity(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles) const {
  CH_TIME("ItoSolver::depositEnergyDensity(EBAMRCellData, ParticleContainer)");
  if(m_verbosity > 5) {
    pout() << m_name + "::depositEnergyDensity(EBAMRCellData, ParticleContainer)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1     );
  CH_assert(!a_particles.isCellSorted());  

  this->depositEnergyDensity(a_phi, a_particles, m_deposition);
}

void ItoSolver::depositEnergyDensity(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles, const DepositionType a_deposition) const {
  CH_TIME("ItoSolver::depositEnergyDensity(EBAMRCellData, ParticleContainer, DepositionType)");
  if(m_verbosity > 5) {
    pout() << m_name + "::depositEnergyDensity(EBAMRCellData, ParticleContainer, DepositionType)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1     );
  CH_assert(!a_particles.isCellSorted());    

  this->depositParticles<ItoParticle, &ItoParticle::totalEnergy>(a_phi, a_particles, a_deposition);
}

void ItoSolver::computeAverageMobility(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles) const {
  CH_TIME("ItoSolver::computeAverageMobility(EBAMRCellData, ParticleContainer)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeAverageMobility(EBAMRCellData, ParticleContainer)" << endl;
  }

  CH_assert( a_phi[0]->nComp() == 1    );
  CH_assert(!a_particles.isCellSorted());

  // TLDR: We compute the average mobility as the average mobility of all particles (mass-weighted). We do this by depositing the particle conductivity and then
  //       dividing by the mass. 
  DataOps::setValue(a_phi,     0.0);
  DataOps::setValue(m_scratch, 0.0);

  // Need scratch storage to deposit into (can't use m_scratch)
  EBAMRCellData mass;
  m_amr->allocate(mass, m_realm, m_phase, m_nComp);
  
  this->depositParticles<ItoParticle, &ItoParticle::conductivity>(a_phi, a_particles, m_deposition);  // Deposit mass*mu
  this->depositParticles<ItoParticle, &ItoParticle::mass>        (mass,  a_particles, m_deposition);  // Deposit mass

  // Make averageMobility = mass*mu/mass. If there is no mass then set the value to zero. 
  constexpr Real zero = 0.0;
  
  DataOps::divideFallback(a_phi, mass, zero);
}

void ItoSolver::computeAverageDiffusion(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles) const {
  CH_TIME("ItoSolver::computeAverageDiffusion(EBAMRCellData, ParticleContainer)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeAverageDiffusion(EBAMRCellData, ParticleContainer)" << endl;
  }

  CH_assert( a_phi[0]->nComp() == 1    );
  CH_assert(!a_particles.isCellSorted());

  // TLDR: We compute the average mobility as the average mobility of all particles (mass-weighted). We do this by depositing the particle conductivity and then
  //       dividing by the mass. 
  DataOps::setValue(a_phi,     0.0);
  DataOps::setValue(m_scratch, 0.0);

  // Need scratch storage to deposit into (can't use m_scratch)
  EBAMRCellData mass;
  m_amr->allocate(mass, m_realm, m_phase, m_nComp);
  
  this->depositParticles<ItoParticle, &ItoParticle::diffusivity>(a_phi, a_particles, m_deposition);  // Deposit mass*D
  this->depositParticles<ItoParticle, &ItoParticle::mass>       (mass,  a_particles, m_deposition);  // Deposit mass


  // Make averageMobility = mass*mu/mass. If there is no mass then set the value to zero. 
  constexpr Real zero = 0.0;
  
  DataOps::divideFallback(a_phi, mass, zero);  
}

void ItoSolver::computeAverageEnergy(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles) const {
  CH_TIME("ItoSolver::computeAverageEnergy(EBAMRCellData, ParticleContainer)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeAverageEnergy(EBAMRCellData, ParticleContainer)" << endl;
  }

  CH_assert( a_phi[0]->nComp() == 1    );
  CH_assert(!a_particles.isCellSorted());

  // TLDR: We compute the average mobility as the average mobility of all particles (mass-weighted). We do this by depositing the particle conductivity and then
  //       dividing by the mass. 
  DataOps::setValue(a_phi,     0.0);
  DataOps::setValue(m_scratch, 0.0);

  // Need scratch storage to deposit into (can't use m_scratch)
  EBAMRCellData mass;
  m_amr->allocate(mass, m_realm, m_phase, m_nComp);  
  
  this->depositParticles<ItoParticle, &ItoParticle::totalEnergy>(a_phi, a_particles, m_deposition);  // Deposit mass*energy
  this->depositParticles<ItoParticle, &ItoParticle::mass>       (mass,  a_particles, m_deposition);  // Deposit mass

  // Make averageMobility = mass*mu/mass. If there is no mass then set the value to zero. 
  constexpr Real zero = 0.0;
  
  DataOps::divideFallback(a_phi, mass, zero);    
}

void ItoSolver::depositParticles() {
  CH_TIME("ItoSolver::depositParticles");
  if(m_verbosity > 5) {
    pout() << m_name + "::depositParticles" << endl;
  }

  this->depositParticles(WhichContainer::Bulk);
}

void ItoSolver::depositParticles(const WhichContainer a_container) {
  CH_TIME("ItoSolver::depositParticles(container)");
  if(m_verbosity > 5) {
    pout() << m_name + "::depositParticles(container)" << endl;
  }

  this->depositParticles<ItoParticle, &ItoParticle::mass>(m_phi, m_particleContainers.at(a_container), m_deposition);
}

void ItoSolver::redistributeAMR(EBAMRCellData& a_phi) const {
  CH_TIME("ItoSolver::redistributeAMR");
  if(m_verbosity > 5) {
    pout() << m_name + "::redistributeAMR" << endl;
  }

  if(m_useRedistribution) {
    this->depositNonConservative(m_depositionNC, a_phi);              // Compute m_depositionNC = sum(kappa*Wc)/sum(kappa)
    this->depositHybrid(a_phi, m_massDiff, m_depositionNC);           // Compute hybrid deposition, including mass differnce
    this->incrementRedist(m_massDiff);                                  // Increment level redistribution register

    // Do the redistribution magic
    const bool ebcf = m_amr->getEbCf();
    if(ebcf) { // Mucho stuff to do here...
      this->coarseFineIncrement(m_massDiff);       // Compute C2F, F2C, and C2C mass transfers
      this->levelRedist(a_phi);           // Level redistribution. Weights is a dummy parameter
      this->coarseFineRedistribution(a_phi);     // Do the coarse-fine redistribution
    }
    else{ // Very simple, redistribute this level.
      this->levelRedist(a_phi);
    }
  }
}

void ItoSolver::depositNonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC) const {
  CH_TIME("ItoSolver::depositNonConservative");
  if(m_verbosity > 5) {
    pout() << m_name + "::depositNonConservative" << endl;
  }

  const std::string cur_Realm = a_depositionNC.getRealm();

  if(m_blendConservation) {
    const IrregAmrStencil<NonConservativeDivergenceStencil>& stencils = m_amr->getNonConservativeDivergenceStencils(cur_Realm, m_phase);
    stencils.apply(a_depositionNC, a_depositionKappaC);
  }
  else{
    DataOps::setValue(a_depositionNC, 0.0);
  }
}

void ItoSolver::depositHybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_massDifference, const EBAMRIVData& a_depositionNC) const {
  CH_TIME("ItoSolver::depositHybrid");
  if(m_verbosity > 5) {
    pout() << m_name + "::depositHybrid" << endl;
  }

  const std::string cur_Realm = a_depositionH.getRealm();

  const int comp  = 0;
  const int ncomp = 1;


  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(cur_Realm)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(cur_Realm, m_phase)[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      EBCellFAB& divH               = (*a_depositionH[lvl])[dit()];  // On input, this contains kappa*depositionWeights
      BaseIVFAB<Real>& deltaM       = (*a_massDifference[lvl])[dit()];
      const BaseIVFAB<Real>& divNC  = (*a_depositionNC[lvl])[dit()]; 
      const EBISBox& ebisbox        = ebisl[dit()];

      VoFIterator& vofit = (*m_amr->getVofIterator(cur_Realm, m_phase)[lvl])[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit) {
	const VolIndex& vof = vofit();
	const Real kappa    = ebisbox.volFrac(vof);
	const Real dc       = divH(vof, comp);
	const Real dnc      = divNC(vof, comp);

	// Note that if dc - kappa*dnc can be negative, i.e. we may end up STEALING mass
	// from other cells. This is why there is a flag m_blendConservation which always
	// gives positive definite results. 
	divH(vof, comp)   = dc + (1.0-kappa)*dnc;        // On output, contains hybrid divergence
	deltaM(vof, comp) = (1-kappa)*(dc - kappa*dnc);  // Remember, dc already scaled by kappa.
      }
    }
  }
}


void ItoSolver::incrementRedist(const EBAMRIVData& a_massDifference) const {
  CH_TIME("ItoSolver::incrementRedist");
  if(m_verbosity > 5) {
    pout() << m_name + "::incrementRedist" << endl;
  }

  const std::string cur_Realm = a_massDifference.getRealm();

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(cur_Realm)[lvl];
    
    EBLevelRedist& level_redist = *(m_amr->getLevelRedist(cur_Realm, m_phase)[lvl]);
    level_redist.setToZero();

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      level_redist.increment((*a_massDifference[lvl])[dit()], dit(), interv);
    }
  }
}

void ItoSolver::levelRedist(EBAMRCellData& a_phi) const {
  CH_TIME("ItoSolver::levelRedist");
  if(m_verbosity > 5) {
    pout() << m_name + "::levelRedist" << endl;
  }

  const std::string cur_Realm = a_phi.getRealm();

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++) {
    EBLevelRedist& level_redist = *(m_amr->getLevelRedist(cur_Realm, m_phase)[lvl]);
    level_redist.redistribute(*a_phi[lvl], interv);
    level_redist.setToZero();
  }
}

void ItoSolver::coarseFineIncrement(const EBAMRIVData& a_massDifference) const {
  CH_TIME("ItoSolver::coarseFineIncrement");
  if(m_verbosity > 5) {
    pout() << m_name + "::coarseFineIncrement" << endl;
  }

  const std::string cur_Realm = a_massDifference.getRealm();

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(0,0);

  for (int lvl = 0; lvl <= finest_level; lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(cur_Realm)[lvl];

    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->getFineToCoarRedist(cur_Realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->getCoarToFineRedist(cur_Realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->getCoarToCoarRedist(cur_Realm, m_phase)[lvl];

    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < 0;

    if(has_coar) {
      fine2coar_redist->setToZero();

    }
    if(has_fine) {
      coar2fine_redist->setToZero();
      coar2coar_redist->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      if(has_coar) {
	fine2coar_redist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
      }

      if(has_fine) {
	coar2fine_redist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
	coar2coar_redist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
      }
    }
  }
}

void ItoSolver::coarseFineRedistribution(EBAMRCellData& a_phi) const {
  CH_TIME("ItoSolver::coarseFineRedistribution");
  if(m_verbosity > 5) {
    pout() << m_name + "::coarseFineRedistribution" << endl;
  }

  const std::string cur_Realm = a_phi.getRealm();

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);


  for (int lvl = 0; lvl <= finest_level; lvl++) {
    const Real dx       = m_amr->getDx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->getCoarToFineRedist(cur_Realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->getCoarToCoarRedist(cur_Realm, m_phase)[lvl];
    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->getFineToCoarRedist(cur_Realm, m_phase)[lvl];
    if(has_coar) {
      fine2coar_redist->redistribute(*a_phi[lvl-1], interv);
      fine2coar_redist->setToZero();
    }

    if(has_fine) {
      coar2fine_redist->redistribute(*a_phi[lvl+1], interv);
      coar2coar_redist->redistribute(*a_phi[lvl],   interv);

      coar2fine_redist->setToZero();
      coar2coar_redist->setToZero();
    }
  }
}



bool ItoSolver::isMobile() const{
  CH_TIME("ItoSolver::isMobile");
  
  return m_isMobile;
}
  

bool ItoSolver::isDiffusive() const{
  CH_TIME("ItoSolver::isDiffusive");
  
  return m_isDiffusive;
}

void ItoSolver::preRegrid(const int a_lbase, const int a_oldFinestLevel) {
  CH_TIME("ItoSolver::preRegrid");
  if(m_verbosity > 5) {
    pout() << m_name + "::preRegrid" << endl;
  }

  CH_assert(a_lbase >= 0);

  // TLDR: Do a pre-regrid operation for all particle containers owned by the ItoSolver. 

  for (auto& container : m_particleContainers) {
    ParticleContainer<ItoParticle>& particles = container.second;

    particles.preRegrid(a_lbase);
  }
}

ParticleContainer<ItoParticle>& ItoSolver::getParticles(const WhichContainer a_container) {
  CH_TIME("ItoSolver::getParticles");
  if(m_verbosity > 5) {
    pout() << m_name + "::getParticles" << endl;
  }
  
  return m_particleContainers.at(a_container);
}

const ParticleContainer<ItoParticle>& ItoSolver::getParticles(const WhichContainer a_container) const {
  CH_TIME("ItoSolver::getParticles");
  if(m_verbosity > 5) {
    pout() << m_name + "::getParticles" << endl;
  }
  
  return m_particleContainers.at(a_container);
}

EBAMRCellData& ItoSolver::getPhi() {
  CH_TIME("ItoSolver::getPhi");
  if(m_verbosity > 5) {
    pout() << m_name + "::getPhi" << endl;
  }

  return m_phi;
}

EBAMRCellData& ItoSolver::getVelocityFunction() {
  CH_TIME("ItoSolver::getVelocityFunction");
  if(m_verbosity > 5) {
    pout() << m_name + "::getVelocityFunction" << endl;
  }

  return m_velocityFunction;
}

const EBAMRCellData& ItoSolver::getVelocityFunction() const {
  CH_TIME("ItoSolver::getVelocityFunction");
  if(m_verbosity > 5) {
    pout() << m_name + "::getVelocityFunction" << endl;
  }

  return m_velocityFunction;
}

EBAMRCellData& ItoSolver::getDiffusionFunction() {
  CH_TIME("ItoSolver::getDiffusionFunction");
  if(m_verbosity > 5) {
    pout() << m_name + "::getDiffusionFunction" << endl;
  }

  return m_diffusionFunction;
}

const EBAMRCellData& ItoSolver::getDiffusionFunction() const {
  CH_TIME("ItoSolver::getDiffusionFunction");
  if(m_verbosity > 5) {
    pout() << m_name + "::getDiffusionFunction" << endl;
  }

  return m_diffusionFunction;
}

EBAMRCellData& ItoSolver::getScratch() {
  CH_TIME("ItoSolver::getScratch");
  if(m_verbosity > 5) {
    pout() << m_name + "::getScratch" << endl;
  }

  return m_scratch;
}

const EBAMRCellData& ItoSolver::getScratch() const {
  CH_TIME("ItoSolver::getScratch");
  if(m_verbosity > 5) {
    pout() << m_name + "::getScratch" << endl;
  }

  return m_scratch;
}

EBAMRCellData& ItoSolver::getMobilityFunction() {
  CH_TIME("ItoSolver::getMobilityFunction");
  if(m_verbosity > 5) {
    pout() << m_name + "::getMobilityFunction" << endl;
  }

  return m_mobilityFunction;
}

const EBAMRCellData& ItoSolver::getMobilityFunction() const {
  CH_TIME("ItoSolver::getMobilityFunction");
  if(m_verbosity > 5) {
    pout() << m_name + "::getMobilityFunction" << endl;
  }

  return m_mobilityFunction;
}

void ItoSolver::setDiffusionFunction(const Real a_diffusionCoefficient) {
  CH_TIME("ItoSolver::setDiffusionFunction");
  if(m_verbosity > 5) {
    pout() << m_name + "::setDiffusionFunction" << endl;
  }

  DataOps::setValue(m_diffusionFunction, a_diffusionCoefficient);
}

void ItoSolver::setVelocityFunction(const RealVect a_velocity) {
  CH_TIME("ItoSolver::setVelocityFunction");
  if(m_verbosity > 5) {
    pout() << m_name + "::setVelocityFunction" << endl;
  }

  for (int dir = 0; dir < SpaceDim; dir++) {
    DataOps::setValue(m_velocityFunction, a_velocity[dir], dir);
  }
}

void ItoSolver::setParticleMobility(const Real a_mobility) {
  CH_TIME("ItoSolver::setParticleMobility");
  if(m_verbosity > 5) {
    pout() << m_name + "::setParticleMobility" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::Bulk);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<ItoParticle>& particlesDit = particles[lvl][dit()].listItems();

      for (ListIterator<ItoParticle> lit(particlesDit); lit.ok(); ++lit) {
	lit().mobility() = a_mobility;
      }
    }
  }
}

void ItoSolver::interpolateVelocities() {
  CH_TIME("ItoSolver::interpolateVelocities");
  if(m_verbosity > 5) {
    pout() << m_name + "::interpolateVelocities" << endl;
  }

  if(m_isMobile) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
	this->interpolateVelocities(lvl, dit());
      }
    }
  }
}

void ItoSolver::interpolateVelocities(const int a_lvl, const DataIndex& a_dit) {
  CH_TIME("ItoSolver::interpolateVelocities");
  if(m_verbosity > 5) {
    pout() << m_name + "::interpolateVelocities" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  if(m_isMobile) {
    const EBCellFAB& velo_func = (*m_velocityFunction[a_lvl])[a_dit];
    const EBISBox& ebisbox     = velo_func.getEBISBox();
    const RealVect dx          = m_amr->getDx()[a_lvl]*RealVect::Unit;
    const RealVect origin      = m_amr->getProbLo();
    const Box box              = m_amr->getGrids(m_realm)[a_lvl][a_dit];

    List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();

    // This interpolates the velocity function on to the particle velocities
    EBParticleMesh meshInterp(box, ebisbox, dx, origin);
    meshInterp.interpolate<ItoParticle, &ItoParticle::velocity>(particleList, velo_func, m_deposition, m_forceIrregInterpolationNGP);

    // Go through the particles and set their velocities to velo_func*mobility
    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
      ItoParticle& p = lit();
      p.velocity() *= p.mobility();
    }
  }
}

void ItoSolver::interpolateMobilities() {
  CH_TIME("ItoSolver::interpolateMobilities()");
  if(m_verbosity > 5) {
    pout() << m_name + "::interpolateMobilities()" << endl;
  }

  if(m_isMobile) {

    switch(m_mobilityInterpolation){
    case WhichMobilityInterpolation::Velocity:
      {
	DataOps::vectorLength(m_scratch, m_velocityFunction); // Compute |E| (or whatever other function you've decided to provide).
	m_amr->averageDown(m_scratch, m_realm, m_phase);
	m_amr->interpGhost(m_scratch, m_realm, m_phase);
      }
    default: // Do nothing
      break;
    }

    // Call the level version.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
	this->interpolateMobilities(lvl, dit());
      }
    }
  }
}

void ItoSolver::interpolateMobilities(const int a_lvl, const DataIndex& a_dit) {
  CH_TIME("ItoSolver::interpolateMobilities(lvl, dit)");
  if(m_verbosity > 5) {
    pout() << m_name + "::interpolateMobilities(lvl, dit)" << endl;
  }

  CH_assert(m_isMobile);

  switch(m_mobilityInterpolation){
  case WhichMobilityInterpolation::Direct:
    this->interpolateMobilitiesMu(a_lvl, a_dit);
    break;
  case WhichMobilityInterpolation::Velocity:
    this->interpolateMobilitiesVel(a_lvl, a_dit);
    break;
  default:
    MayDay::Error("ItoSolver::interpolateMobilities(int, DataIndex) - logic bust");
  }
}

void ItoSolver::interpolateMobilitiesMu(const int a_lvl, const DataIndex& a_dit) {
  CH_TIME("ItoSolver::interpolateMobilitiesMu");
  if(m_verbosity > 5) {
    pout() << m_name + "::interpolateMobilitiesMu" << endl;
  }

  CH_assert(m_isMobile);
  CH_assert(m_mobilityInterpolation == WhichMobilityInterpolation::Direct);

  // TLDR: This will compute the particle mobility by interpolating a scalar mobility field (stored on the mesh) to the particle positions.

  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  const EBCellFAB& mobilityFunction = (*m_mobilityFunction[a_lvl])[a_dit];
  const EBISBox&   ebisbox          = mobilityFunction.getEBISBox();
  const RealVect   dx               = m_amr->getDx()[a_lvl]*RealVect::Unit;
  const RealVect   probLo           = m_amr->getProbLo();
  const Box        box              = m_amr->getGrids(m_realm)[a_lvl][a_dit];

  List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();

  // Interpolate onto the mobility field
  EBParticleMesh meshInterp(box, ebisbox, dx, probLo);
  meshInterp.interpolate<ItoParticle, &ItoParticle::mobility>(particleList, mobilityFunction, m_deposition, m_forceIrregInterpolationNGP);
}

void ItoSolver::interpolateMobilitiesVel(const int a_lvl, const DataIndex& a_dit) {
  CH_TIME("ItoSolver::interpolateMobilitiesVel");
  if(m_verbosity > 5) {
    pout() << m_name + "::interpolateMobilitiesVel" << endl;
  }

  CH_assert(m_isMobile);
  CH_assert(m_mobilityInterpolation == WhichMobilityInterpolation::Velocity);

  // TLDR: This function computes the particle mobilities by interpolating mu*V to the particle position and then setting
  //       the mobility as mu = [mu*V(Xp)]/V(Xp). We happen to know that |V| is already stored in m_scratch. 

  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  const EBCellFAB& mobilityFunction = (*m_mobilityFunction[a_lvl])[a_dit];
  const EBISBox&   ebisbox          = mobilityFunction.getEBISBox();
  const RealVect   dx               = m_amr->getDx()[a_lvl]*RealVect::Unit;
  const RealVect   probLo           = m_amr->getProbLo();
  const Box        box              = m_amr->getGrids(m_realm)[a_lvl][a_dit];

  EBCellFAB& scratch = (*m_scratch[a_lvl])[a_dit];

  List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();
  EBParticleMesh meshInterp(box, ebisbox, dx, probLo);
    
  // First, interpolate |V| to the particle position, it will be stored on m_tmp. 
  meshInterp.interpolate<ItoParticle, &ItoParticle::tmp>(particleList, scratch, m_deposition, m_forceIrregInterpolationNGP);    

  // Secondly, let m_scratch hold mu*|V| and interpolate that to the particle mobility field. 
  scratch *= mobilityFunction;
  meshInterp.interpolate<ItoParticle, &ItoParticle::mobility>(particleList, scratch, m_deposition, m_forceIrregInterpolationNGP);

  // We now have ItoParticle::tmp = |V(Xp)| and ItoParticle::mobility = |mu*V|(Xp). Now let the mobility be mu(Xp) = |mu*V|(Xp)/|V|(Xp)
  for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
    ItoParticle& p = lit();

    p.mobility() *= 1./p.tmp();
  }
}

void ItoSolver::updateMobilities() {
  CH_TIME("ItoSolver::updateMobilities");
  if(m_verbosity > 5) {
    pout() << m_name + "::updateMobilities" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      this->updateMobilities(lvl, dit());
    }
  }
}

void ItoSolver::updateMobilities(const int a_level, const DataIndex a_dit) {
  CH_TIME("ItoSolver::updateMobilities(lvl, dit)");
  if(m_verbosity > 5) {
    pout() << m_name + "::updateMobilities(lvl, dit)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  if(m_isMobile) {
    List<ItoParticle>& particleList = particles[a_level][a_dit].listItems();

    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
      ItoParticle& p = lit();
      
      p.mobility() = m_species->mobility(p.energy());
    }
  }
}

void ItoSolver::interpolateDiffusion() {
  CH_TIME("ItoSolver::interpolateDiffusion");
  if(m_verbosity > 5) {
    pout() << m_name + "::interpolateDiffusion" << endl;
  }

  if(m_isDiffusive) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
	this->interpolateDiffusion(lvl, dit());
      }
    }
  }
}

void ItoSolver::interpolateDiffusion(const int a_lvl, const DataIndex& a_dit) {
  CH_TIME("ItoSolver::interpolateDiffusion");
  if(m_verbosity > 5) {
    pout() << m_name + "::interpolateDiffusion" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  if(m_isDiffusive) {
    const EBCellFAB& dco_cell  = (*m_diffusionFunction[a_lvl])[a_dit];
    const EBISBox& ebisbox     = dco_cell.getEBISBox();
    const RealVect dx          = m_amr->getDx()[a_lvl]*RealVect::Unit;
    const RealVect origin      = m_amr->getProbLo();
    const Box box              = m_amr->getGrids(m_realm)[a_lvl][a_dit];

    List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();

    EBParticleMesh meshInterp(box, ebisbox,dx, origin);
    meshInterp.interpolate<ItoParticle, &ItoParticle::diffusion>(particleList, dco_cell, m_deposition, m_forceIrregInterpolationNGP);    
  }
}

void ItoSolver::updateDiffusion() {
  CH_TIME("ItoSolver::updateDiffusion");
  if(m_verbosity > 5) {
    pout() << m_name + "::updateDiffusion" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      this->updateDiffusion(lvl, dit());
    }
  }
}

void ItoSolver::updateDiffusion(const int a_level, const DataIndex a_dit) {
  CH_TIME("ItoSolver::updateDiffusion(lvl, dit)");
  if(m_verbosity > 5) {
    pout() << m_name + "::updateDiffusion(lvl, dit)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::Bulk);

  if(m_isDiffusive) {
    List<ItoParticle>& particleList = particles[a_level][a_dit].listItems();

    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
      ItoParticle& p = lit();
      
      p.mobility() = m_species->diffusion(p.energy());
    }
  }
}



Real ItoSolver::computeDt() const {
  CH_TIME("ItoSolver::computeDt(allAMRlevels)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeDt(allAMRlevels)" << endl;
  }

  Real dt = 1.E99;
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real levelDt = this->computeDt(lvl);
    dt = Min(dt, levelDt);
  }

  return dt;
}

Real ItoSolver::computeDt(const int a_lvl) const{
  CH_TIME("ItoSolver::computeDt(lvl)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeDt(lvl)" << endl;
  }

  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const Real dx = m_amr->getDx()[a_lvl];
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const Real boxDt = this->computeDt(a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS) {
    MayDay::Error("ItoSolver::computeDt(lvl) - communication error on norm");
  }
  dt = tmp;
#endif  

  return dt;
}
  

Real ItoSolver::computeDt(const int a_lvl, const DataIndex a_dit, const Real a_dx) const{
  CH_TIME("ItoSolver::computeDt(lvl, dit, dx)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeDt(lvl, dit, dx)" << endl;
  }

  Real dt = 1.E99;

  const ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::Bulk);

  const List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();
  
  ListIterator<ItoParticle> lit(particleList);
  if(m_isMobile && !m_isDiffusive) {
    for (lit.rewind(); lit.ok(); ++lit) {
      const ItoParticle& p = particleList[lit];
      const RealVect& v = p.velocity();

      const int maxDir  = v.maxDir(true);
      const Real vMax   = Abs(v[maxDir]);
      const Real thisDt = (vMax > 0.0) ? a_dx/vMax : 1.E99;

      dt = Min(dt, thisDt);
    }
  }
  else if(!m_isMobile && m_isDiffusive) {
    for (lit.rewind(); lit.ok(); ++lit) {
      const ItoParticle& p = particleList[lit];

      const Real D      = p.diffusion();
      const Real thisDt = (D > 0.0) ? a_dx*a_dx/(2.0*D) : 1.E99;

      dt = Min(dt, thisDt);
    }
  }
  else if(m_isMobile && m_isDiffusive) {
    for (lit.rewind(); lit.ok(); ++lit) {
      const ItoParticle& p = particleList[lit];
      
      const RealVect& v = p.velocity();
      const int maxDir = v.maxDir(true);
      const Real vMax  = Abs(v[maxDir]);
      const Real D     = p.diffusion();

      const Real dtAdvect = (vMax > 0.0) ? a_dx/vMax         : 1.E99;
      const Real dtDiffus = (D > 0.0)    ? a_dx*a_dx/(2.0*D) : 1.E99;

      const Real thisDt = 1./(1./dtAdvect + 1./dtDiffus);
      
      dt = Min(dt, thisDt);
    }
  }

  return dt;
}

Real ItoSolver::computeMinDt(const Real a_maxCellsToMove) const {
  CH_TIME("ItoSolver::computeMinDt(allAMRlevels, maxCellsToMove)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeMinDt(allAMRlevels, maxCellsToMove)" << endl;
  }

  Real dt = 1.E99;
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real levelDt = this->computeMinDt(a_maxCellsToMove, lvl);
    dt = Min(dt, levelDt);
  }

  return dt;
}

Real ItoSolver::computeMinDt(const Real a_maxCellsToMove, const int a_lvl) const{
  CH_TIME("ItoSolver::computeMinDt(maxCellsToMove, lvl)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeMinDt(maxCellsToMove, lvl)" << endl;
  }

  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const Real dx = m_amr->getDx()[a_lvl];
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const Real boxDt = this->computeMinDt(a_maxCellsToMove, a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS) {
    MayDay::Error("ItoSolver::computeDt(lvl) - communication error on norm");
  }
  dt = tmp;
#endif  

  return dt;
}

Real ItoSolver::computeMinDt(const Real a_maxCellsToMove, const int a_lvl, const DataIndex a_dit, const Real a_dx) const{
  CH_TIME("ItoSolver::computeMinDt(maxCellsToMove, lvl, dit, dx)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeMinDt(maxCellsToMove, lvl, dit, dx)" << endl;
  }

  Real dt = 1.E99;

  const ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::Bulk);

  const Real dMax  = a_maxCellsToMove*a_dx;
  const Real dMax2 = dMax*dMax;
  const Real W0    = m_normalDistributionTruncation;
  const Real W02   = m_normalDistributionTruncation*m_normalDistributionTruncation;

  const List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();
  
  ListIterator<ItoParticle> lit(particleList);
  if(m_isMobile && !m_isDiffusive) {
    for (lit.rewind(); lit; ++lit) {
      const ItoParticle& p = particleList[lit];
      const RealVect& v = p.velocity();

      const int maxDir = v.maxDir(true);
      const Real thisDt = dMax/Abs(v[maxDir]);

      dt = Min(dt, thisDt);
    }
  }
  else if(!m_isMobile && m_isDiffusive) {
    for (lit.rewind(); lit; ++lit) {
      const ItoParticle& p = particleList[lit];
      
      const Real thisDt = dMax2/(2.0*p.diffusion()*W02);
      dt = Min(dt, thisDt);
    }
  }
  else if(m_isMobile && m_isDiffusive) {
    for (lit.rewind(); lit; ++lit) {
      const ItoParticle& p = particleList[lit];
      
      const RealVect& v = p.velocity();
      const int maxDir = v.maxDir(true);
      const Real vMax  = Abs(v[maxDir]);
      const Real D     = p.diffusion();

      
      if(vMax > 0.0) {
	const Real a = vMax;
	const Real b = W0*sqrt(2.0*D);
	const Real c = dMax;
	
	const Real A = a*a;
	const Real B = -(b*b + 2*a*c);
	const Real C = c*c;

	const Real thisDt = (-B - sqrt(B*B - 4.*A*C))/(2.*A);
	dt = Min(dt, thisDt);
      }
      else{
	if(D > 0.0) {
	  const Real thisDt = dMax2/(2.0*D*W02);
	  dt = Min(dt, thisDt);
	}
      }
    }
  }

  return dt;

}

Real ItoSolver::computeMinDriftDt(const Real a_maxCellsToMove) const{
  CH_TIME("ItoSolver::computeMinDriftDt(allAMRlevels, maxCellsToMove)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeMinDriftDt(allAMRlevels, maxCellsToMove)" << endl;
  }

  Vector<Real> dt = this->computeDriftDt(a_maxCellsToMove);

  Real minDt = dt[0];
  for (int lvl = 0; lvl < dt.size(); lvl++) {
    minDt = Min(minDt, dt[lvl]);
  }

  return minDt;
}

Vector<Real> ItoSolver::computeDriftDt(const Real a_maxCellsToMove) const {
  CH_TIME("ItoSolver::computeDriftDt(amr, maxCellsToMove)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeDriftDt(amr, maxCellsToMove)" << endl;
  }

  Vector<Real> dt = this->computeDriftDt();
  for (int lvl = 0; lvl < dt.size(); lvl++) {
    dt[lvl] = dt[lvl]*a_maxCellsToMove;
  }

  return dt;

}

Real ItoSolver::computeAdvectiveDt() const {
  CH_TIME("ItoSolver::computeAdvectiveDt");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeAdvectiveDt()" << endl;
  }
  
  Real minDt = 1.E99;
  const Vector<Real> levelDts = this->computeDriftDt();

  for (int lvl = 0; lvl < levelDts.size(); lvl++) {
    minDt = Min(levelDts[lvl], minDt);
  }

  return minDt;
}

Vector<Real> ItoSolver::computeDriftDt() const {
  CH_TIME("ItoSolver::computeDriftDt(amr)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeDriftDt(amr)" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  Vector<Real> dt(1 + finest_level, 1.2345E67);

  for (int lvl = 0; lvl <= finest_level; lvl++) {
    dt[lvl] = this->computeDriftDt(lvl);
  }

  return dt;
}

Real ItoSolver::computeDriftDt(const int a_lvl) const {
  CH_TIME("ItoSolver::computeDriftDt(level)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeDriftDt(level)" << endl;
  }

  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const RealVect dx = m_amr->getDx()[a_lvl]*RealVect::Unit;
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const Real boxDt = this->computeDriftDt(a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS) {
    MayDay::Error("ItoSolver::compute_drift_level(lvl) - communication error on norm");
  }
  dt = tmp;
#endif  

  return dt;
}

Real ItoSolver::computeDriftDt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const{
  CH_TIME("ItoSolver::computeDriftDt(level, dataindex, dx)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeDriftDt(level, dataindex, dx)" << endl;
  }

  const ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  constexpr Real safety = 1.E-10;

  const List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();

  Real dt = 1.E99;

  if(m_isMobile) {
    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
      const ItoParticle& p = particleList[lit];
      const RealVect& v     = p.velocity();
      const int maxDir      = v.maxDir(true);
      const Real thisDt     = a_dx[maxDir]/(safety + Abs(v[maxDir]));

      dt = Min(dt, thisDt);
    }
  }

  return dt;
}

Real ItoSolver::computeMinDiffusionDt(const Real a_maxCellsToHop) const{
  CH_TIME("ItoSolver::computeMinDiffusionDt(min, maxCellsToHop)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeMinDiffusionDt(min, maxCellsToHop)" << endl;
  }

  Vector<Real> dt = this->computeDiffusionDt(a_maxCellsToHop);
  Real minDt = dt[0];
  for (int lvl = 0; lvl < dt.size(); lvl++) {
    minDt = Min(minDt, dt[lvl]);
  }

  return minDt;
}

Vector<Real> ItoSolver::computeDiffusionDt(const Real a_maxCellsToHop) const{
  CH_TIME("ItoSolver::computeMinDiffusionDt(maxCellsToHop)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeMinDiffusionDt(maxCellsToHop)" << endl;
  }

  const Real factor  = a_maxCellsToHop/m_normalDistributionTruncation;
  const Real factor2 = factor*factor;
  
  Vector<Real> dt = this->computeDiffusionDt();
  for (int lvl = 0; lvl < dt.size(); lvl++) {
    dt[lvl] = dt[lvl]*factor2;
  }

  return dt;
}

Real ItoSolver::computeDiffusiveDt() const {
  CH_TIME("ItoSolver::computeDiffusiveDt");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeDiffusiveDt" << endl;
  }

  Real minDt = 1.E99;
  const Vector<Real> levelDts = this->computeDiffusionDt();

  for (int lvl = 0; lvl < levelDts.size(); lvl++) {
    minDt = Min(levelDts[lvl], minDt);
  }

  return minDt;
}

Vector<Real> ItoSolver::computeDiffusionDt() const{
  CH_TIME("ItoSolver::computeDiffusionDt");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeDiffusionDt" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();
    
  Vector<Real> dt(1 + finest_level, 1.2345E6);

  for (int lvl = 0; lvl <= finest_level; lvl++) {
    dt[lvl] = this->computeDiffusionDt(lvl);
  }

  return dt;
}

Real ItoSolver::computeDiffusionDt(const int a_lvl) const{
  CH_TIME("ItoSolver::computeDiffusionDt(level)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeDiffusionDt(level)" << endl;
  }

  
  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const RealVect dx = m_amr->getDx()[a_lvl]*RealVect::Unit;
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const Real boxDt = this->computeDiffusionDt(a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS) {
    MayDay::Error("ItoSolver::compute_diffusion_level(lvl) - communication error on norm");
  }
  dt = tmp;
#endif

  return dt;
}

Real ItoSolver::computeDiffusionDt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const{
  CH_TIME("ItoSolver::computeDiffusionDt(level, dataindex, dx)");
  if(m_verbosity > 5) {
    pout() << m_name + "::computeDiffusionDt(level, dataindex, dx)" << endl;
  }

  const ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::Bulk);
  
  const List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();

  Real dt = 1.E99;

  if(m_isDiffusive) {
    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
      const ItoParticle& p = particleList[lit];
    
      const Real thisDt = a_dx[0]*a_dx[0]/(2.0*p.diffusion());
    
      dt = Min(dt, thisDt);
    }
  }

  return dt;
}

void ItoSolver::remap() {
  CH_TIME("ItoSolver::remap");
  if(m_verbosity > 5) {
    pout() << m_name + "::remap" << endl;
  }

  this->remap(WhichContainer::Bulk);

}

void ItoSolver::remap(const WhichContainer a_container) {
  CH_TIME("ItoSolver::remap(container)");
  if(m_verbosity > 5) {
    pout() << m_name + "::remap(container)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);
  
  particles.remap();
}

DepositionType ItoSolver::getDeposition() const {
  return m_deposition;
}

phase::which_phase ItoSolver::getPhase() const{
  return m_phase;
}

void ItoSolver::sortParticlesByCell(const WhichContainer a_container) {
  CH_TIME("ItoSolver::sortParticlesByCell(container)");
  if(m_verbosity > 5) {
    pout() << m_name + "::sortParticlesByCell(container)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);

  particles.sortParticlesByCell();
}


void ItoSolver::sortParticlesByPatch(const WhichContainer a_container) {
  CH_TIME("ItoSolver::sortParticlesByPatch(container)");
  if(m_verbosity > 5) {
    pout() << m_name + "::sortParticlesByPatch(container)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);

  particles.sortParticlesByPatch();
}

void ItoSolver::makeSuperparticles(const WhichContainer a_container, const int a_particlesPerPatch) {
  CH_TIME("ItoSolver::makeSuperparticles(int)");
  if(m_verbosity > 5) {
    pout() << m_name + "::makeSuperparticles(int)" << endl;
  }
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->makeSuperparticles(a_container, a_particlesPerPatch, lvl);
  }

}

void ItoSolver::makeSuperparticles(const WhichContainer a_container, const int a_particlesPerPatch, const int a_level) {
  CH_TIME("ItoSolver::makeSuperparticles(int, level)");
  if(m_verbosity > 5) {
    pout() << m_name + "::makeSuperparticles(int, level)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  
  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    this->makeSuperparticles(a_container, a_particlesPerPatch, a_level, dit());
  }
}

void ItoSolver::makeSuperparticles(const WhichContainer a_container, const int a_particlesPerCell, const int a_level, const DataIndex a_dit) {
  CH_TIME("ItoSolver::makeSuperparticles(int, level, patch)");
  if(m_verbosity > 5) {
    pout() << m_name + "::makeSuperparticles(int, level, patch)" << endl;
  }

  constexpr int comp = 0;  

  // This are the particles in the box we're currently looking at.
  ParticleContainer<ItoParticle>& particles     = this->getParticles(a_container);  
  BinFab<ItoParticle>&            cellParticles = particles.getCellParticles(a_level, a_dit);

  // Iterate over particles
  const Box cellBox  = m_amr->getGrids(m_realm)[a_level].get(a_dit);  
  for (BoxIterator bit(cellBox); bit.ok(); ++bit) {
    const IntVect iv = bit();
  
    List<ItoParticle>& particles = cellParticles(iv, comp);

    if(particles.length() > 0) {
      this->mergeBVH(particles, a_particlesPerCell);
    }
  }
}

void ItoSolver::mergeBVH(List<ItoParticle>& a_particles, const int a_particlesPerCell) {
  CH_TIME("ItoSolver::mergeBVH");

#if ITO_DEBUG
  Real mass_before = 0.0;
  Real energy_before = 0.0;
  
  Real mass_after = 0.0;
  Real energy_after = 0.0;
  RealVect center_before = RealVect::Zero;
  RealVect center_after = RealVect::Zero;
#endif
  
  std::vector<PointMass> pointMasses(a_particles.length());
  int i = 0;
  Real mass = 0.0;
  for (ListIterator<ItoParticle> lit(a_particles); lit.ok(); ++lit, i++) {
    const ItoParticle& p = lit();
    pointMasses[i].define(p.position(), p.mass(), p.energy()); // Note, p.energy() is average energy and not total energy. 
    mass += p.mass();
    //    i++;

#if ITO_DEBUG
    mass_before   += p.mass();
    center_before += p.mass()*p.position();
    energy_before += p.mass()*p.energy();
#endif
  }

#if ITO_DEBUG
  center_before *= 1./mass_before;
  if(mass_before   < 1.0)            MayDay::Abort("ItoSolver::mergeBVH - bad initial mass!");
  if(mass_before   != mass_before)   MayDay::Abort("ItoSolver::mergeBVH - initial mass is NaN");
  if(energy_before != energy_before) MayDay::Abort("ItoSolver::mergeBVH - initial energy is NaN");
#endif
  
  // 2. Build the BVH tree and get the leaves of the tree
  const int dir = (m_directionKD < 0) ? m_uniformDistribution0d(m_rng) : m_directionKD;
  m_mergeTree.define(pointMasses, mass);
  m_mergeTree.buildTree(dir, a_particlesPerCell);
  const std::vector<std::shared_ptr<ItoMerge::Node<PointMass> > >& leaves = m_mergeTree.getLeaves();

  // 3. Clear particles in this cell and add new ones.
  a_particles.clear();
  for (int i = 0; i < leaves.size(); i++) {
    PointMass pointMass(leaves[i]->getData());
    ItoParticle p(pointMass.mass(), pointMass.pos(), RealVect::Zero, 0.0, 0.0, pointMass.energy());
    a_particles.add(p);

#if ITO_DEBUG
    mass_after   += p.mass();
    energy_after += p.mass()*p.energy();
    center_after += p.mass()*p.position();
#endif
  }

#if ITO_DEBUG
  center_after *= 1./mass_after;

  constexpr Real BVH_ERR_TOL = 1.E-6;
  const bool break_energy = Abs(energy_before/energy_after - 1.0) > BVH_ERR_TOL;
  const bool break_center = (center_before/center_after - RealVect::Unit).vectorLength() > BVH_ERR_TOL;
  const bool break_mass   = Abs(mass_before/mass_after - 1.0) > BVH_ERR_TOL;
  if(break_mass)   pout() << "ItoSolver::mergeBVH failed. Mass before = "   << mass_before   << "\t Mass after = "   << mass_after << endl;
  if(break_center) pout() << "ItoSolver::mergeBVH failed. center before = " << center_before << "\t center after = " << center_after << endl;
  if(break_energy) pout() << "ItoSolver::mergeBVH failed. Energy before = " << energy_before << "\t Energy after = " << energy_after << endl;
#endif
}

void ItoSolver::clear(const WhichContainer a_container) {
  CH_TIME("ItoSolver::clear(string)");
  if(m_verbosity > 5) {
    pout() << m_name + "::clear(string)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);

  this->clear(particles);
}

void ItoSolver::clear(ParticleContainer<ItoParticle>& a_particles) const {
  CH_TIME("ItoSolver::clear(ParticleContainer)");
  if(m_verbosity > 5) {
    pout() << m_name + "::clear(ParticleContainer)" << endl;
  }

  this->clear(a_particles.getParticles());
}

void ItoSolver::clear(AMRParticles<ItoParticle>& a_particles) const {
  CH_TIME("ItoSolver::clear");
  if(m_verbosity > 5) {
    pout() << m_name + "::clear" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    a_particles[lvl]->clear();
  }
}

#include <CD_NamespaceFooter.H>
