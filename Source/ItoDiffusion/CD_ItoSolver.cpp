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
#include <ParmParse.H>
#include <ParticleIO.H>

// Our includes
#include <CD_SimpleItoParticle.H>
#include <CD_ItoSolver.H>
#include <CD_DataOps.H>
#include <CD_ParticleOps.H>
#include <CD_NamespaceHeader.H>

#define ITO_DEBUG 0

ItoSolver::ItoSolver(){
  m_name      = "ItoSolver";
  m_className = "ItoSolver";

  m_WhichMobilityInterpolation = WhichMobilityInterpolation::mobility;
}

ItoSolver::~ItoSolver(){

}

std::string ItoSolver::getName(){
  return m_name;
}

const std::string ItoSolver::getRealm() const{
  return m_realm;
}

void ItoSolver::setRealm(const std::string a_realm){
  m_realm = a_realm;

  // This is for later, in case we want to move averaging/interpolating onto a different Realm. 
  m_fluid_Realm = m_realm;
}

RefCountedPtr<ItoSpecies>& ItoSolver::getSpecies(){
  return m_species;
}

void ItoSolver::parseOptions(){
  CH_TIME("ItoSolver::parseOptions");
  if(m_verbosity > 5){
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

void ItoSolver::parseRuntimeOptions(){
  CH_TIME("ItoSolver::parseRuntimeOptions");
  if(m_verbosity > 5){
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

void ItoSolver::parseSuperParticles(){
  CH_TIME("ItoSolver::parseSuperParticles");
  if(m_verbosity > 5){
    pout() << m_name + "::parseSuperParticles" << endl;
  }

  // Seed the RNG
  ParmParse pp(m_className.c_str());
  pp.get("kd_direction", m_kd_direction);

  m_kd_direction = min(m_kd_direction, SpaceDim-1);
}

void ItoSolver::parseRng(){
  CH_TIME("ItoSolver::parseRng");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRng" << endl;
  }

  // Seed the RNG
  ParmParse pp(m_className.c_str());
  pp.get("seed",       m_seed_rng);
  pp.get("normal_max", m_normal_max);
  if(m_seed_rng < 0) { // Random seed if input < 0
    m_seed_rng = std::chrono::system_clock::now().time_since_epoch().count();
  }
  
  m_rng     = RAN::mt19937_64(m_seed_rng);
  m_udist01 = RAN::uniform_real_distribution<Real>( 0.0, 1.0);
  m_udist11 = RAN::uniform_real_distribution<Real>(-1.0, 1.0);
  m_gauss01 = RAN::normal_distribution<Real>(0.0, 1.0);
  m_udist0d = RAN::uniform_int_distribution<int>(0, SpaceDim-1);
}

void ItoSolver::parsePlotVariables(){
  CH_TIME("McPhoto::parsePlotVariables");
  if(m_verbosity > 5){
    pout() << m_name + "::parsePlotVariables" << endl;
  }

  m_plotPhi = false;
  m_plotVelocity = false;
  m_plotDiffusionCoefficient = false;
  m_plot_particles        = false;
  m_plot_eb_particles     = false;
  m_plot_domain_particles = false;
  m_plot_source_particles = false;
  m_plot_energy_density   = false;
  m_plot_average_energy   = false;

  ParmParse pp(m_className.c_str());
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++){
    if(     str[i] == "phi")            m_plotPhi              = true;
    else if(str[i] == "vel")            m_plotVelocity              = true;
    else if(str[i] == "dco")            m_plotDiffusionCoefficient              = true;
    else if(str[i] == "part")           m_plot_particles        = true;
    else if(str[i] == "eb_part")        m_plot_eb_particles     = true;
    else if(str[i] == "dom_part")       m_plot_domain_particles = true;
    else if(str[i] == "src_part")       m_plot_source_particles = true;
    else if(str[i] == "energy_density") m_plot_energy_density   = true;
    else if(str[i] == "average_energy") m_plot_average_energy   = true;
  }
}

void ItoSolver::parseDeposition(){
  CH_TIME("ItoSolver::parseRng");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRng" << endl;
  }

  ParmParse pp(m_className.c_str());
  std::string str;

  // Deposition for particle-mesh operations
  pp.get("deposition", str);
  if(str == "ngp"){
    m_deposition = DepositionType::NGP;
  }
  else if(str == "cic"){
    m_deposition = DepositionType::CIC;
  }
  else if(str == "tsc"){
    m_deposition = DepositionType::TSC;
  }
  else if(str == "w4"){
    m_deposition = DepositionType::W4;
  }
  else{
    MayDay::Abort("ItoSolver::parseDeposition - unknown interpolant requested");
  }

  // Deposition for plotting only
  pp.get("plot_deposition", str);

  if(str == "ngp"){
    m_plot_deposition = DepositionType::NGP;
  }
  else if(str == "cic"){
    m_plot_deposition = DepositionType::CIC;
  }
  else if(str == "tsc"){
    m_plot_deposition = DepositionType::TSC;
  }
  else if(str == "w4"){
    m_plot_deposition = DepositionType::W4;
  }
  else{
    MayDay::Abort("ItoSolver::parseDeposition - unknown interpolant requested");
  }

  // Mobility interpolation.
  pp.get("WhichMobilityInterpolation",str);

  if(str == "mobility"){
    m_WhichMobilityInterpolation = WhichMobilityInterpolation::mobility;
  }
  else if(str == "velocity"){
    m_WhichMobilityInterpolation = WhichMobilityInterpolation::velocity;
  }
  else{
    MayDay::Abort("ItoSolver::parseDeposition - unknown interpolation method for mobility");
  }

  pp.get("irr_ngp_deposition", m_irreg_ngp_deposition);
  pp.get("irr_ngp_interp",     m_irreg_ngp_interpolation);
}

void ItoSolver::parseBisectStep(){
  CH_TIME("ItoSolver::parseBisectStep");
  if(m_verbosity > 5){
    pout() << m_name + "::parseBisectStep" << endl;
  }

  ParmParse pp(m_className.c_str());
  pp.get("bisect_step", m_bisect_step);
}

void ItoSolver::parsePvrBuffer(){
  CH_TIME("ItoSolver::parsePvrBuffer");
  if(m_verbosity > 5){
    pout() << m_name + "::parsePvrBuffer" << endl;
  }

  ParmParse pp(m_className.c_str());
  pp.get("pvr_buffer",  m_pvr_buffer);
  pp.get("halo_buffer", m_halo_buffer);

  std::string str;
  pp.get("halo_deposition", str);
  if(str == "ngp"){
    m_ngp_halo = true;
  }
  else if(str == "native"){
    m_ngp_halo = false;
  }
  else{
    MayDay::Abort("ItoSolver::parsePvrBuffer - unknown argument to 'halo_deposition'");
  }
}

void ItoSolver::parseDiffusionHop(){
  CH_TIME("ItoSolver::parseDiffusionHop");
  if(m_verbosity > 5){
    pout() << m_name + "::parseDiffusionHop" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("max_diffusion_hop", m_max_diffusion_hop);
}

void ItoSolver::parseRedistribution(){
  CH_TIME("ItoSolver::parseRedistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRedistribution" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("redistribute", m_redistribute);
}

void ItoSolver::parseDivergenceComputation(){
  CH_TIME("ItoSolver::parseDivergenceComputation");
  if(m_verbosity > 5){
    pout() << m_name + "::parseDivergenceComputation" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("blend_conservation", m_blendConservation);
}

void ItoSolver::parseCheckpointing(){
  CH_TIME("ItoSolver::parseCheckpointing");
  if(m_verbosity > 5){
    pout() << m_name + "::parseCheckpointing" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("checkpointing", str);
  pp.get("ppc_restart", m_ppc_restart);
  if(str == "particles"){
    m_checkpointing = WhichCheckpoint::particles;
  }
  else if(str == "numbers"){
    m_checkpointing = WhichCheckpoint::numbers;
  }
  else{
    MayDay::Abort("ItoSolver::parseCheckpointing - unknown checkpointing method requested");
  }
}

Vector<std::string> ItoSolver::getPlotVariableNames() const {
  CH_TIME("ItoSolver::getPlotVariableNames");
  if(m_verbosity > 5){
    pout() << m_name + "::getPlotVariableNames" << endl;
  }

  Vector<std::string> names(0);
  if(m_plotPhi) {
    names.push_back(m_name + " phi");
  }
  if(m_plotDiffusionCoefficient && m_isDiffusive){
    names.push_back(m_name + " diffusion_coefficient");
  }
  if(m_plotVelocity && m_isMobile){
    names.push_back("x-Velocity " + m_name);
    names.push_back("y-Velocity " + m_name);
    if(SpaceDim == 3){
      names.push_back("z-Velocity " + m_name);
    }
  }
  if(m_plot_particles)         names.push_back(m_name + " particles");
  if(m_plot_eb_particles)      names.push_back(m_name + " eb_particles");
  if(m_plot_domain_particles)  names.push_back(m_name + " domain_particles");
  if(m_plot_source_particles)  names.push_back(m_name + " source_particles");
  if(m_plot_covered_particles) names.push_back(m_name + " covered_particles");
  if(m_plot_energy_density)    names.push_back(m_name + " energy * phi");
  if(m_plot_average_energy)    names.push_back(m_name + " average_energy");

  return names;
}

int ItoSolver::getNumberOfPlotVariables() const {
  CH_TIME("ItoSolver::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }

  int num_plotvars = 0;
  
  if(m_plotPhi)                num_plotvars += 1;
  if(m_plotDiffusionCoefficient && m_isDiffusive) num_plotvars += 1;
  if(m_plotVelocity && m_isMobile)    num_plotvars += SpaceDim;
  if(m_plot_particles)          num_plotvars = num_plotvars + 1;
  if(m_plot_eb_particles)       num_plotvars = num_plotvars + 1;
  if(m_plot_domain_particles)   num_plotvars = num_plotvars + 1;
  if(m_plot_source_particles)   num_plotvars = num_plotvars + 1;
  if(m_plot_covered_particles)  num_plotvars = num_plotvars + 1;
  if(m_plot_energy_density)     num_plotvars += 1;
  if(m_plot_average_energy)     num_plotvars += 1;

  return num_plotvars;
}

int ItoSolver::getPVRBuffer() const {
  CH_TIME("ItoSolver::getPVRBuffer");
  if(m_verbosity > 5){
    pout() << m_name + "::getPVRBuffer" << endl;
  }

  return m_pvr_buffer;
}

int ItoSolver::getHaloBuffer() const {
  CH_TIME("ItoSolver::getHaloBuffer");
  if(m_verbosity > 5){
    pout() << m_name + "::getHaloBuffer" << endl;
  }

  return m_halo_buffer;
}

void ItoSolver::setPVRBuffer(const int a_buffer) {
  CH_TIME("ItoSolver::setPVRBuffer");
  if(m_verbosity > 5){
    pout() << m_name + "::setPVRBuffer" << endl;
  }

  m_pvr_buffer = a_buffer;
}

void ItoSolver::setHalobuffer(const int a_buffer)  {
  CH_TIME("ItoSolver::setHalobuffer");
  if(m_verbosity > 5){
    pout() << m_name + "::setHalobuffer" << endl;
  }

  m_halo_buffer = a_buffer;
}


size_t ItoSolver::getNumParticles(const WhichContainer a_container, const bool a_local) const{
  CH_TIME("ItoSolver::getNumParticles(string, bool)");
  if(m_verbosity > 5){
    pout() << m_name + "::getNumParticles(string, bool)" << endl;
  }

  size_t N;

  const ParticleContainer<ItoParticle>& particles = m_ParticleContainers.at(a_container);

  if(a_local){
    N = particles.getNumberOfValidParticesLocal();
  }
  else{
    N = particles.getNumberOfValidParticesGlobal();
  }

  return N;
}

void ItoSolver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry> a_computationalGeometry){
  CH_TIME("ItoSolver::setComputationalGeometry");
  if(m_verbosity > 5){
    pout() << m_name + "::setComputationalGeometry" << endl;
  }
  
  m_computationalGeometry = a_computationalGeometry;
}

void ItoSolver::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  CH_TIME("ItoSolver::setAmr");
  if(m_verbosity > 5){
    pout() << m_name + "::setAmr" << endl;
  }

  m_amr = a_amr;
}

void ItoSolver::registerOperators(){
  CH_TIME("ItoSolver::registerOperators");
  if(m_verbosity > 5){
    pout() << m_name + "::registerOperators" << endl;
  }

  if(m_amr.isNull()){
    MayDay::Abort("CdrSolver::registerOperators - need to set AmrMesh!");
  }
  else{
    m_amr->registerOperator(s_eb_coar_ave,     m_realm, m_phase);
    m_amr->registerOperator(s_eb_fill_patch,   m_realm, m_phase);
    m_amr->registerOperator(s_eb_mg_interp,    m_realm, m_phase);
    m_amr->registerOperator(s_eb_copier,       m_realm, m_phase);
    m_amr->registerOperator(s_eb_ghostcloud,   m_realm, m_phase);
    m_amr->registerOperator(s_noncons_div,  m_realm, m_phase);
    
    if(m_redistribute)    m_amr->registerOperator(s_eb_redist,  m_realm, m_phase);
    if(m_halo_buffer > 0) m_amr->registerMask(s_particle_halo, m_halo_buffer, m_realm);
  }
}

void ItoSolver::setPhase(phase::which_phase a_phase){
  CH_TIME("ItoSolver::setPhase");
  if(m_verbosity > 5){
    pout() << m_name + "::setPhase" << endl;
  }

  m_phase = a_phase;
}

void ItoSolver::setVerbosity(const int a_verbosity){
  CH_TIME("ItoSolver::setVerbosity");
  m_verbosity = a_verbosity;
  if(m_verbosity > 5){
    pout() << m_name + "::setVerbosity" << endl;
  }
}

void ItoSolver::setTime(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("ItoSolver::setTime");
  if(m_verbosity > 5){
    pout() << m_name + "::setTime" << endl;
  }

  m_timeStep = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void ItoSolver::initialData(){
  CH_TIME("ItoSolver::initialData");
  if(m_verbosity > 5){
    pout() << m_name + "::initialData" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_ParticleContainers.at(WhichContainer::bulk);

  particles.clearParticles();

  // Add particles, remove the ones that are inside the EB, and then depsit
  particles.addParticles(m_species->getInitialParticles());
  this->removeCoveredParticles(particles, EbRepresentation::ImplicitFunction, 0.0);
  this->depositParticles(m_phi, particles, m_deposition);
}

void ItoSolver::computeLoads(Vector<long int>& a_loads, const DisjointBoxLayout& a_dbl, const int a_level){
  CH_TIME("ItoSolver::computeLoads");
  if(m_verbosity > 5){
    pout() << m_name + "::computeLoads" << endl;
  }

  const ParticleContainer<ItoParticle>& particles = m_ParticleContainers.at(WhichContainer::bulk);

  a_loads.resize(a_dbl.size(),0);
  
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit){
    const int numPart = particles[a_level][dit()].numItems();
    a_loads[dit().intCode()] = numPart;
  }

  // Gather loads globally
#ifdef CH_MPI
  int count = a_loads.size();
  Vector<long int> tmp(count);
  MPI_Allreduce(&(a_loads[0]),&(tmp[0]), count, MPI_LONG, MPI_SUM, Chombo_MPI::comm);
  a_loads = tmp;
#endif
}

void ItoSolver::removeCoveredParticles(const EbRepresentation a_representation, const Real a_tol){
  CH_TIME("ItoSolver::removeCoveredParticles(EbRepresentation, tolerance)");
  if(m_verbosity > 5){
    pout() << m_name + "::removeCoveredParticles(EbRepresentation, tolerance)" << endl;
  }

  this->removeCoveredParticles(WhichContainer::bulk, a_representation, a_tol);
}


void ItoSolver::removeCoveredParticles(const WhichContainer a_container, const EbRepresentation a_representation, const Real a_tol){
  CH_TIME("ItoSolver::removeCoveredParticles(container, EbRepresentation, tolerance)");
  if(m_verbosity > 5){
    pout() << m_name + "::removeCoveredParticles(container, EbRepresentation, tolerance)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);
  
  this->removeCoveredParticles(particles, a_representation, a_tol);
}

void ItoSolver::removeCoveredParticles(ParticleContainer<ItoParticle>& a_particles, const EbRepresentation a_representation, const Real a_tol){
  CH_TIME("ItoSolver::removeCoveredParticles(particles, EbRepresentation)");
  if(m_verbosity > 5){
    pout() << m_name + "::removeCoveredParticles(particles, EbRepresentation)" << endl;
  }

  switch(a_representation){
  case EbRepresentation::ImplicitFunction:
    this->removeCoveredParticles_if(a_particles, a_tol);
    break;
  case EbRepresentation::Discrete:
    this->removeCoveredParticles_discrete(a_particles);
    break;
  case EbRepresentation::Voxel:
    this->removeCoveredParticles_voxels(a_particles);
    break;
  default:
    MayDay::Abort("ItoSolver::removeCoveredParticles - unsupported EB representation requested");
  }
}

void ItoSolver::removeCoveredParticles_if(ParticleContainer<ItoParticle>& a_particles, const Real a_tol){
  CH_TIME("ItoSolver::removeCoveredParticles_if(particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::removeCoveredParticles_if(particles)" << endl;
  }

  const RefCountedPtr<BaseIF>& func = (m_phase == phase::gas) ? m_computationalGeometry->getGasImplicitFunction() : m_computationalGeometry->getSolidImplicitFunction();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    const Real tol               = a_tol*dx;

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      
      List<ItoParticle>& particles = a_particles[lvl][dit()].listItems();

      // Check if particles are outside the implicit function. 
      for (ListIterator<ItoParticle> lit(particles); lit.ok(); ++lit){
	ItoParticle& p = lit();

	const Real f = func->value(p.position());
	if(f > tol) particles.remove(lit);
      }
    }
  }
}


void ItoSolver::transferCoveredParticles_if(ParticleContainer<ItoParticle>& a_src, ParticleContainer<ItoParticle>& a_dst, const Real a_tol){
  CH_TIME("ItoSolver::transferCoveredParticles_if(container, container, tolerance)");
  if(m_verbosity > 5){
    pout() << m_name + "::transferCoveredParticles_if(container, container, tolerance)" << endl;
  }

  const RefCountedPtr<BaseIF>& func = (m_phase == phase::gas) ? m_computationalGeometry->getGasImplicitFunction() : m_computationalGeometry->getSolidImplicitFunction();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    const Real tol               = a_tol*dx;

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      
      List<ItoParticle>& src = a_src[lvl][dit()].listItems();
      List<ItoParticle>& dst = a_dst[lvl][dit()].listItems();

      // Check if particles are outside the implicit function. 
      for (ListIterator<ItoParticle> lit(src); lit.ok(); ++lit){
	ItoParticle& p = lit();

	const Real f = func->value(p.position());
	if(f > tol) {
	  dst.transfer(lit);
	}
      }
    }
  }
}

void ItoSolver::removeCoveredParticles_discrete(ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::removeCoveredParticles_discrete(particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::removeCoveredParticles_discrete(particles)" << endl;
  }

  const RealVect prob_lo = m_amr->getProbLo();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->getDx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      
      List<ItoParticle>& particles = a_particles[lvl][dit()].listItems();

      if(ebisbox.isAllCovered()){ // Box is all covered, remove everything. 
	particles.clear();
      }
      else if(ebisbox.isAllRegular()){ // Do nothing
      }
      else{
	for (ListIterator<ItoParticle> lit(particles); lit.ok(); ++lit){
	  ItoParticle& p = lit();

	  
	  const RealVect rv  = (p.position() - prob_lo)/dx;
	  const IntVect iv   = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));


	  if(ebisbox.isCovered(iv)){
	    particles.remove(lit);
	  }
	  else if(ebisbox.isIrregular(iv)){
	    const VolIndex vof = VolIndex(iv, 0);
	    const RealVect xc  = ebisbox.bndryCentroid(vof);
	    const RealVect nc  = ebisbox.normal(vof);

	    const Real proj    = PolyGeom::dot(p.position() - xc, nc);

	    if(proj < 0.0){
	      particles.remove(lit);
	    }
	  }
	}
      }
    }
  }
}

void ItoSolver::removeCoveredParticles_voxels(ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::removeCoveredParticles_voxels(particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::removeCoveredParticles_voxels(particles)" << endl;
  }

  const RealVect prob_lo = m_amr->getProbLo();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->getDx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      
      List<ItoParticle>& particles = a_particles[lvl][dit()].listItems();

      if(ebisbox.isAllCovered()){ // Box is all covered, remove everything. 
	particles.clear();
      }
      else if(ebisbox.isAllRegular()){ // Do nothing
      }
      else{
	for (ListIterator<ItoParticle> lit(particles); lit.ok(); ++lit){
	  ItoParticle& p = lit();
	  
	  const RealVect rv  = (p.position() - prob_lo)/dx;
	  const IntVect iv   = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

	  if(ebisbox.isCovered(iv)){
	    particles.remove(lit);
	  }
	}
      }
    }
  }
}

void ItoSolver::transferCoveredParticles(const EbRepresentation a_representation, const Real a_tol){
  CH_TIME("ItoSolver::transferCoveredParticles(EbRepresentation, tol)");
  if(m_verbosity > 5){
    pout() << m_name + "::transferCoveredParticles(EbRepresentation, tol)" << endl;
  }

  this->transferCoveredParticles(WhichContainer::bulk, WhichContainer::covered, a_representation, a_tol);
}

void ItoSolver::transferCoveredParticles(const WhichContainer a_containerFrom, const WhichContainer a_containerTo, const EbRepresentation a_representation, const Real a_tol){
  CH_TIME("ItoSolver::transferCoveredParticles(EbRepresentation, string, string, tol)");
  if(m_verbosity > 5){
    pout() << m_name + "::transferCoveredParticles(EbRepresentation, string, string, tol)" << endl;
  }

  ParticleContainer<ItoParticle>& src = this->getParticles(a_containerFrom);
  ParticleContainer<ItoParticle>& dst = this->getParticles(a_containerTo);

  this->transferCoveredParticles(src, dst, a_representation, a_tol);
}

void ItoSolver::transferCoveredParticles(ParticleContainer<ItoParticle>& a_containerFrom,
					 ParticleContainer<ItoParticle>& a_containerTo,
					 const EbRepresentation           a_representation,
					 const Real                        a_tol){
  CH_TIME("ItoSolver::transferCoveredParticles(EbRepresentation, container, container, tol)");
  if(m_verbosity > 5){
    pout() << m_name + "::transferCoveredParticles(EbRepresentation, container, container, tol)" << endl;
  }

  switch(a_representation){
  case EbRepresentation::ImplicitFunction:
    this->transferCoveredParticles_if(a_containerFrom, a_containerTo, a_tol);
    break;
  default:
    MayDay::Abort("ItoSolver::intersectParticles - unsupported EB representation requested");
  }
}

void ItoSolver::intersectParticles(const EbRepresentation a_representation, const bool a_delete){
  CH_TIME("ItoSolver::intersectParticles(EbRepresentation, bool)");
  if(m_verbosity > 5){
    pout() << m_name + "::intersectParticles(EbRepresentation, bool)" << endl;
  }

  this->intersectParticles(WhichContainer::bulk, WhichContainer::eb, WhichContainer::domain, a_representation, a_delete);
}

void ItoSolver::intersectParticles(const WhichContainer       a_particles,
				   const WhichContainer       a_eb_particles,
				   const WhichContainer       a_dom_particles,
				   const EbRepresentation a_representation,				     
				   const bool              a_delete){
  CH_TIME("ItoSolver::intersectParticles(string, string, string, bool, EbRepresentation)");
  if(m_verbosity > 5){
    pout() << m_name + "::intersectParticles(string, string, string, bool, EbRepresentation)" << endl;
  }

  ParticleContainer<ItoParticle>& particles     = this->getParticles(a_particles);
  ParticleContainer<ItoParticle>& eb_particles  = this->getParticles(a_eb_particles);
  ParticleContainer<ItoParticle>& dom_particles = this->getParticles(a_dom_particles);

  this->intersectParticles(particles, eb_particles, dom_particles, a_representation, a_delete);
}


void ItoSolver::intersectParticles(ParticleContainer<ItoParticle>& a_particles,
				   ParticleContainer<ItoParticle>& a_eb_particles,
				   ParticleContainer<ItoParticle>& a_dom_particles,
				   const EbRepresentation           a_representation,
				   const bool                        a_delete){
  CH_TIME("ItoSolver::intersectParticles(container, container, container, EbRepresentation, bool)");
  if(m_verbosity > 5){
    pout() << m_name + "::intersectParticles(container, container, container, EbRepresentation, bool)" << endl;
  }

  switch(a_representation){
  case EbRepresentation::ImplicitFunction:
    this->intersectParticlesIF(a_particles, a_eb_particles, a_dom_particles, a_delete);
    break;
  default:
    MayDay::Abort("ItoSolver::intersectParticles - unsupported EB representation requested");
  }
}

void ItoSolver::intersectParticlesIF(ParticleContainer<ItoParticle>& a_particles,
				     ParticleContainer<ItoParticle>& a_eb_particles,
				     ParticleContainer<ItoParticle>& a_domain_particles,
				     const bool                        a_delete){
  CH_TIME("ItoSolver::intersectParticlesIF(container, container, container, bool)");
  if(m_verbosity > 5){
    pout() << m_name + "::intersectParticlesIF(container, container, container, bool)" << endl;
  }

  const RealVect prob_lo = m_amr->getProbLo();
  const RealVect prob_hi = m_amr->getProbHi();
  const Real     SAFETY  = 1.E-6;

  // This is the implicit function used for intersection tests
  const RefCountedPtr<BaseIF>& impfunc = m_amr->getBaseImplicitFunction(m_phase);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    for (DataIterator dit = m_amr->getGrids(m_realm)[lvl]; dit.ok(); ++dit){

      const Real     dx      = m_amr->getDx()[lvl];
      const EBISBox& ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[lvl][dit()];

      List<ItoParticle>& particles    = a_particles[lvl][dit()].listItems();
      List<ItoParticle>& ebParticles  = a_eb_particles[lvl][dit()].listItems();
      List<ItoParticle>& domParticles = a_domain_particles[lvl][dit()].listItems();

      ebParticles.clear();
      domParticles.clear();

      for (ListIterator<ItoParticle> lit(particles); lit.ok(); ++lit){
	ItoParticle& p = lit();

	const RealVect newPos  = p.position();
	const RealVect oldPos  = p.oldPosition();
	const RealVect path    = newPos - oldPos;
	const Real     pathLen = path.vectorLength();

	// Check if we should check of different types of boundary intersections. These are checp initial tests that allow
	// us to skip intersection tests for some Photons.
	bool checkEB  = false;
	bool checkDom = false;

	if(!impfunc.isNull()){
	  checkEB = true;
	}
	for (int dir = 0; dir < SpaceDim; dir++){
	  if(newPos[dir] < prob_lo[dir] || newPos[dir] > prob_hi[dir]){ // Checks if we crossed a domain boundary. 
	    checkDom = true; 
	  }
	}

	// Must do intersection test on at least one of these. 
	if(checkEB || checkDom){ 
	  Real dom_s = 1.E99;
	  Real eb_s  = 1.E99;

	  bool contact_domain = false;
	  bool contact_eb     = false;
	      
	  if(checkDom) contact_domain = ParticleOps::domainIntersection(oldPos, newPos, path, prob_lo, prob_hi, dom_s);
#if 0
	  if(checkEB)  contact_eb     = ParticleOps::ebIntersectionBisect(impfunc, oldPos, newPos, pathLen, dx, eb_s);
#else
	  if(checkEB)  contact_eb     = ParticleOps::ebIntersectionRaycast(impfunc, oldPos, newPos, 1.E-10*dx, eb_s);
#endif
	  
	  if(contact_eb || contact_domain){ // Particle trajectory crossed something. 
	    if(eb_s < dom_s){ // It was the EB first. 
	      p.position() = oldPos + eb_s*path;
	      if(a_delete){
		ebParticles.transfer(lit);
	      }
	      else{
		ebParticles.add(lit());
	      }
	    }
	    else{ // It was the domain boundary. 
	      p.position() = oldPos + Max(0.0,dom_s-SAFETY)*path;
	      if(a_delete){
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

void ItoSolver::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("ItoSolver::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }

  // Reallocate mesh data
  m_amr->reallocate(m_phi,         m_phase, a_lmin);
  m_amr->reallocate(m_scratch,       m_phase, a_lmin);
  m_amr->reallocate(m_mobility_func, m_phase, a_lmin);
  m_amr->reallocate(m_depositionNC,  m_phase, a_lmin);
  m_amr->reallocate(m_massDiff,      m_phase, a_lmin);
  
  // Only allocate memory if we actually have a mobile solver
  if(m_isMobile){
    m_amr->reallocate(m_velo_func, m_phase, a_lmin);
  }
  else{ 
    m_amr->allocatePointer(m_velo_func);
  }

  // Only allocate memory if we actually a diffusion solver
  if(m_isDiffusive){
    m_amr->reallocate(m_faceCenteredDiffusionCoefficient_cell, m_phase, a_lmin);
  }
  else{
    m_amr->allocatePointer(m_faceCenteredDiffusionCoefficient_cell);
  }

  // Particle data regrids
  const Vector<DisjointBoxLayout>& grids = m_amr->getGrids(m_realm);
  const Vector<ProblemDomain>& domains   = m_amr->getDomains();
  const Vector<Real>& dx                 = m_amr->getDx();
  const Vector<int>& ref_rat             = m_amr->getRefinementRatios();

  for (auto& container : m_ParticleContainers){
    ParticleContainer<ItoParticle>& particles = container.second;
    particles.regrid(grids, domains, dx, ref_rat, a_lmin, a_newFinestLevel);
  }
}

void ItoSolver::setSpecies(RefCountedPtr<ItoSpecies> a_species){
  CH_TIME("ItoSolver::setSpecies");
  if(m_verbosity > 5){
    pout() << m_name + "::setSpecies" << endl;
  }
  
  m_species   = a_species;
  m_name      = a_species->getName();
  m_isDiffusive = m_species->isDiffusive();
  m_isMobile    = m_species->isMobile();
}

void ItoSolver::allocateInternals(){
  CH_TIME("ItoSolver::allocateInternals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocateInternals" << endl;
  }
  
  const int ncomp = 1;

  m_amr->allocate(m_phi,         m_realm, m_phase, ncomp);
  m_amr->allocate(m_scratch,       m_realm, m_phase, ncomp);
  m_amr->allocate(m_mobility_func, m_realm, m_phase, ncomp);
  m_amr->allocate(m_depositionNC,  m_realm, m_phase, ncomp);
  m_amr->allocate(m_massDiff,      m_realm, m_phase, ncomp);

  // Only allocate memory for velocity if we actually have a mobile solver
  if(m_isMobile){
    m_amr->allocate(m_velo_func, m_realm, m_phase, SpaceDim);
  }
  else{ 
    m_amr->allocatePointer(m_velo_func);
  }

  // Only allocate memory if we actually have a diffusion solver
  if(m_isDiffusive){
    m_amr->allocate(m_faceCenteredDiffusionCoefficient_cell, m_realm, m_phase, 1);
  }
  else{
    m_amr->allocatePointer(m_faceCenteredDiffusionCoefficient_cell);
  }

  m_ParticleContainers.emplace(WhichContainer::bulk,    ParticleContainer<ItoParticle>());
  m_ParticleContainers.emplace(WhichContainer::eb,      ParticleContainer<ItoParticle>());
  m_ParticleContainers.emplace(WhichContainer::domain,  ParticleContainer<ItoParticle>());
  m_ParticleContainers.emplace(WhichContainer::source,  ParticleContainer<ItoParticle>());
  m_ParticleContainers.emplace(WhichContainer::covered, ParticleContainer<ItoParticle>());
  m_ParticleContainers.emplace(WhichContainer::scratch, ParticleContainer<ItoParticle>());

  for (auto& container : m_ParticleContainers){
    m_amr->allocate(container.second, m_pvr_buffer, m_realm);
  }
}

#ifdef CH_USE_HDF5  
void ItoSolver::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("ItoSolver::writeCheckpointLevel");
  if(m_verbosity > 5){
    pout() << m_name + "::writeCheckpointLevel" << endl;
  }

  // Write state. 
  write(a_handle, *m_phi[a_level], m_name);

  // Write particles.
  if(m_checkpointing == WhichCheckpoint::particles){
    this->writeCheckPointLevelParticles(a_handle, a_level);
  }
  else{ // In this case we need to write the number of physical particles in a grid cell. 
    this->writeCheckPointLevelFluid(a_handle, a_level);
  }
}
#endif

#ifdef CH_USE_HDF5  
void ItoSolver::writeCheckPointLevelParticles(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("ItoSolver::writeCheckPointLevelParticles");
  if(m_verbosity > 5){
    pout() << m_name + "::writeCheckPointLevelParticles" << endl;
  }

  const int halo        = 0;
  const std::string str = m_name + "_particles";

  ParticleContainer<SimpleItoParticle> RealmParticles;
  m_amr->allocate(RealmParticles,  m_pvr_buffer, m_realm);

  const ParticleContainer<ItoParticle>& myParticles = this->getParticles(WhichContainer::bulk);

  // Make ItoParticle into SimpleItoParticle. This saves a shitload of disk space. 
  for (DataIterator dit(m_amr->getGrids(m_realm)[a_level]); dit.ok(); ++dit){
    List<SimpleItoParticle>& other_particles = (RealmParticles[a_level])[dit()].listItems();

    other_particles.clear();
      
    for (ListIterator<ItoParticle> lit(myParticles[a_level][dit()].listItems()); lit.ok(); ++lit){
      other_particles.append(SimpleItoParticle(lit().mass(), lit().position(), lit().energy()));
    }
  }

  writeParticlesToHDF(a_handle, RealmParticles[a_level], str);
}
#endif

#ifdef CH_USE_HDF5  
void ItoSolver::writeCheckPointLevelFluid(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("ItoSolver::writeCheckPointLevelFluid");
  if(m_verbosity > 5){
    pout() << m_name + "::writeCheckPointLevelFluid" << endl;
  }

  const std::string str = m_name + "_particles";

  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[a_level];
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
    
  const int comp     = 0;
  const int ncomp    = 1;
  const RealVect dx  = m_amr->getDx()[a_level]*RealVect::Unit;
  const RealVect plo = m_amr->getProbLo();
  const IntVect ghos = IntVect::Zero;

  // Relevant container
  const ParticleContainer<ItoParticle>& particles = m_ParticleContainers.at(WhichContainer::bulk);

  // Make something that can hold the particle numbers (stored as a Real)
  EBCellFactory fact(ebisl);
  LevelData<EBCellFAB> particleNumbers(dbl, ncomp, ghos, fact);
  EBLevelDataOps::setVal(particleNumbers, 0.0);

  // Now go through the grid and add the number of particles in each cell
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    BaseFab<Real>& pNum = particleNumbers[dit()].getSingleValuedFAB(); // No multivalued cells please. 

    // Get cell particles
    BinFab<ItoParticle> pCel(dbl.get(dit()), dx, plo);
    pCel.addItems(particles[a_level][dit()].listItems());
				 
    for (BoxIterator bit(dbl.get(dit())); bit.ok(); ++bit){
      pNum(bit(), comp) = 0.0;
      for (ListIterator<ItoParticle> lit(pCel(bit(), comp)); lit.ok(); ++lit){
	pNum(bit()) += lit().mass();
      }
    }
  }

  write(a_handle, particleNumbers, str);
}
#endif

#ifdef CH_USE_HDF5  
void ItoSolver::readCheckpointLevel(HDF5Handle& a_handle, const int a_level){
  CH_TIME("ItoSolver::readCheckpointLevel");
  if(m_verbosity > 5){
    pout() << m_name + "::readCheckpointLevel" << endl;
  }

  // Read state vector
  read<EBCellFAB>(a_handle, *m_phi[a_level], m_name, m_amr->getGrids(m_realm)[a_level], Interval(0,0), false);

  // Read particles.
  const std::string str = m_name + "_particles";
  if(m_checkpointing == WhichCheckpoint::particles){

    // Get SimpleItoParticles from data file
    Vector<RefCountedPtr<ParticleData<SimpleItoParticle> > > simpleParticles;
    m_amr->allocate(simpleParticles, m_realm);
    readParticlesFromHDF(a_handle, *simpleParticles[a_level], str);

    ParticleContainer<ItoParticle>& particles = m_ParticleContainers.at(WhichContainer::bulk);

    // Make SimpleItoParticles into ito_ particles
    for (DataIterator dit(m_amr->getGrids(m_realm)[a_level]); dit.ok(); ++dit){
      List<ItoParticle>& particlesDit = particles[a_level][dit()].listItems();
      
      for (ListIterator<SimpleItoParticle> lit((*simpleParticles[a_level])[dit()].listItems()); lit.ok(); ++lit){
	particlesDit.append(ItoParticle(lit().mass(), lit().position(), RealVect::Zero, 0.0, 0.0, lit().energy()));
      }
    }
  }
  else{

    // Read particles per cell
    EBAMRCellData ppc;
    m_amr->allocate(ppc, m_realm, m_phase, 1, 0);
    read<EBCellFAB>(a_handle, *ppc[a_level], str, m_amr->getGrids(m_realm)[a_level], Interval(0,0), false);
    
    // Restart particles
    this->restartParticles(*ppc[a_level], a_level);
  }
    
}
#endif

void ItoSolver::restartParticles(LevelData<EBCellFAB>& a_num_particles, const int a_level){
  CH_TIME("ItoSolver::restartParticles");
  if(m_verbosity > 5){
    pout() << m_name + "::restartParticles" << endl;
  }

  const Real dx          = m_amr->getDx()[a_level];
  const RealVect prob_lo = m_amr->getProbLo();

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];

  ParticleContainer<ItoParticle>& particles = m_ParticleContainers.at(WhichContainer::bulk);

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box           = dbl.get(dit());
    const EBISBox& ebisbox   = m_amr->getEBISLayout(m_realm, m_phase)[a_level][dit()];
    const BaseFab<Real>& ppc = a_num_particles[dit()].getSingleValuedFAB();

    // Clear just to be safe
    List<ItoParticle>& myParticles = particles[a_level][dit()].listItems();
    myParticles.clear();

    // Do regular cells
    for (BoxIterator bit(box); bit.ok(); ++bit){
      if(ebisbox.isRegular(bit())){

	// Compute weights and remainder
	const unsigned long long numParticles = (unsigned long long) llround(ppc(bit()));
	unsigned long long w, N, r;
	DataOps::computeParticleWeights(w, N, r, numParticles, m_ppc_restart);

	const RealVect minLo = -0.5*RealVect::Unit;
	const RealVect minHi =  0.5*RealVect::Unit;
	const RealVect norma = RealVect::Zero;
	const RealVect centr = RealVect::Zero;
	const RealVect pos   = prob_lo + (RealVect(bit()) + 0.5*RealVect::Unit)*dx;
	const Real kappa     = 1.0;


	// Now add N partices. If r > 0 we add another one with weight w + r
	for (unsigned long long i = 0; i < N; i++){
	  const RealVect particlePosition = this->randomPosition(pos, minLo, minHi, centr, norma, dx, kappa);
	  const Real particleWeight       = (Real) w;
	  
	  myParticles.add(ItoParticle(particleWeight, particlePosition));
	}

	if(r > 0){
	  const RealVect particlePosition = this->randomPosition(pos, minLo, minHi, centr, norma, dx, kappa);
	  const Real particleWeight       = (Real) (w+r);

	  myParticles.add(ItoParticle(particleWeight, particlePosition));
	}
      }
    }

    // Do irregular cells
    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_level])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const IntVect iv    = vof.gridIndex();
      const RealVect cent = ebisbox.bndryCentroid(vof);
      const RealVect norm = ebisbox.normal(vof);
      const RealVect pos  = prob_lo + dx*(RealVect(iv) + 0.5*RealVect::Unit);
      const Real kappa    = ebisbox.volFrac(vof);


      // Compute a small box that encloses the cut-cell volume
      RealVect minLo = -0.5*RealVect::Unit;
      RealVect minHi =  0.5*RealVect::Unit;
      if(kappa < 1.0){
	DataOps::computeMinValidBox(minLo, minHi, norm, cent);
      }

      // Compute weights and remainder
      const unsigned long long numParticles = (unsigned long long) llround(ppc(iv));
      unsigned long long w, N, r;
      DataOps::computeParticleWeights(w, N, r, numParticles, m_ppc_restart);

      // Now add N partices. If r > 0 we add another one with weight w + r
      for (unsigned long long i = 0; i < N; i++){
	const RealVect particlePosition = this->randomPosition(pos, minLo, minHi, cent, norm, dx, kappa);
	const Real particleWeight       = (Real) w;
	  
	myParticles.add(ItoParticle(particleWeight, particlePosition));
      }

      if(r > 0){
	const RealVect particlePosition = this->randomPosition(pos, minLo, minHi, cent, norm, dx, kappa);
	const Real particleWeight       = (Real) (w+r);

	myParticles.add(ItoParticle(particleWeight, particlePosition));
      }
    }
  }
}

void ItoSolver::writePlotData(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("ItoSolver::writePlotData");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotData" << endl;
  }

  // Write phi
  if(m_plotPhi){
    const Interval src(0, 0);
    const Interval dst(a_comp, a_comp);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      if(m_realm == a_output.getRealm()){
	m_phi[lvl]->localCopyTo(src, *a_output[lvl], dst);
      }
      else{
	m_phi[lvl]->copyTo(src, *a_output[lvl], dst);
      }
    }
    //DataOps::setCoveredValue(a_output, a_comp, 0.0);
    a_comp++;
  }

  // Plot diffusion coefficient
  if(m_plotDiffusionCoefficient && m_isDiffusive){
    const Interval src(0, 0);
    const Interval dst(a_comp, a_comp);
    
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      if(m_realm == a_output.getRealm()){
	m_faceCenteredDiffusionCoefficient_cell[lvl]->localCopyTo(src, *a_output[lvl], dst);
      }
      else{
	m_faceCenteredDiffusionCoefficient_cell[lvl]->copyTo(src, *a_output[lvl], dst);
      }
    }
    DataOps::setCoveredValue(a_output, a_comp, 0.0);
    a_comp++;
  }

  // Write velocities
  if(m_plotVelocity && m_isMobile){
    const int ncomp = SpaceDim;
    const Interval src(0, ncomp-1);
    const Interval dst(a_comp, a_comp + ncomp-1);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      if(m_realm == a_output.getRealm()){
	m_velo_func[lvl]->localCopyTo(src, *a_output[lvl], dst);
      }
      else{
	m_velo_func[lvl]->copyTo(src, *a_output[lvl], dst);
      }
    }

    for (int c = 0; c < SpaceDim; c++){
      DataOps::setCoveredValue(a_output, a_comp + c, 0.0);
    }

    a_comp += ncomp;
  }

  if(m_plot_particles){
    this->depositParticles(m_scratch, m_ParticleContainers.at(WhichContainer::bulk), m_plot_deposition);
    this->writeData(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_eb_particles){
    this->depositParticles(m_scratch, m_ParticleContainers.at(WhichContainer::eb), m_plot_deposition);
    this->writeData(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_domain_particles){
    this->depositParticles(m_scratch, m_ParticleContainers.at(WhichContainer::domain), m_plot_deposition);
    this->writeData(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_source_particles){
    this->depositParticles(m_scratch, m_ParticleContainers.at(WhichContainer::source), m_plot_deposition);
    this->writeData(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_energy_density){
    this->depositEnergyDensity(m_scratch, m_ParticleContainers.at(WhichContainer::bulk), m_plot_deposition);
    this->writeData(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_average_energy){
    this->sortParticlesByCell(WhichContainer::bulk);
    this->computeAverageEnergy(m_scratch, m_ParticleContainers.at(WhichContainer::bulk));
    this->writeData(a_output, a_comp, m_scratch,  false);
    this->sortParticlesByPatch(WhichContainer::bulk);
  }
}

void ItoSolver::writeData(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp){
  CH_TIME("ItoSolver::writeData");
  if(m_verbosity > 5){
    pout() << m_name + "::writeData" << endl;
  }

  const int comp = 0;
  const int ncomp = a_data[0]->nComp();

  const Interval src_interv(0, ncomp-1);
  const Interval dst_interv(a_comp, a_comp + ncomp - 1);

  // Copy data onto scratch
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, ncomp);
  DataOps::copy(scratch, a_data);

  // Interp if we should
  if(a_interp){
    m_amr->interpToCentroids(scratch, m_realm, phase::gas);
  }

  m_amr->averageDown(scratch, m_realm, m_phase);
  m_amr->interpGhost(scratch, m_realm, m_phase);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    if(m_realm == a_output.getRealm()){
      scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else{
      scratch[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }

  DataOps::setCoveredValue(a_output, a_comp, 0.0);

  a_comp += ncomp;
}

void ItoSolver::setMassToConductivity(ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::unsetMassToConductivity");
  if(m_verbosity > 5){
    pout() << m_name + "::unsetMassToConductivity" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ItoParticle>& particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<ItoParticle> lit(particles); lit.ok(); ++lit){
	ItoParticle& p = lit();
	p.tmp()   = p.mass();
	p.mass() *= p.mobility();
      }
    }
  }
}

void ItoSolver::unsetMassToConductivity(ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::unsetMassToConductivity");
  if(m_verbosity > 5){
    pout() << m_name + "::unsetMassToConductivity" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ItoParticle>& particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<ItoParticle> lit(particles); lit.ok(); ++lit){
	ItoParticle& p = lit();
	p.mass() = p.tmp();
      }
    }
  }
}

void ItoSolver::setMassToDiffusivity(ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::unsetMassToDiffusivity");
  if(m_verbosity > 5){
    pout() << m_name + "::unsetMassToDiffusivity" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ItoParticle>& particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<ItoParticle> lit(particles); lit.ok(); ++lit){
	ItoParticle& p = lit();
	p.tmp()   = p.mass();
	p.mass() *= p.diffusion();
      }
    }
  }
}

void ItoSolver::unsetMassToDiffusivity(ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::unsetMassToDiffusivity");
  if(m_verbosity > 5){
    pout() << m_name + "::unsetMassToDiffusivity" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ItoParticle>& particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<ItoParticle> lit(particles); lit.ok(); ++lit){
	ItoParticle& p = lit();
	p.mass() = p.tmp();
      }
    }
  }
}

void ItoSolver::setMassToEnergy(ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::setMassToEnergy");
  if(m_verbosity > 5){
    pout() << m_name + "::setMassToEnergy" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ItoParticle>& particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<ItoParticle> lit(particles); lit.ok(); ++lit){
	ItoParticle& p = lit();
	p.tmp()   = p.mass();
	p.mass() *= p.energy();
      }
    }
  }
}

void ItoSolver::unsetMassToEnergy(ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::unsetMassToEnergy");
  if(m_verbosity > 5){
    pout() << m_name + "::unsetMassToEnergy" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ItoParticle>& particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<ItoParticle> lit(particles); lit.ok(); ++lit){
	ItoParticle& p = lit();
	p.mass() = p.tmp();
      }
    }
  }
}

void ItoSolver::depositConductivity(){
  CH_TIME("ItoSolver::depositConductivity()");
  if(m_verbosity > 5){
    pout() << m_name + "::depositConductivity()" << endl;
  }

  this->depositConductivity(m_phi, m_ParticleContainers.at(WhichContainer::bulk));
}

void ItoSolver::depositConductivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::depositConductivity(state, particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::depositConductivity(state, particles)" << endl;
  }

  this->depositConductivity(a_phi, a_particles, m_deposition);
}

void ItoSolver::depositConductivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles, const DepositionType a_deposition){
  CH_TIME("ItoSolver::depositConductivity(state, particles, deposition_type)");
  if(m_verbosity > 5){
    pout() << m_name + "::depositConductivity(state, particles, deposition_type)" << endl;
  }

  this->setMassToConductivity(a_particles);                 // Make mass = mass*mu
  this->depositParticles(a_phi, a_particles, a_deposition); // Deposit mass*mu
  this->unsetMassToConductivity(a_particles);               // Make mass = mass/mu
}

void ItoSolver::depositDiffusivity(){
  CH_TIME("ItoSolver::depositDiffusivity()");
  if(m_verbosity > 5){
    pout() << m_name + "::depositDiffusivity()" << endl;
  }

  this->depositDiffusivity(m_phi, m_ParticleContainers.at(WhichContainer::bulk));
}

void ItoSolver::depositDiffusivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::depositDiffusivity(state, particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::depositDiffusivity(state, particles)" << endl;
  }

  this->depositDiffusivity(a_phi, a_particles, m_deposition);
}

void ItoSolver::depositDiffusivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles, const DepositionType a_deposition){
  CH_TIME("ItoSolver::depositDiffusivity(state, particles, deposition_type)");
  if(m_verbosity > 5){
    pout() << m_name + "::depositDiffusivity(state, particles, deposition_type)" << endl;
  }

  this->setMassToDiffusivity(a_particles);                                  // Make mass = mass*D
  this->depositParticles(a_phi, a_particles, a_deposition); // Deposit mass*D
  this->unsetMassToDiffusivity(a_particles);                                // Make mass = mass/D
}

void ItoSolver::depositEnergyDensity(){
  CH_TIME("ItoSolver::depositEnergyDensity()");
  if(m_verbosity > 5){
    pout() << m_name + "::depositEnergyDensity()" << endl;
  }

  this->depositEnergyDensity(m_phi, m_ParticleContainers.at(WhichContainer::bulk));
}

void ItoSolver::depositEnergyDensity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::depositEnergyDensity(state, particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::depositEnergyDensity(state, particles)" << endl;
  }

  this->depositEnergyDensity(a_phi, a_particles, m_deposition);
}

void ItoSolver::depositEnergyDensity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles, const DepositionType a_deposition){
  CH_TIME("ItoSolver::depositEnergyDensity(state, particles, deposition_type)");
  if(m_verbosity > 5){
    pout() << m_name + "::depositEnergyDensity(state, particles, deposition_type)" << endl;
  }

  this->setMassToEnergy(a_particles);                       // Make mass = mass*E
  this->depositParticles(a_phi, a_particles, a_deposition); // Deposit mass*E
  this->unsetMassToEnergy(a_particles);                     // Make mass = mass/E
}


void ItoSolver::computeAverageMobility(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::computeAverageMobility(state, particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeAverageMobility(state, particles)" << endl;
  }

  constexpr int comp = 0;

  DataOps::setValue(a_phi, 0.0);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const BinFab<ItoParticle>& particles = a_particles.getCellParticles(lvl, dit());

      BaseFab<Real>& state = (*a_phi[lvl])[dit()].getSingleValuedFAB();

      for (BoxIterator bit(dbl[dit()]); bit.ok(); ++bit){
	const IntVect iv = bit();

	const List<ItoParticle>& listParticles = particles(iv, comp);

	if(listParticles.length() > 0){
	  Real M  = 0.0;
	  Real MU = 0.0;
	  for (ListIterator<ItoParticle> lit(listParticles); lit.ok(); ++lit){
	    const ItoParticle& p = lit();
	    const Real m  = p.mass();
	    const Real mu = p.mobility();

	    M  += m;
	    MU += m*mu;
	  }

	  state(iv, comp) = MU/M;
	}
      }
    }
  }
}

void ItoSolver::computeAverageDiffusion(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::computeAverageDiffusion(state, particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeAverageDiffusion(state, particles)" << endl;
  }

  constexpr int comp = 0;

  DataOps::setValue(a_phi, 0.0);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const BinFab<ItoParticle>& particles = a_particles.getCellParticles(lvl, dit());

      BaseFab<Real>& state = (*a_phi[lvl])[dit()].getSingleValuedFAB();

      for (BoxIterator bit(dbl[dit()]); bit.ok(); ++bit){
	const IntVect iv = bit();

	const List<ItoParticle>& listParticles = particles(iv, comp);

	if(listParticles.length() > 0){
	  Real M = 0.0;
	  Real D = 0.0;
	  for (ListIterator<ItoParticle> lit(listParticles); lit.ok(); ++lit){
	    const ItoParticle& p = lit();
	    const Real m = p.mass();
	    const Real d = p.diffusion();

	    M += m;
	    D += m*d;
	  }

	  state(iv, comp) = D/M;
	}
      }
    }
  }
}

void ItoSolver::computeAverageEnergy(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::computeAverageEnergy(state, particles)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeAverageEnergy(state, particles)" << endl;
  }

  constexpr int comp = 0;

  DataOps::setValue(a_phi, 0.0);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const BinFab<ItoParticle>& particles = a_particles.getCellParticles(lvl, dit());

      BaseFab<Real>& state = (*a_phi[lvl])[dit()].getSingleValuedFAB();

      for (BoxIterator bit(dbl[dit()]); bit.ok(); ++bit){
	const IntVect iv = bit();

	const List<ItoParticle>& listParticles = particles(iv, comp);

	if(listParticles.length() > 0){
	  Real M = 0.0;
	  Real E = 0.0;
	  for (ListIterator<ItoParticle> lit(listParticles); lit.ok(); ++lit){
	    const ItoParticle& p = lit();
	    const Real m = p.mass();
	    const Real e = p.energy();

	    M += m;
	    E += m*e;
	  }

	  state(iv, comp) = E/M;
	}
      }
    }
  }
}

void ItoSolver::depositParticles(){
  CH_TIME("ItoSolver::depositParticles");
  if(m_verbosity > 5){
    pout() << m_name + "::depositParticles" << endl;
  }

  this->depositParticles(WhichContainer::bulk);
}

void ItoSolver::depositParticles(const WhichContainer a_container){
  CH_TIME("ItoSolver::depositParticles(container)");
  if(m_verbosity > 5){
    pout() << m_name + "::depositParticles(container)" << endl;
  }

  this->depositParticles(m_phi, m_ParticleContainers.at(a_container), m_deposition);
}

void ItoSolver::depositNonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC){
  CH_TIME("ItoSolver::depositNonConservative");
  if(m_verbosity > 5){
    pout() << m_name + "::depositNonConservative" << endl;
  }

  const std::string cur_Realm = a_depositionNC.getRealm();

  if(m_blendConservation){
    const IrregAmrStencil<NonConservativeDivergenceStencil>& stencils = m_amr->getNonConservativeDivergenceStencils(cur_Realm, m_phase);
    stencils.apply(a_depositionNC, a_depositionKappaC);
  }
  else{
    DataOps::setValue(a_depositionNC, 0.0);
  }
}

void ItoSolver::depositHybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_massDifference, const EBAMRIVData& a_depositionNC){
  CH_TIME("ItoSolver::depositHybrid");
  if(m_verbosity > 5){
    pout() << m_name + "::depositHybrid" << endl;
  }

  const std::string cur_Realm = a_depositionH.getRealm();

  const int comp  = 0;
  const int ncomp = 1;


  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(cur_Realm)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(cur_Realm, m_phase)[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      EBCellFAB& divH               = (*a_depositionH[lvl])[dit()];  // On input, this contains kappa*depositionWeights
      BaseIVFAB<Real>& deltaM       = (*a_massDifference[lvl])[dit()];
      const BaseIVFAB<Real>& divNC  = (*a_depositionNC[lvl])[dit()]; 
      const EBISBox& ebisbox        = ebisl[dit()];

      VoFIterator& vofit = (*m_amr->getVofIterator(cur_Realm, m_phase)[lvl])[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
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


void ItoSolver::incrementRedist(const EBAMRIVData& a_massDifference){
  CH_TIME("ItoSolver::incrementRedist");
  if(m_verbosity > 5){
    pout() << m_name + "::incrementRedist" << endl;
  }

  const std::string cur_Realm = a_massDifference.getRealm();

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(cur_Realm)[lvl];
    
    EBLevelRedist& level_redist = *(m_amr->getLevelRedist(cur_Realm, m_phase)[lvl]);
    level_redist.setToZero();

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      level_redist.increment((*a_massDifference[lvl])[dit()], dit(), interv);
    }
  }
}

void ItoSolver::levelRedist(EBAMRCellData& a_phi){
  CH_TIME("ItoSolver::levelRedist");
  if(m_verbosity > 5){
    pout() << m_name + "::levelRedist" << endl;
  }

  const std::string cur_Realm = a_phi.getRealm();

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    EBLevelRedist& level_redist = *(m_amr->getLevelRedist(cur_Realm, m_phase)[lvl]);
    level_redist.redistribute(*a_phi[lvl], interv);
    level_redist.setToZero();
  }
}

void ItoSolver::coarseFineIncrement(const EBAMRIVData& a_massDifference){
  CH_TIME("ItoSolver::coarseFineIncrement");
  if(m_verbosity > 5){
    pout() << m_name + "::coarseFineIncrement" << endl;
  }

  const std::string cur_Realm = a_massDifference.getRealm();

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(0,0);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(cur_Realm)[lvl];

    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->getFineToCoarRedist(cur_Realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->getCoarToFineRedist(cur_Realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->getCoarToCoarRedist(cur_Realm, m_phase)[lvl];

    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < 0;

    if(has_coar){
      fine2coar_redist->setToZero();

    }
    if(has_fine){
      coar2fine_redist->setToZero();
      coar2coar_redist->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      if(has_coar){
	fine2coar_redist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
      }

      if(has_fine){
	coar2fine_redist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
	coar2coar_redist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
      }
    }
  }
}


void ItoSolver::coarseFineRedistribution(EBAMRCellData& a_phi){
  CH_TIME("ItoSolver::coarseFineRedistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::coarseFineRedistribution" << endl;
  }

  const std::string cur_Realm = a_phi.getRealm();

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);


  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx       = m_amr->getDx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->getCoarToFineRedist(cur_Realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->getCoarToCoarRedist(cur_Realm, m_phase)[lvl];
    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->getFineToCoarRedist(cur_Realm, m_phase)[lvl];
    if(has_coar){
      fine2coar_redist->redistribute(*a_phi[lvl-1], interv);
      fine2coar_redist->setToZero();
    }

    if(has_fine){
      coar2fine_redist->redistribute(*a_phi[lvl+1], interv);
      coar2coar_redist->redistribute(*a_phi[lvl],   interv);

      coar2fine_redist->setToZero();
      coar2coar_redist->setToZero();
    }
  }
}

void ItoSolver::depositWeights(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::depositWeights");
  if(m_verbosity > 5){
    pout() << m_name + "::depositWeights" << endl;
  }

  this->depositParticles(a_phi, m_ParticleContainers.at(WhichContainer::bulk), DepositionType::NGP);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    DataOps::scale(*a_phi[lvl], pow(m_amr->getDx()[lvl], SpaceDim));
  }
}

void ItoSolver::addParticles(ListBox<ItoParticle>& a_part, const int a_lvl, const DataIndex a_dit, const bool a_destructive){
  CH_TIME("ItoSolver::addParticles(lvl, dit)");
  if(m_verbosity > 5){
    pout() << m_name + "::addParticles(lvl, dit)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_ParticleContainers.at(WhichContainer::bulk);

  ListBox<ItoParticle>& my_particles = particles[a_lvl][a_dit];

  if(a_destructive){
    my_particles.addItemsDestructive(a_part.listItems());
  }
  else{
    my_particles.addItems(a_part.listItems());
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

void ItoSolver::preRegrid(const int a_base, const int a_oldFinestLevel){
  CH_TIME("ItoSolver::preRegrid");
  if(m_verbosity > 5){
    pout() << m_name + "::preRegrid" << endl;
  }

  for (auto& container : m_ParticleContainers){
    ParticleContainer<ItoParticle>& particles = container.second;

    particles.preRegrid(a_base);
  }
}

ParticleContainer<ItoParticle>& ItoSolver::getParticles(const WhichContainer a_container){
  return m_ParticleContainers.at(a_container);
}

const ParticleContainer<ItoParticle>& ItoSolver::getParticles(const WhichContainer a_container) const {
  return m_ParticleContainers.at(a_container);
}

EBAMRCellData& ItoSolver::getPhi(){
  CH_TIME("ItoSolver::getPhi");
  if(m_verbosity > 5){
    pout() << m_name + "::getPhi" << endl;
  }

  return m_phi;
}

EBAMRCellData& ItoSolver::getVelocityFunction(){
  CH_TIME("ItoSolver::getVelocityFunction");
  if(m_verbosity > 5){
    pout() << m_name + "::getVelocityFunction" << endl;
  }

  return m_velo_func;
}

EBAMRCellData& ItoSolver::getDiffusionFunction(){
  CH_TIME("ItoSolver::getDiffusionFunction");
  if(m_verbosity > 5){
    pout() << m_name + "::getDiffusionFunction" << endl;
  }

  return m_faceCenteredDiffusionCoefficient_cell;
}

EBAMRCellData& ItoSolver::getScratch(){
  CH_TIME("ItoSolver::getScratch");
  if(m_verbosity > 5){
    pout() << m_name + "::getScratch" << endl;
  }

  return m_scratch;
}

EBAMRCellData& ItoSolver::getMobilityFunction(){
  CH_TIME("ItoSolver::getMobilityFunction");
  if(m_verbosity > 5){
    pout() << m_name + "::getMobilityFunction" << endl;
  }

  return m_mobility_func;
}

void ItoSolver::setDiffusionFunction(const Real a_diffusionCoefficient){
  CH_TIME("ItoSolver::setDiffusionFunction");
  if(m_verbosity > 5){
    pout() << m_name + "::setDiffusionFunction" << endl;
  }

  DataOps::setValue(m_faceCenteredDiffusionCoefficient_cell, a_diffusionCoefficient);
}

void ItoSolver::setVelocityFunction(const RealVect a_vel){
  CH_TIME("ItoSolver::setVelocityFunction");
  if(m_verbosity > 5){
    pout() << m_name + "::setVelocityFunction" << endl;
  }

  for (int comp = 0; comp < SpaceDim; comp++){
    DataOps::setValue(m_velo_func, a_vel[comp], comp);
  }
}

void ItoSolver::set_mobility(const Real a_mobility){
  CH_TIME("ItoSolver::set_mobility");
  if(m_verbosity > 5){
    pout() << m_name + "::set_mobility" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::bulk);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<ItoParticle>& particlesDit = particles[lvl][dit()].listItems();

      for (ListIterator<ItoParticle> lit(particlesDit); lit.ok(); ++lit){
	lit().mobility() = a_mobility;
      }
    }
  }
}

void ItoSolver::interpolateVelocities(){
  CH_TIME("ItoSolver::interpolateVelocities");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolateVelocities" << endl;
  }

  if(m_isMobile){
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	this->interpolateVelocities(lvl, dit());
      }
    }
  }
}

void ItoSolver::interpolateVelocities(const int a_lvl, const DataIndex& a_dit){
  CH_TIME("ItoSolver::interpolateVelocities");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolateVelocities" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_ParticleContainers.at(WhichContainer::bulk);

  if(m_isMobile){
    const EBCellFAB& velo_func = (*m_velo_func[a_lvl])[a_dit];
    const EBISBox& ebisbox     = velo_func.getEBISBox();
    const FArrayBox& vel_fab   = velo_func.getFArrayBox();
    const RealVect dx          = m_amr->getDx()[a_lvl]*RealVect::Unit;
    const RealVect origin      = m_amr->getProbLo();
    const Box box              = m_amr->getGrids(m_realm)[a_lvl][a_dit];

    List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();

    // This interpolates the velocity function on to the particle velocities
    EbParticleInterp meshInterp(box, ebisbox, dx, origin);
    meshInterp.interpolate<ItoParticle, &ItoParticle::velocity>(particleList, vel_fab, m_deposition, m_irreg_ngp_interpolation);

    // Go through the particles and set their velocities to velo_func*mobility
    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit){
      ItoParticle& p = lit();
      p.velocity() *= p.mobility();
    }
  }
}

void ItoSolver::interpolateMobilities(){
  CH_TIME("ItoSolver::interpolateMobilities()");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolateMobilities()" << endl;
  }

  if(m_isMobile){

    if(m_WhichMobilityInterpolation == WhichMobilityInterpolation::velocity){
      DataOps::vectorLength(m_scratch, m_velo_func); // Compute |E| (or whatever other function you've decided to provide).
      m_amr->averageDown(m_scratch, m_realm, m_phase);
      m_amr->interpGhost(m_scratch, m_realm, m_phase);
    }
      
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	this->interpolateMobilities(lvl, dit());
      }
    }
  }
}

void ItoSolver::interpolateMobilities(const int a_lvl, const DataIndex& a_dit){
  CH_TIME("ItoSolver::interpolateMobilities(lvl, dit)");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolateMobilities(lvl, dit)" << endl;
  }


  if(m_WhichMobilityInterpolation == WhichMobilityInterpolation::mobility){
    this->interpolateMobilitiesMu(a_lvl, a_dit);
  }
  else if (m_WhichMobilityInterpolation == WhichMobilityInterpolation::velocity){
    this->interpolateMobilitiesVel(a_lvl, a_dit);
  }

}

void ItoSolver::interpolateMobilitiesMu(const int a_lvl, const DataIndex& a_dit){
  CH_TIME("ItoSolver::interpolateMobilitiesMu(lvl, dit)");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolateMobilitiesMu(lvl, dit)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_ParticleContainers.at(WhichContainer::bulk);

  if(m_isMobile){
    const EBCellFAB& mob_func  = (*m_mobility_func[a_lvl])[a_dit];
    const EBISBox& ebisbox     = mob_func.getEBISBox();
    const FArrayBox& mob_fab   = mob_func.getFArrayBox();
    const RealVect dx          = m_amr->getDx()[a_lvl]*RealVect::Unit;
    const RealVect origin      = m_amr->getProbLo();
    const Box box              = m_amr->getGrids(m_realm)[a_lvl][a_dit];

    List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();
    EbParticleInterp meshInterp(box, ebisbox, dx, origin);
    
    meshInterp.interpolate<ItoParticle, &ItoParticle::mobility>(particleList, mob_fab, m_deposition, m_irreg_ngp_interpolation);
  }
}

void ItoSolver::interpolateMobilitiesVel(const int a_lvl, const DataIndex& a_dit){
  CH_TIME("ItoSolver::interpolateMobilitiesVel(lvl, dit)");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolateMobilitiesVel(lvl, dit)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_ParticleContainers.at(WhichContainer::bulk);

  if(m_isMobile){
    const EBCellFAB& mob_func  = (*m_mobility_func[a_lvl])[a_dit];
    const EBISBox& ebisbox     = mob_func.getEBISBox();
    const FArrayBox& mob_fab   = mob_func.getFArrayBox();
    const RealVect dx          = m_amr->getDx()[a_lvl]*RealVect::Unit;
    const RealVect origin      = m_amr->getProbLo();
    const Box box              = m_amr->getGrids(m_realm)[a_lvl][a_dit];

    FArrayBox& scratch = (*m_scratch[a_lvl])[a_dit].getFArrayBox();

    List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();
    EbParticleInterp meshInterp(box, ebisbox, dx, origin);
    
    // First, interpolate |E| to the particle position, it will be stored on m_tmp. 
    meshInterp.interpolate<ItoParticle, &ItoParticle::mobility>(particleList, scratch, m_deposition, m_irreg_ngp_interpolation);    

    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit){
      ItoParticle& p = lit();
      
      p.tmp() = p.mobility();
    }

    // This interpolates mu*|E| to the particle position and stores it on the mobility. After that, we compute mu_p = (mu*E)/E
    scratch *= mob_fab;
    meshInterp.interpolate<ItoParticle, &ItoParticle::mobility>(particleList, scratch, m_deposition, m_irreg_ngp_interpolation);    
    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit){
      ItoParticle& p = lit();

      p.mobility() *= 1./p.tmp();
    }
  }
}

void ItoSolver::updateMobilities(){
  CH_TIME("ItoSolver::updateMobilities");
  if(m_verbosity > 5){
    pout() << m_name + "::updateMobilities" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      this->updateMobilities(lvl, dit());
    }
  }
}

void ItoSolver::updateMobilities(const int a_level, const DataIndex a_dit){
  CH_TIME("ItoSolver::updateMobilities(lvl, dit)");
  if(m_verbosity > 5){
    pout() << m_name + "::updateMobilities(lvl, dit)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_ParticleContainers.at(WhichContainer::bulk);

  if(m_isMobile){
    List<ItoParticle>& particleList = particles[a_level][a_dit].listItems();

    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit){
      ItoParticle& p = lit();
      
      p.mobility() = m_species->mobility(p.energy());
    }
  }
}

void ItoSolver::interpolateDiffusion(){
  CH_TIME("ItoSolver::interpolateDiffusion");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolateDiffusion" << endl;
  }

  if(m_isDiffusive){
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	this->interpolateDiffusion(lvl, dit());
      }
    }
  }
}

void ItoSolver::interpolateDiffusion(const int a_lvl, const DataIndex& a_dit){
  CH_TIME("ItoSolver::interpolateDiffusion");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolateDiffusion" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_ParticleContainers.at(WhichContainer::bulk);

  if(m_isDiffusive){
    const EBCellFAB& dco_cell  = (*m_faceCenteredDiffusionCoefficient_cell[a_lvl])[a_dit];
    const EBISBox& ebisbox     = dco_cell.getEBISBox();
    const FArrayBox& dco_fab   = dco_cell.getFArrayBox();
    const RealVect dx          = m_amr->getDx()[a_lvl]*RealVect::Unit;
    const RealVect origin      = m_amr->getProbLo();
    const Box box              = m_amr->getGrids(m_realm)[a_lvl][a_dit];

    List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();

    EbParticleInterp meshInterp(box, ebisbox,dx, origin);
    meshInterp.interpolate<ItoParticle, &ItoParticle::diffusion>(particleList, dco_fab, m_deposition, m_irreg_ngp_interpolation);    
  }
}

void ItoSolver::updateDiffusion(){
  CH_TIME("ItoSolver::updateDiffusion");
  if(m_verbosity > 5){
    pout() << m_name + "::updateDiffusion" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      this->updateDiffusion(lvl, dit());
    }
  }
}

void ItoSolver::updateDiffusion(const int a_level, const DataIndex a_dit){
  CH_TIME("ItoSolver::updateDiffusion(lvl, dit)");
  if(m_verbosity > 5){
    pout() << m_name + "::updateDiffusion(lvl, dit)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::bulk);

  if(m_isDiffusive){
    List<ItoParticle>& particleList = particles[a_level][a_dit].listItems();

    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit){
      ItoParticle& p = lit();
      
      p.mobility() = m_species->diffusion(p.energy());
    }
  }
}



Real ItoSolver::computeDt() const {
  CH_TIME("ItoSolver::computeDt(allAMRlevels)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDt(allAMRlevels)" << endl;
  }

  Real dt = 1.E99;
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const Real levelDt = this->computeDt(lvl);
    dt = Min(dt, levelDt);
  }

  return dt;
}

Real ItoSolver::computeDt(const int a_lvl) const{
  CH_TIME("ItoSolver::computeDt(lvl)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDt(lvl)" << endl;
  }

  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const Real dx = m_amr->getDx()[a_lvl];
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Real boxDt = this->computeDt(a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("ItoSolver::computeDt(lvl) - communication error on norm");
  }
  dt = tmp;
#endif  

  return dt;
}
  

Real ItoSolver::computeDt(const int a_lvl, const DataIndex a_dit, const Real a_dx) const{
  CH_TIME("ItoSolver::computeDt(lvl, dit, dx)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDt(lvl, dit, dx)" << endl;
  }

  Real dt = 1.E99;

  const ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::bulk);

  const List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();
  
  ListIterator<ItoParticle> lit(particleList);
  if(m_isMobile && !m_isDiffusive){
    for (lit.rewind(); lit.ok(); ++lit){
      const ItoParticle& p = particleList[lit];
      const RealVect& v = p.velocity();

      const int maxDir  = v.maxDir(true);
      const Real vMax   = Abs(v[maxDir]);
      const Real thisDt = (vMax > 0.0) ? a_dx/vMax : 1.E99;

      dt = Min(dt, thisDt);
    }
  }
  else if(!m_isMobile && m_isDiffusive){
    for (lit.rewind(); lit.ok(); ++lit){
      const ItoParticle& p = particleList[lit];

      const Real D      = p.diffusion();
      const Real thisDt = (D > 0.0) ? a_dx*a_dx/(2.0*D) : 1.E99;

      dt = Min(dt, thisDt);
    }
  }
  else if(m_isMobile && m_isDiffusive){
    for (lit.rewind(); lit.ok(); ++lit){
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
  if(m_verbosity > 5){
    pout() << m_name + "::computeMinDt(allAMRlevels, maxCellsToMove)" << endl;
  }

  Real dt = 1.E99;
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const Real levelDt = this->computeMinDt(a_maxCellsToMove, lvl);
    dt = Min(dt, levelDt);
  }

  return dt;
}

Real ItoSolver::computeMinDt(const Real a_maxCellsToMove, const int a_lvl) const{
  CH_TIME("ItoSolver::computeMinDt(maxCellsToMove, lvl)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeMinDt(maxCellsToMove, lvl)" << endl;
  }

  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const Real dx = m_amr->getDx()[a_lvl];
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Real boxDt = this->computeMinDt(a_maxCellsToMove, a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("ItoSolver::computeDt(lvl) - communication error on norm");
  }
  dt = tmp;
#endif  

  return dt;
}

Real ItoSolver::computeMinDt(const Real a_maxCellsToMove, const int a_lvl, const DataIndex a_dit, const Real a_dx) const{
  CH_TIME("ItoSolver::computeMinDt(maxCellsToMove, lvl, dit, dx)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeMinDt(maxCellsToMove, lvl, dit, dx)" << endl;
  }

  Real dt = 1.E99;

  const ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::bulk);

  const Real dMax  = a_maxCellsToMove*a_dx;
  const Real dMax2 = dMax*dMax;
  const Real W0    = m_normal_max;
  const Real W02   = m_normal_max*m_normal_max;

  const List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();
  
  ListIterator<ItoParticle> lit(particleList);
  if(m_isMobile && !m_isDiffusive){
    for (lit.rewind(); lit; ++lit){
      const ItoParticle& p = particleList[lit];
      const RealVect& v = p.velocity();

      const int maxDir = v.maxDir(true);
      const Real thisDt = dMax/Abs(v[maxDir]);

      dt = Min(dt, thisDt);
    }
  }
  else if(!m_isMobile && m_isDiffusive){
    for (lit.rewind(); lit; ++lit){
      const ItoParticle& p = particleList[lit];
      
      const Real thisDt = dMax2/(2.0*p.diffusion()*W02);
      dt = Min(dt, thisDt);
    }
  }
  else if(m_isMobile && m_isDiffusive){
    for (lit.rewind(); lit; ++lit){
      const ItoParticle& p = particleList[lit];
      
      const RealVect& v = p.velocity();
      const int maxDir = v.maxDir(true);
      const Real vMax  = Abs(v[maxDir]);
      const Real D     = p.diffusion();

      
      if(vMax > 0.0){
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
	if(D > 0.0){
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
  if(m_verbosity > 5){
    pout() << m_name + "::computeMinDriftDt(allAMRlevels, maxCellsToMove)" << endl;
  }

  Vector<Real> dt = this->computeDriftDt(a_maxCellsToMove);

  Real minDt = dt[0];
  for (int lvl = 0; lvl < dt.size(); lvl++){
    minDt = Min(minDt, dt[lvl]);
  }

  return minDt;
}

Vector<Real> ItoSolver::computeDriftDt(const Real a_maxCellsToMove) const {
  CH_TIME("ItoSolver::computeDriftDt(amr, maxCellsToMove)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDriftDt(amr, maxCellsToMove)" << endl;
  }

  Vector<Real> dt = this->computeDriftDt();
  for (int lvl = 0; lvl < dt.size(); lvl++){
    dt[lvl] = dt[lvl]*a_maxCellsToMove;
  }

  return dt;

}

Real ItoSolver::computeAdvectiveDt() const {
  CH_TIME("ItoSolver::computeAdvectiveDt");
  if(m_verbosity > 5){
    pout() << m_name + "::computeAdvectiveDt()" << endl;
  }
  
  Real minDt = 1.E99;
  const Vector<Real> levelDts = this->computeDriftDt();

  for (int lvl = 0; lvl < levelDts.size(); lvl++){
    minDt = Min(levelDts[lvl], minDt);
  }

  return minDt;
}

Vector<Real> ItoSolver::computeDriftDt() const {
  CH_TIME("ItoSolver::computeDriftDt(amr)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDriftDt(amr)" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  Vector<Real> dt(1 + finest_level, 1.2345E67);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    dt[lvl] = this->computeDriftDt(lvl);
  }

  return dt;
}

Real ItoSolver::computeDriftDt(const int a_lvl) const {
  CH_TIME("ItoSolver::computeDriftDt(level)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDriftDt(level)" << endl;
  }

  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const RealVect dx = m_amr->getDx()[a_lvl]*RealVect::Unit;
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Real boxDt = this->computeDriftDt(a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("ItoSolver::compute_drift_level(lvl) - communication error on norm");
  }
  dt = tmp;
#endif  

  return dt;
}

Real ItoSolver::computeDriftDt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const{
  CH_TIME("ItoSolver::computeDriftDt(level, dataindex, dx)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDriftDt(level, dataindex, dx)" << endl;
  }

  const ParticleContainer<ItoParticle>& particles = m_ParticleContainers.at(WhichContainer::bulk);

  constexpr Real safety = 1.E-10;

  const List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();

  Real dt = 1.E99;

  if(m_isMobile){
    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit){
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
  if(m_verbosity > 5){
    pout() << m_name + "::computeMinDiffusionDt(min, maxCellsToHop)" << endl;
  }

  Vector<Real> dt = this->computeDiffusionDt(a_maxCellsToHop);
  Real minDt = dt[0];
  for (int lvl = 0; lvl < dt.size(); lvl++){
    minDt = Min(minDt, dt[lvl]);
  }

  return minDt;
}

Vector<Real> ItoSolver::computeDiffusionDt(const Real a_maxCellsToHop) const{
  CH_TIME("ItoSolver::computeMinDiffusionDt(maxCellsToHop)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeMinDiffusionDt(maxCellsToHop)" << endl;
  }

  const Real factor  = a_maxCellsToHop/m_normal_max;
  const Real factor2 = factor*factor;
  
  Vector<Real> dt = this->computeDiffusionDt();
  for (int lvl = 0; lvl < dt.size(); lvl++){
    dt[lvl] = dt[lvl]*factor2;
  }

  return dt;
}

Real ItoSolver::computeDiffusiveDt() const {
  CH_TIME("ItoSolver::computeDiffusiveDt");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDiffusiveDt" << endl;
  }

  Real minDt = 1.E99;
  const Vector<Real> levelDts = this->computeDiffusionDt();

  for (int lvl = 0; lvl < levelDts.size(); lvl++){
    minDt = Min(levelDts[lvl], minDt);
  }

  return minDt;
}

Vector<Real> ItoSolver::computeDiffusionDt() const{
  CH_TIME("ItoSolver::computeDiffusionDt");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDiffusionDt" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();
    
  Vector<Real> dt(1 + finest_level, 1.2345E6);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    dt[lvl] = this->computeDiffusionDt(lvl);
  }

  return dt;
}

Real ItoSolver::computeDiffusionDt(const int a_lvl) const{
  CH_TIME("ItoSolver::computeDiffusionDt(level)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDiffusionDt(level)" << endl;
  }

  
  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const RealVect dx = m_amr->getDx()[a_lvl]*RealVect::Unit;
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Real boxDt = this->computeDiffusionDt(a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("ItoSolver::compute_diffusion_level(lvl) - communication error on norm");
  }
  dt = tmp;
#endif

  return dt;
}

Real ItoSolver::computeDiffusionDt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const{
  CH_TIME("ItoSolver::computeDiffusionDt(level, dataindex, dx)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDiffusionDt(level, dataindex, dx)" << endl;
  }

  const ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::bulk);
  
  const List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();

  Real dt = 1.E99;

  if(m_isDiffusive){
    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit){
      const ItoParticle& p = particleList[lit];
    
      const Real thisDt = a_dx[0]*a_dx[0]/(2.0*p.diffusion());
    
      dt = Min(dt, thisDt);
    }
  }

  return dt;
}

void ItoSolver::remap(){
  CH_TIME("ItoSolver::remap");
  if(m_verbosity > 5){
    pout() << m_name + "::remap" << endl;
  }

  this->remap(WhichContainer::bulk);

}

void ItoSolver::remap(const WhichContainer a_container){
  CH_TIME("ItoSolver::remap(container)");
  if(m_verbosity > 5){
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

void ItoSolver::sortParticlesByCell(const WhichContainer a_container){
  CH_TIME("ItoSolver::sortParticlesByCell(container)");
  if(m_verbosity > 5){
    pout() << m_name + "::sortParticlesByCell(container)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);

  particles.sortParticlesByCell();
}


void ItoSolver::sortParticlesByPatch(const WhichContainer a_container){
  CH_TIME("ItoSolver::sortParticlesByPatch(container)");
  if(m_verbosity > 5){
    pout() << m_name + "::sortParticlesByPatch(container)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);

  particles.sortParticlesByPatch();
}

void ItoSolver::makeSuperparticles(const WhichContainer a_container, const int a_particlesPerPatch){
  CH_TIME("ItoSolver::makeSuperparticles(int)");
  if(m_verbosity > 5){
    pout() << m_name + "::makeSuperparticles(int)" << endl;
  }
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    this->makeSuperparticles(a_container, a_particlesPerPatch, lvl);
  }

}

void ItoSolver::makeSuperparticles(const WhichContainer a_container, const int a_particlesPerPatch, const int a_level){
  CH_TIME("ItoSolver::makeSuperparticles(int, level)");
  if(m_verbosity > 5){
    pout() << m_name + "::makeSuperparticles(int, level)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::bulk);

  const DisjointBoxLayout& dbl = particles.getGrids()[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    this->makeSuperparticles(a_container, a_particlesPerPatch, a_level, dit());
  }
}

void ItoSolver::makeSuperparticles(const WhichContainer a_container, const int a_particlesPerCell, const int a_level, const DataIndex a_dit){
  CH_TIME("ItoSolver::makeSuperparticles(int, level, patch)");
  if(m_verbosity > 5){
    pout() << m_name + "::makeSuperparticles(int, level, patch)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);
  
  const int comp = 0;
  //  const Box box  = m_amr->getGrids(m_realm)[a_level].get(a_dit);

  const Box box  = particles.getGrids()[a_level][a_dit];

  // This are the particles in the box we're currently looking at. 
  BinFab<ItoParticle>& cellParticles = particles.getCellParticles(a_level, a_dit);

  // Iterate over particles
  for (BoxIterator bit(box); bit.ok(); ++bit){
    const IntVect iv = bit();
  
    List<ItoParticle>& particles = cellParticles(iv, comp);

    if(particles.length() > 0){
      this->mergeBVH(particles, a_particlesPerCell);
    }
  }
}

void ItoSolver::mergeBVH(List<ItoParticle>& a_particles, const int a_particlesPerCell){
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
  for (ListIterator<ItoParticle> lit(a_particles); lit.ok(); ++lit, i++){
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
  const int dir = (m_kd_direction < 0) ? m_udist0d(m_rng) : m_kd_direction;
  m_tree.define(pointMasses, mass);
  m_tree.buildTree(dir, a_particlesPerCell);
  const std::vector<std::shared_ptr<ItoMerge::Node<PointMass> > >& leaves = m_tree.getLeaves();

  // 3. Clear particles in this cell and add new ones.
  a_particles.clear();
  for (int i = 0; i < leaves.size(); i++){
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

void ItoSolver::clear(const WhichContainer a_container){
  CH_TIME("ItoSolver::clear(string)");
  if(m_verbosity > 5){
    pout() << m_name + "::clear(string)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);

  this->clear(particles);
}

void ItoSolver::clear(ParticleContainer<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::clear(ParticleContainer)");
  if(m_verbosity > 5){
    pout() << m_name + "::clear(ParticleContainer)" << endl;
  }

  this->clear(a_particles.getParticles());
}

void ItoSolver::clear(AMRParticles<ItoParticle>& a_particles){
  CH_TIME("ItoSolver::clear");
  if(m_verbosity > 5){
    pout() << m_name + "::clear" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    a_particles[lvl]->clear();
  }
}

RealVect ItoSolver::randomPosition(const RealVect a_pos,
				   const RealVect a_lo,
				   const RealVect a_hi,
				   const RealVect a_bndryCentroid,
				   const RealVect a_bndryNormal,
				   const Real     a_dx,
				   const Real     a_kappa) {

  RealVect pos;
  if(a_kappa < 1.0){ // Rejection sampling. 
    pos = this->randomPosition(a_lo, a_hi, a_bndryCentroid, a_bndryNormal);
  }
  else{ // Regular cell. Get a position. 
    pos = this->randomPosition(a_lo, a_hi);
  }

  pos = a_pos + pos*a_dx;

  return pos;
}

RealVect ItoSolver::randomPosition(const RealVect a_lo,
				   const RealVect a_hi,
				   const RealVect a_bndryCentroid,
				   const RealVect a_bndryNormal) {
  RealVect pos = this->randomPosition(a_lo, a_hi);
  bool valid   = PolyGeom::dot(pos-a_bndryCentroid, a_bndryNormal) >= 0.0;

  while(!valid){
    pos    = this->randomPosition(a_lo, a_hi);
    valid = PolyGeom::dot(pos-a_bndryCentroid, a_bndryNormal) >= 0.0;
  }

  return pos;
}

RealVect ItoSolver::randomPosition(const RealVect a_lo, const RealVect a_hi) {

  RealVect pos = RealVect::Unit;

  for (int dir = 0; dir < SpaceDim; dir++){
    pos[dir] = a_lo[dir] + 0.5*(1.0 + m_udist11(m_rng))*(a_hi[dir] - a_lo[dir]);
  }

  return pos;
}

#include <CD_NamespaceFooter.H>
