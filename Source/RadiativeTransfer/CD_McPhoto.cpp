/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_McPhoto.cpp
  @brief  Implementation of CD_McPhoto.H
  @author Robert Marskar
*/

// Std includes
#include <time.h>
#include <chrono>

// Chombo includes
#include <PolyGeom.H>
#include <EBAlias.H>
#include <EBLevelDataOps.H>
#include <BoxIterator.H>
#include <EBArith.H>
#include <ParmParse.H>
#include <ParticleIO.H>

// Our includes
#include <CD_McPhoto.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_PolyUtils.H>
#include <CD_ParticleOps.H>
#include <CD_NamespaceHeader.H>

#define MC_PHOTO_DEBUG 0

McPhoto::McPhoto(){
  m_name      = "McPhoto";
  m_className = "McPhoto";

  m_stationary = false;
}

McPhoto::~McPhoto(){

}

bool McPhoto::advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zerophi){
  CH_TIME("McPhoto::advance");
  if(m_verbosity > 5){
    pout() << m_name + "::advance" << endl;
  }

  // Note: This routine does an on-the-fly generation of Photons based on the contents in a_source. Pure particle methods
  //       will probably fill m_sourcePhotons themselves, and this routine will probably not be used with a pure particle
  //       method. If you find yourself calling this routine with a pure particle method, you're probably something wrong. 

  // If stationary, do a safety cleanout first. Then generate new Photons
  if(m_instantaneous){
    this->clear(m_photons);
  }

  // Generate new Photons and add them to m_photons
  this->clear(m_sourcePhotons);
  this->generatePhotons(m_sourcePhotons, a_source, a_dt);
  m_photons.addParticles(m_sourcePhotons);

  // Advance stationary or transient
  if(m_instantaneous){
    this->advancePhotonsStationary(m_bulkPhotons, m_ebPhotons, m_domainPhotons, m_photons);
  }
  else{
    this->advancePhotonsTransient(m_bulkPhotons, m_ebPhotons, m_domainPhotons, m_photons, a_dt);
    this->remap(m_photons);                                       
  }

  // Deposit volumetric Photons. 
  this->depositPhotons(a_phi, m_bulkPhotons, m_deposition);
  
  return true;
}

bool McPhoto::isInstantaneous(){
  return m_instantaneous;
}

void McPhoto::parseOptions(){
  CH_TIME("McPhoto::parseOptions");
  if(m_verbosity > 5){
    pout() << m_name + "::parseOptions" << endl;
  }
  
  this->parseRng();
  this->parsePseudoPhotons();
  this->parsePhotoGeneration();
  this->parseSourceType();
  this->parseDeposition();
  this->parseBisectStep();
  this->parseDomainBc();
  this->parsePvrBuffer();
  this->parsePlotVariables();
  this->parseInstantaneous();
  this->parseDivergenceComputation();
}

void McPhoto::parseRuntimeOptions(){
  CH_TIME("McPhoto::parseRuntimeOptions");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }
  
  this->parsePseudoPhotons();
  this->parsePhotoGeneration();
  this->parseSourceType();
  this->parseDeposition();
  this->parseBisectStep();
  this->parseDomainBc();
  this->parsePvrBuffer();
  this->parsePlotVariables();
  this->parseInstantaneous();
  this->parseDivergenceComputation();
}

void McPhoto::parseDivergenceComputation(){
  CH_TIME("McPhoto::parseDivergenceComputation");
  if(m_verbosity > 5){
    pout() << m_name + "::parseDivergenceComputation" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("blend_conservation", m_blendConservation);
}

void McPhoto::parseRng(){
  CH_TIME("McPhoto::parseRng");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRng" << endl;
  }

  // Seed the RNG
  ParmParse pp(m_className.c_str());
  pp.get("seed", m_seed);
  pp.get("poiss_exp_swap", m_poiss_exp_swap);
  if(m_seed < 0) {
    m_seed = std::chrono::system_clock::now().time_since_epoch().count();
    m_seed += procID(); 
  }
  m_rng = new std::mt19937_64(m_seed);

  m_udist01 = new uniform_real_distribution<Real>( 0.0, 1.0);
  m_udist11 = new uniform_real_distribution<Real>(-1.0, 1.0);
}

void McPhoto::parseInstantaneous(){
  CH_TIME("McPhoto::parseInstantaneous");
  if(m_verbosity > 5){
    pout() << m_name + "::parseInstantaneous" << endl;
  }
  
  ParmParse pp(m_className.c_str());

  pp.get("instantaneous", m_instantaneous);
}

void McPhoto::parsePseudoPhotons(){
  CH_TIME("McPhoto::parsePseudoPhotons");
  if(m_verbosity > 5){
    pout() << m_name + "::parsePseudoPhotons" << endl;
  }
  
  ParmParse pp(m_className.c_str());
  
  pp.get("max_photons", m_max_photons);
  if(m_max_photons <= 0){ // = -1 => no restriction
    m_max_photons = 99999999;
  }
}

void McPhoto::parsePhotoGeneration(){
  CH_TIME("McPhoto::parsePhotoGeneration");
  if(m_verbosity > 5){
    pout() << m_name + "::parsePhotoGeneration" << endl;
  }
  
  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("photon_generation", str);

  if(str == "deterministic"){
    m_photogen = PhotonGeneration::Deterministic;
  }
  else if(str == "stochastic"){
    m_photogen = PhotonGeneration::Stochastic;
  }
  else{
    MayDay::Abort("McPhoto::set_PhotonGeneration - unknown Photon generation type requested");
  }
}

void McPhoto::parseSourceType(){
  CH_TIME("McPhoto::parseSourceType");
  if(m_verbosity > 5){
    pout() << m_name + "::parseSourceType" << endl;
  }
  
  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("source_type", str);

  if(str == "number"){
    m_src_type = SourceType::Number;
  }
  else if(str == "volume"){
    m_src_type = SourceType::PerVol;
  }
  else if(str == "volume_rate"){
    m_src_type = SourceType::PerVolSecond;
  }
  else if(str == "rate"){
    m_src_type = SourceType::PerSecond;
  }
  else{
    MayDay::Abort("McPhoto::setSourceType - unknown source type requested");
  }
}

void McPhoto::parseDeposition(){
  CH_TIME("McPhoto::parseDeposition");
  if(m_verbosity > 5){
    pout() << m_name + "::parseDeposition" << endl;
  }
  
  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("deposition", str);

  m_deposit_numbers = false;
  if(str == "num"){
    m_deposition = DepositionType::NGP;
    m_deposit_numbers = true;
  }
  else if(str == "ngp"){
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
    MayDay::Abort("McPhoto::set_deposition_type - unknown interpolant requested");
  }

  pp.get("plot_deposition", str);
  m_plotNumbers = false;
  if(str == "num"){
    m_plot_deposition = DepositionType::NGP;
    m_plotNumbers = true;
  }
  else if(str == "ngp"){
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
    MayDay::Abort("McPhoto::set_deposition_type - unknown interpolant requested");
  }
}

void McPhoto::parseBisectStep(){
  CH_TIME("McPhoto::parseBisectStep");
  if(m_verbosity > 5){
    pout() << m_name + "::parseBisectStep" << endl;
  }
  
  ParmParse pp(m_className.c_str());
  pp.get("bisect_step", m_bisect_step);
}

void McPhoto::parseDomainBc(){
  CH_TIME("McPhoto::parseDomainBc");
  if(m_verbosity > 5){
    pout() << m_name + "::parseDomainBc" << endl;
  }
  
  m_domainbc.resize(2*SpaceDim);
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const Side::LoHiSide side = sit();
      const int idx = domainBcMap(dir, side);

      ParmParse pp(m_className.c_str());
      std::string str_dir;
      if(dir == 0){
	str_dir = "x";
      }
      else if(dir == 1){
	str_dir = "y";
      }
      else if(dir == 2){
	str_dir = "z";
      }

      std::string sidestr = (side == Side::Lo) ? "_low" : "_high";
      std::string bc_string = "bc_" + str_dir + sidestr;
      std::string type;
      pp.get(bc_string.c_str(), type);
      if(type == "outflow"){
	m_domainbc[idx] = wallbc::outflow;
      }
      else if(type == "symmetry"){
	m_domainbc[idx] = wallbc::symmetry;
      }
      else if(type == "wall"){
	m_domainbc[idx] = wallbc::wall;
      }
      else {
	std::string error = "McPhoto::setDomainBc - unsupported boundary condition requested: " + bc_string;
	MayDay::Abort(error.c_str());
      }
    }
  }
}

void McPhoto::parsePvrBuffer(){
  CH_TIME("McPhoto::parsePvrBuffer");
  if(m_verbosity > 5){
    pout() << m_name + "::parsePvrBuffer" << endl;
  }
  
  ParmParse pp(m_className.c_str());
  pp.get("pvr_buffer",   m_pvr_buffer);
  pp.get("halo_buffer",  m_halo_buffer);
}

void McPhoto::parsePlotVariables(){
  CH_TIME("McPhoto::parsePlotVariables");
  if(m_verbosity > 5){
    pout() << m_name + "::parsePlotVariables" << endl;
  }

  m_plotPhi       = false;
  m_plotSource       = false;
  m_plot_phot      = false;
  m_plot_bulk_phot = false;
  m_plot_eb_phot   = false;
  m_plot_dom_phot  = false;
  m_plot_source_phot  = false;

  ParmParse pp(m_className.c_str());
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++){
    if(     str[i] == "phi")       m_plotPhi       = true;
    else if(str[i] == "src")       m_plotSource       = true;
    else if(str[i] == "phot")      m_plot_phot      = true;
    else if(str[i] == "bulk_phot") m_plot_bulk_phot = true;
    else if(str[i] == "eb_phot")   m_plot_eb_phot   = true;
    else if(str[i] == "dom_phot")  m_plot_dom_phot  = true;
    else if(str[i] == "src_phot")  m_plot_source_phot  = true;
  }
}

void McPhoto::clear(){
  CH_TIME("McPhoto::clear()");
  if(m_verbosity > 5){
    pout() << m_name + "::clear()" << endl;
  }

  this->clear(m_photons);
}

void McPhoto::clear(ParticleContainer<Photon>& a_photons){
  CH_TIME("McPhoto::clear(ParticleContainer)");
  if(m_verbosity > 5){
    pout() << m_name + "::clear(ParticleContainer)" << endl;
  }

  this->clear(a_photons.getParticles());
}

void McPhoto::clear(AMRParticles<Photon>& a_photons){
  CH_TIME("McPhoto::clear(AMRParticles)");
  if(m_verbosity > 5){
    pout() << m_name + "::clear(AMRParticles)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    a_photons[lvl]->clear();
  }
}
  
void McPhoto::allocateInternals(){
  CH_TIME("McPhoto::allocateInternals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocateInternals" << endl;
  }

  const int ncomp  = 1;

  // Allocate mesh data
  m_amr->allocate(m_phi,        m_realm, m_phase, ncomp); 
  m_amr->allocate(m_source,       m_realm, m_phase, ncomp);
  m_amr->allocate(m_scratch,      m_realm, m_phase, ncomp);
  m_amr->allocate(m_depositionNC, m_realm, m_phase, ncomp);
  m_amr->allocate(m_massDiff,     m_realm, m_phase, ncomp);

  // Allocate particle data holders
  m_amr->allocate(m_photons,        m_pvr_buffer, m_realm);
  m_amr->allocate(m_bulkPhotons,   m_pvr_buffer, m_realm);
  m_amr->allocate(m_ebPhotons,     m_pvr_buffer, m_realm);
  m_amr->allocate(m_domainPhotons, m_pvr_buffer, m_realm);
  m_amr->allocate(m_sourcePhotons, m_pvr_buffer, m_realm);
}

void McPhoto::preRegrid(const int a_lmin, const int a_oldFinestLevel){
  CH_TIME("McPhoto::preRegrid");
  if(m_verbosity > 5){
    pout() << m_name + "::pre_grid" << endl;
  }

  m_photons.preRegrid(a_lmin);         // TLDR: This moves Photons from l >= a_lmin to Max(a_lmin-1,0)
  m_bulkPhotons.preRegrid(a_lmin);    // TLDR: This moves Photons from l >= a_lmin to Max(a_lmin-1,0)
  m_ebPhotons.preRegrid(a_lmin);      // TLDR: This moves Photons from l >= a_lmin to Max(a_lmin-1,0)
  m_domainPhotons.preRegrid(a_lmin);  // TLDR: This moves Photons from l >= a_lmin to Max(a_lmin-1,0)
  m_sourcePhotons.preRegrid(a_lmin);  // TLDR: This moves Photons from l >= a_lmin to Max(a_lmin-1,0)
}

void McPhoto::deallocateInternals(){
  CH_TIME("McPhoto::deallocateInternals");
  if(m_verbosity > 5){ 
    pout() << m_name + "::deallocateInternals" << endl;
  }

  // Don't deallocate, instead, reallocate. 
  // m_amr->deallocate(m_phi);
  // m_amr->deallocate(m_source);
  // m_amr->deallocate(m_scratch);
  // m_amr->deallocate(m_depositionNC);
  // m_amr->deallocate(m_massDiff);
}

void McPhoto::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("McPhoto::regrid");
  if(m_verbosity > 5){ 
    pout() << m_name + "::regrid" << endl;
  }

  // Mesh data regrids
  m_amr->reallocate(m_phi,        m_phase, a_lmin); 
  m_amr->reallocate(m_source,       m_phase, a_lmin);
  m_amr->reallocate(m_scratch,      m_phase, a_lmin);
  m_amr->reallocate(m_depositionNC, m_phase, a_lmin);
  m_amr->reallocate(m_massDiff,     m_phase, a_lmin);

  // Particle data regrids
  const Vector<DisjointBoxLayout>& grids = m_amr->getGrids(m_realm);
  const Vector<ProblemDomain>& domains   = m_amr->getDomains();
  const Vector<Real>& dx                 = m_amr->getDx();
  const Vector<int>& ref_rat             = m_amr->getRefinementRatios();

  m_photons.regrid(       grids, domains, dx, ref_rat, a_lmin, a_newFinestLevel);
  m_bulkPhotons.regrid(  grids, domains, dx, ref_rat, a_lmin, a_newFinestLevel);
  m_ebPhotons.regrid(    grids, domains, dx, ref_rat, a_lmin, a_newFinestLevel);
  m_domainPhotons.regrid(grids, domains, dx, ref_rat, a_lmin, a_newFinestLevel);
  m_sourcePhotons.regrid(grids, domains, dx, ref_rat, a_lmin, a_newFinestLevel);

  // Deposit
  this->depositPhotons();
}

void McPhoto::sortPhotonsByCell(){
  CH_TIME("McPhoto::sortPhotonsByCell()");
  if(m_verbosity > 5){
    pout() << m_name + "::sortPhotonsByCell()" << endl;
  }

  m_photons.sortParticlesByCell();
}

void McPhoto::sortPhotonsByPatch(){
  CH_TIME("McPhoto::sortPhotonsByPatch()");
  if(m_verbosity > 5){
    pout() << m_name + "::sortPhotonsByPatch()" << endl;
  }

  m_photons.sortParticlesByPatch();
}

void McPhoto::sortBulkPhotonsByCell(){
  CH_TIME("McPhoto::sortBulkPhotonsByCell()");
  if(m_verbosity > 5){
    pout() << m_name + "::sortBulkPhotonsByCell()" << endl;
  }

  m_bulkPhotons.sortParticlesByCell();
}

void McPhoto::sortBulkPhotonsByPatch(){
  CH_TIME("McPhoto::sortBulkPhotonsByPatch()");
  if(m_verbosity > 5){
    pout() << m_name + "::sortBulkPhotonsByPatch()" << endl;
  }

  m_bulkPhotons.sortParticlesByPatch();
}

void McPhoto::sortSourcePhotonsByCell(){
  CH_TIME("McPhoto::sortSourcePhotonsByCell()");
  if(m_verbosity > 5){
    pout() << m_name + "::sortSourcePhotonsByCell()" << endl;
  }

  m_sourcePhotons.sortParticlesByCell();
}

void McPhoto::sortSourcePhotonsByPatch(){
  CH_TIME("McPhoto::sortSourcePhotonsByPatch()");
  if(m_verbosity > 5){
    pout() << m_name + "::sortSourcePhotonsByPatch()" << endl;
  }

  m_sourcePhotons.sortParticlesByPatch();
}

void McPhoto::sortDomainPhotonsByCell(){
  CH_TIME("McPhoto::sortDomainPhotonsByCell()");
  if(m_verbosity > 5){
    pout() << m_name + "::sortDomainPhotonsByCell()" << endl;
  }

  m_domainPhotons.sortParticlesByCell();
}

void McPhoto::sortDomainPhotonsByPatch(){
  CH_TIME("McPhoto::sortDomainPhotonsByPatch()");
  if(m_verbosity > 5){
    pout() << m_name + "::sortDomainPhotonsByPatch()" << endl;
  }

  m_domainPhotons.sortParticlesByPatch();
}

void McPhoto::sortEbPhotonsByCell(){
  CH_TIME("McPhoto::sortEbPhotonsByCell()");
  if(m_verbosity > 5){
    pout() << m_name + "::sortEbPhotonsByCell()" << endl;
  }

  m_ebPhotons.sortParticlesByCell();
}

void McPhoto::sortEbPhotonsByPatch(){
  CH_TIME("McPhoto::sortEbPhotonsByPatch()");
  if(m_verbosity > 5){
    pout() << m_name + "::sortEbPhotonsByPatch()" << endl;
  }

  m_ebPhotons.sortParticlesByPatch();
}

void McPhoto::registerOperators(){
  CH_TIME("McPhoto::registerOperators");
  if(m_verbosity > 5){
    pout() << m_name + "::registerOperators" << endl;
  }

  if(m_amr.isNull()){
    MayDay::Abort("McPhoto::registerOperators - need to set AmrMesh!");
  }
  else{
    m_amr->registerOperator(s_eb_coar_ave,     m_realm, m_phase);
    m_amr->registerOperator(s_eb_fill_patch,   m_realm, m_phase);
    m_amr->registerOperator(s_eb_mg_interp,    m_realm, m_phase);
    m_amr->registerOperator(s_eb_redist,       m_realm, m_phase);
    m_amr->registerOperator(s_noncons_div,  m_realm, m_phase);
    m_amr->registerOperator(s_eb_copier,       m_realm, m_phase);
    if(m_pvr_buffer <= 0){
      m_amr->registerOperator(s_eb_ghostcloud, m_realm, m_phase);
    }
    if(m_halo_buffer > 0){
      m_amr->registerMask(s_particle_halo, m_halo_buffer, m_realm);
    }
  }
}

void McPhoto::computeBoundaryFlux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi){
  CH_TIME("McPhoto::computeBoundaryFlux");
  if(m_verbosity > 5){
    pout() << m_name + "::computeBoundaryFlux" << endl;
  }

  DataOps::setValue(a_ebFlux, 0.0);
}

void McPhoto::computeDomainFlux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_phi){
  CH_TIME("McPhoto::computeDomainFlux");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDomainFlux" << endl;
  }
  
  DataOps::setValue(a_domainflux, 0.0);
}

void McPhoto::computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_phi){
  const std::string str = "McPhoto::computeFlux - Fluid flux can't be computed with discrete Photons. Calling this is an error";
  MayDay::Abort(str.c_str());
}

void McPhoto::computeDensity(EBAMRCellData& a_isotropic, const EBAMRCellData& a_phi){
  MayDay::Abort("McPhoto::computeDensity - Calling this is an error");
}

void McPhoto::writePlotFile(){
  CH_TIME("McPhoto::writePlotFile");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotFile" << endl;
  }

  MayDay::Abort("McPhoto::writePlotFile - not implemented");
}

void McPhoto::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("McPhoto::writeCheckpointLevel");
  if(m_verbosity > 5){
    pout() << m_name + "::writeCheckpointLevel" << endl;
  }

  // Write state vector
  write(a_handle, *m_phi[a_level], m_name);

  // Write particles. Must be implemented.
  std::string str = m_name + "_particles";
  writeParticlesToHDF(a_handle, m_photons[a_level], str);
}

void McPhoto::readCheckpointLevel(HDF5Handle& a_handle, const int a_level){
  CH_TIME("McPhoto::readCheckpointLevel");
  if(m_verbosity > 5){
    pout() << m_name + "::readCheckpointLevel" << endl;
  }

  // Read state vector
  read<EBCellFAB>(a_handle, *m_phi[a_level], m_name, m_amr->getGrids(m_realm)[a_level], Interval(0,0), false);

  // Read particles. Should be implemented
  std::string str = m_name + "_particles";
  readParticlesFromHDF(a_handle, m_photons[a_level], str);
}

Vector<std::string> McPhoto::getPlotVariableNames() const {
  CH_TIME("McPhoto::getPlotVariableNames");
  if(m_verbosity > 5){
    pout() << m_name + "::getPlotVariableNames" << endl;
  }
  
  Vector<std::string> names(0);
  
  if(m_plotPhi)       names.push_back(m_name + " phi");
  if(m_plotSource)       names.push_back(m_name + " source");
  if(m_plot_phot)      names.push_back(m_name + " Photons");
  if(m_plot_bulk_phot) names.push_back(m_name + " bulkPhotons");
  if(m_plot_eb_phot)   names.push_back(m_name + " ebPhotons");
  if(m_plot_dom_phot)  names.push_back(m_name + " domainPhotons");
  if(m_plot_source_phot)  names.push_back(m_name + " sourcePhotons");

  return names;
}

int McPhoto::getNumberOfPlotVariables() const{
  CH_TIME("McPhoto::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }

  int num_output = 0;

  if(m_plotPhi)       num_output = num_output + 1;
  if(m_plotSource)       num_output = num_output + 1;
  if(m_plot_phot)      num_output = num_output + 1;
  if(m_plot_bulk_phot) num_output = num_output + 1;
  if(m_plot_eb_phot)   num_output = num_output + 1;
  if(m_plot_dom_phot)  num_output = num_output + 1;
  if(m_plot_source_phot)  num_output = num_output + 1;

  return num_output;
}

int McPhoto::queryGhost() const {
  return 3;
}

int McPhoto::randomPoisson(const Real a_mean){
  if(a_mean < m_poiss_exp_swap){
    std::poisson_distribution<int> pdist(a_mean);
    return pdist(*m_rng);
  }
  else {
    std::normal_distribution<Real> ndist(a_mean, sqrt(a_mean));
    return (int) round(ndist(*m_rng));
  }
}

int McPhoto::domainBcMap(const int a_dir, const Side::LoHiSide a_side) {
  const int iside = (a_side == Side::Lo) ? 0 : 1;

  return 2*a_dir + iside;
}

Real McPhoto::randomExponential(const Real a_mean){
  std::exponential_distribution<Real> dist(a_mean);
  return dist(*m_rng);
}

RealVect McPhoto::randomDirection(){
#if CH_SPACEDIM == 2
  return randomDirection2D();
#else
  return randomDirection3D();
#endif
}

#if CH_SPACEDIM == 2
RealVect McPhoto::randomDirection2D(){
  const Real EPS = 1.E-8;
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = (*m_udist11)(*m_rng);
    x2 = (*m_udist11)(*m_rng);
    r  = x1*x1 + x2*x2;
  }

  return RealVect(x1,x2)/sqrt(r);
}
#endif

#if CH_SPACEDIM==3
RealVect McPhoto::randomDirection3D(){
  const Real EPS = 1.E-8;
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = (*m_udist11)(*m_rng);
    x2 = (*m_udist11)(*m_rng);
    r  = x1*x1 + x2*x2;
  }

  const Real x = 2*x1*sqrt(1-r);
  const Real y = 2*x2*sqrt(1-r);
  const Real z = 1 - 2*r;

  return RealVect(x,y,z);
}
#endif

int McPhoto::getPVRBuffer() const {
  CH_TIME("McPhoto::getPVRBuffer");
  if(m_verbosity > 5){
    pout() << m_name + "::getPVRBuffer" << endl;
  }

  return m_pvr_buffer;
}

int McPhoto::getHaloBuffer() const {
  CH_TIME("McPhoto::getHaloBuffer");
  if(m_verbosity > 5){
    pout() << m_name + "::getHaloBuffer" << endl;
  }

  return m_halo_buffer;
}

void McPhoto::setPVRBuffer(const int a_buffer) {
  CH_TIME("McPhoto::setPVRBuffer");
  if(m_verbosity > 5){
    pout() << m_name + "::setPVRBuffer" << endl;
  }

  m_pvr_buffer = a_buffer;
}

void McPhoto::setHalobuffer(const int a_buffer)  {
  CH_TIME("McPhoto::setHalobuffer");
  if(m_verbosity > 5){
    pout() << m_name + "::setHalobuffer" << endl;
  }

  m_halo_buffer = a_buffer;
}

void McPhoto::generatePhotons(ParticleContainer<Photon>& a_photons, const EBAMRCellData& a_source, const Real a_dt){
  CH_TIME("McPhoto::generatePhotons");
  if(m_verbosity > 5){
    pout() << m_name + "::generatePhotons" << endl;
  }

  const RealVect prob_lo = m_amr->getProbLo();
  const int finest_level = m_amr->getFinestLevel();

  // Put here. 
  AMRParticles<Photon>& Photons = a_photons.getParticles();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& dom     = m_amr->getDomains()[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    const Real vol               = pow(dx, SpaceDim);
    const bool has_coar          = (lvl > 0);

    // If there is a coarser level, remove particles from the overlapping region and regenerate on this level.
    if(has_coar) {
      const AMRPVR& pvr = a_photons.getPVR();
      const int ref_ratio = m_amr->getRefinementRatios()[lvl-1];
      collectValidParticles(Photons[lvl]->outcast(),
			    *Photons[lvl-1],
			    pvr[lvl]->mask(),
			    dx*RealVect::Unit,
			    ref_ratio,
			    false,
			    prob_lo);
      Photons[lvl]->outcast().clear(); 
    }

    // Create new particles on this level. 
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = (*a_source[lvl])[dit()].getEBISBox();
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs   = IntVectSet(box);

      const EBCellFAB& source = (*a_source[lvl])[dit()];
      const FArrayBox& srcFAB = source.getFArrayBox();

      Real sum = srcFAB.sum(0);
      
      // Generate new particles in this box
      List<Photon> particles;
      if(sum > 0){
	
	// Regular cells
	for (BoxIterator bit(box); bit.ok(); ++bit){
	  const IntVect iv   = bit();

	  if(ebisbox.isRegular(iv)){
	    const RealVect pos = prob_lo + (RealVect(iv) + 0.5*RealVect::Unit)*dx;

	    // Number of physical Photons
	    const int num_phys_photons = this->drawPhotons(srcFAB(iv,0), vol, a_dt);

	    // Make superPhotons if we have to
	    if(num_phys_photons > 0){
	      const int num_photons = (num_phys_photons <= m_max_photons) ? num_phys_photons : m_max_photons;
	      const Real weight     = (1.0*num_phys_photons)/num_photons; 
	      
	      for (int i = 0; i < num_photons; i++){
		const RealVect v = Units::c*this->randomDirection();
		particles.append(Photon(pos, v, m_RtSpecies->getKappa(pos), weight));
	      }
	    }
	  }
	}

	// Irregular and multicells
	VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
	for (vofit.reset(); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();
	  const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, prob_lo);
	  const Real kappa    = ebisbox.volFrac(vof);

	  const int num_phys_photons = drawPhotons(source(vof,0), vol, a_dt);
	  
	  if(num_phys_photons > 0){
	    const int num_photons = (num_phys_photons <= m_max_photons) ? num_phys_photons : m_max_photons;
	    const Real weight      = (1.0*num_phys_photons)/num_photons;

	    // Generate computational Photons 
	    for (int i = 0; i < num_photons; i++){
	      const RealVect v = Units::c*this->randomDirection();
	      particles.append(Photon(pos, v, m_RtSpecies->getKappa(pos), weight));
	    }
	  }

	}
      }

      // Add new particles to data holder
      (*Photons[lvl])[dit()].addItemsDestructive(particles);
    }
  }
}

int McPhoto::drawPhotons(const Real a_source, const Real a_volume, const Real a_dt){
  CH_TIME("McPhoto::drawPhotons");
  if(m_verbosity > 5){
    pout() << m_name + "::drawPhotons" << endl;
  }

  int num_photons = 0;

  // Check if we need any type of source term normalization
  Real factor;
  if(m_src_type == SourceType::Number){
    factor = 1.0;
  }
  else if(m_src_type == SourceType::PerVol){
    factor = a_volume;
  }
  else if(m_src_type == SourceType::PerVolSecond){
    factor = a_volume*a_dt;
  }
  else if(m_src_type == SourceType::PerSecond){
    factor = a_dt;
  }

  // Draw a number of Photons with the desired algorithm
  if(m_photogen == PhotonGeneration::Stochastic){
    const Real mean = a_source*factor;
    num_photons = randomPoisson(mean);
  }
  else if(m_photogen == PhotonGeneration::Deterministic){
    num_photons = round(a_source*factor);
  }
  else{
    MayDay::Abort("mc::drawPhotons - unknown generation requested. Aborting...");
  }

  return num_photons;
}

void McPhoto::depositPhotons(){
  CH_TIME("McPhoto::depositPhotons");
  if(m_verbosity > 5){
    pout() << m_name + "::depositPhotons" << endl;
  }

  this->depositPhotons(m_phi, m_photons, m_deposition);
}

void McPhoto::depositPhotons(EBAMRCellData&                    a_phi,
			       const ParticleContainer<Photon>& a_photons,
			       const DepositionType::Which&      a_deposition){
  CH_TIME("McPhoto::depositPhotons(ParticleContainer)");
  if(m_verbosity > 5){
    pout() << m_name + "::depositPhotons(ParticleContainer)" << endl;
  }

  this->depositPhotons(a_phi, a_photons.getParticles(), a_deposition);
}

void McPhoto::depositPhotons(EBAMRCellData&               a_phi,
			       const AMRParticles<Photon>&  a_photons,
			       const DepositionType::Which& a_deposition){
  CH_TIME("McPhoto::depositPhotons(AMRParticles)");
  if(m_verbosity > 5){
    pout() << m_name + "::depositPhotons(AMRParticles)" << endl;
  }

           
  this->depositKappaConservative(a_phi, a_photons, a_deposition); // a_phi contains only weights, i.e. not divided by kappa
  this->depositNonConservative(m_depositionNC, a_phi);              // Compute m_depositionNC = sum(kappa*Wc)/sum(kappa)
  this->depositHybrid(a_phi, m_massDiff, m_depositionNC);           // Compute hybrid deposition, including mass differnce
  this->incrementRedist(m_massDiff);                                 // Increment level redistribution register

  // Do the redistribution magic
  const bool ebcf = m_amr->getEbCf();
  if(ebcf){ // Mucho stuff to do here...
    this->coarseFineIncrement(m_massDiff);       // Compute C2F, F2C, and C2C mass transfers
    this->levelRedist(a_phi);           // Level redistribution. Weights is a dummy parameter
    this->coarseFineRedistribution(a_phi);     // Do the coarse-fine redistribution
  }
  else{ // Very simple, redistribute this level.
    this->levelRedist(a_phi);
  }

  // Average down and interpolate
  m_amr->averageDown(a_phi, m_realm, m_phase);
  m_amr->interpGhost(a_phi, m_realm, m_phase);
}

void McPhoto::depositKappaConservative(EBAMRCellData&              a_phi,
					 const AMRParticles<Photon>& a_photons,
					 const DepositionType::Which a_deposition){
  CH_TIME("McPhoto::depositKappaConservative");
  if(m_verbosity > 5){
    pout() << m_name + "::depositKappaConservative" << endl;
  }

  const int comp = 0;
  const Interval interv(comp, comp);

  const RealVect origin  = m_amr->getProbLo();
  const int finest_level = m_amr->getFinestLevel();

  DataOps::setValue(a_phi,    0.0);
  DataOps::setValue(m_scratch,  0.0);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->getDx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& dom     = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const RefCountedPtr<EBLevelGrid>& eblg = m_amr->getEBLevelGrid(m_realm, m_phase)[lvl];

    const bool has_coar = (lvl > 0);
    const bool has_fine = (lvl < finest_level);

    // 1. If we have a coarser level whose cloud extends beneath this level, interpolate that result here first. 
    if(has_coar && m_pvr_buffer > 0){
      RefCountedPtr<EBMGInterp>& interp = m_amr->getEBMGInterp(m_realm, m_phase)[lvl];
      interp->pwcInterp(*a_phi[lvl], *m_scratch[lvl-1], interv);
    }
    
    // 2. Deposit this levels Photons. Note that this will deposit into ghost cells, which must later
    //    be added to neighboring patches
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      EbParticleInterp interp(box, ebisbox, dx*RealVect::Unit, origin, true);
      interp.deposit((*a_photons[lvl])[dit()].listItems(), (*a_phi[lvl])[dit()].getFArrayBox(), m_deposition);
    }

    // This code adds contributions from ghost cells into the valid region
    const RefCountedPtr<Copier>& reversecopier = m_amr->getReverseCopier(m_realm, m_phase)[lvl];
    LDaddOp<FArrayBox> addOp;
    LevelData<FArrayBox> aliasFAB;
    aliasEB(aliasFAB, *a_phi[lvl]);
    aliasFAB.exchange(Interval(0,0), *reversecopier, addOp);

    // 3. If we have a finer level, copy contributions from this level to the temporary holder that is used for
    //    interpolation of "hanging clouds"
    if(has_fine){
      a_phi[lvl]->localCopyTo(*m_scratch[lvl]);
    }
    else if(m_pvr_buffer <= 0 && has_coar){
      EbGhostCloud& ghostcloud = *(m_amr->getGhostCloud(m_realm, m_phase)[lvl]);
      ghostcloud.addFineGhostsToCoarse(*a_phi[lvl-1], *a_phi[lvl]);
    }
  }
}

void McPhoto::depositNonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC){
  CH_TIME("McPhoto::depositNonConservative");
  if(m_verbosity > 5){
    pout() << m_name + "::depositNonConservative" << endl;
  }

  if(m_blendConservation){
    const IrregAmrStencil<NonConservativeDivergenceStencil>& stencils = m_amr->getNonConservativeDivergenceStencils(m_realm, m_phase);
    stencils.apply(a_depositionNC, a_depositionKappaC);
  }
  else{
    DataOps::setValue(a_depositionNC, 0.0);
  }
}

void McPhoto::depositHybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_massDifference, const EBAMRIVData& a_depositionNC){
  CH_TIME("McPhoto::depositHybrid");
  if(m_verbosity > 5){
    pout() << m_name + "::depositHybrid" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      EBCellFAB& divH               = (*a_depositionH[lvl])[dit()];  // On input, this contains kappa*depositionWeights
      BaseIVFAB<Real>& deltaM       = (*a_massDifference[lvl])[dit()];
      const BaseIVFAB<Real>& divNC  = (*a_depositionNC[lvl])[dit()]; 

      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs   = ebisbox.getIrregIVS(box);

      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const Real kappa    = ebisbox.volFrac(vof);
	const Real dc       = divH(vof, comp);
	const Real dnc      = divNC(vof, comp);

	// Note that if dc - kappa*dnc can be negative, i.e. we may end up STEALING mass
	// from other cells. This is why there is a flag m_blendConservation which always
	// gives positive definite results. 
	divH(vof, comp)   = dc + (1-kappa)*dnc;          // On output, contains hybrid divergence
	deltaM(vof, comp) = (1-kappa)*(dc - kappa*dnc);  // Remember, dc already scaled by kappa. 
      }
    }
  }
}

void McPhoto::incrementRedist(const EBAMRIVData& a_massDifference){
  CH_TIME("McPhoto::incrementRedist");
  if(m_verbosity > 5){
    pout() << m_name + "::incrementRedist" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    
    EBLevelRedist& level_redist = *(m_amr->getLevelRedist(m_realm, m_phase)[lvl]);
    level_redist.setToZero();

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      level_redist.increment((*a_massDifference[lvl])[dit()], dit(), interv);
    }
  }
}

void McPhoto::levelRedist(EBAMRCellData& a_phi){
  CH_TIME("McPhoto::levelRedist");
  if(m_verbosity > 5){
    pout() << m_name + "::levelRedist" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    EBLevelRedist& level_redist = *(m_amr->getLevelRedist(m_realm, m_phase)[lvl]);
    level_redist.redistribute(*a_phi[lvl], interv);
    level_redist.setToZero();
  }
}

void McPhoto::coarseFineIncrement(const EBAMRIVData& a_massDifference){
  CH_TIME("McPhoto::coarseFineIncrement");
  if(m_verbosity > 5){
    pout() << m_name + "::coarseFineIncrement" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(0,0);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->getFineToCoarRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->getCoarToFineRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->getCoarToCoarRedist(m_realm, m_phase)[lvl];

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

void McPhoto::coarseFineRedistribution(EBAMRCellData& a_phi){
  CH_TIME("McPhoto::coarseFineRedistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::coarseFineRedistribution" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx       = m_amr->getDx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->getCoarToFineRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->getCoarToCoarRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->getFineToCoarRedist(m_realm, m_phase)[lvl];
    
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

void McPhoto::advancePhotonsStationary(ParticleContainer<Photon>& a_bulkPhotons,
					  ParticleContainer<Photon>& a_ebPhotons,
					  ParticleContainer<Photon>& a_domainPhotons,
					  ParticleContainer<Photon>& a_photons){
  CH_TIME("McPhoto::advancePhotonsStationary");
  if(m_verbosity > 5){
    pout() << m_name + "::advancePhotonsStationary" << endl;
  }

  // TLDR: This routine iterates over the levels and boxes and does the following
  //
  //       Forall Photons in a_photons: {
  //          1. Draw random absorption position
  //          2. Check if path intersects boundary, either EB or domain
  //          3. Move the Photon to appropriate data holder:
  //                 Path crossed EB   => a_ebPhotons
  //                 Path cross domain => a_domainPhotons
  //                 Absorbe in bulk   => a_bulkPhotons
  //       }
  //
  //       Remap a_bulkPhotons, a_ebPhotons, a_domainPhotons

  // Low and high corners
  const RealVect prob_lo = m_amr->getProbLo();
  const RealVect prob_hi = m_amr->getProbHi();

  // Safety factor
  const Real SAFETY = 1.E-6;

  // This is the implicit function used for intersection tests
  RefCountedPtr<BaseIF> impfunc;
  if(m_phase == phase::gas){
    impfunc = m_computationalGeometry->getGasImplicitFunction();
  }
  else{
    impfunc = m_computationalGeometry->getSolidImplicitFunction();
  }

#if MC_PHOTO_DEBUG // Debug
  int PhotonsBefore = this->countPhotons(a_photons.getParticles());
#endif
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<Photon>& bulkPhotons = a_bulkPhotons[lvl][dit()].listItems();
      List<Photon>& ebPhotons   = a_ebPhotons[lvl][dit()].listItems();
      List<Photon>& domPhotons  = a_domainPhotons[lvl][dit()].listItems();
      List<Photon>& allPhotons  = a_photons[lvl][dit()].listItems();

      // These must be cleared
      bulkPhotons.clear();
      ebPhotons.clear();
      domPhotons.clear();

      // Iterate over the Photons that will be moved. 
      for (ListIterator<Photon> lit(allPhotons); lit.ok(); ++lit){
	Photon& p = lit();

	// Draw a new random absorption position
	const RealVect oldPos  = p.position();
	const RealVect unitV   = p.velocity()/(p.velocity().vectorLength());
	const RealVect newPos  = oldPos + unitV*this->randomExponential(p.kappa());
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
	  if(newPos[dir] < prob_lo[dir] || newPos[dir] > prob_hi[dir]){
	    checkDom = true; 
	  }
	}


	if(!checkEB && !checkDom){ // No intersection test necessary, add to bulk absorption
	  p.position() = newPos;
	  bulkPhotons.add(p);
	}
	else{ // Must do nasty intersect test.
	  Real dom_s = 1.E99;
	  Real eb_s  = 1.E99;

	  bool contact_domain = false;
	  bool contact_eb     = false;

	  // Do intersection tests
	  if(checkDom) contact_domain = ParticleOps::domainIntersection(oldPos, newPos, path, prob_lo, prob_hi, dom_s);
	  //	  if(checkEB)  contact_eb = ParticleOps::ebIntersectionBisect(impfunc, oldPos, newPos, pathLen, m_bisect_step, eb_s);
	  if(checkEB)  contact_eb = ParticleOps::ebIntersectionRaycast(impfunc, oldPos, newPos, 1.E-10*dx, eb_s);

	  // Move the Photon to the data holder where it belongs. 
	  if(!contact_eb && !contact_domain){
	    p.position() = newPos;
	    bulkPhotons.add(p);
	  }
	  else {
	    if(eb_s < dom_s){
	      p.position() = oldPos + eb_s*path;
	      ebPhotons.add(p);
	    }
	    else{
	      p.position() = oldPos + Max(0.0,dom_s-SAFETY)*path;
	      domPhotons.add(p);
	    }
	  }
	}
      }

      // Should clear these out. 
      allPhotons.clear();
    }
  }

  // Remap and clear.
  a_bulkPhotons.remap();
  a_ebPhotons.remap();
  a_domainPhotons.remap();

#if MC_PHOTO_DEBUG // Debug
  int bulkPhotons = this->countPhotons(a_bulkPhotons.getParticles());
  int ebPhotons = this->countPhotons(a_ebPhotons.getParticles());
  int domPhotons = this->countPhotons(a_domainPhotons.getParticles());

  if(procID() == 0){
    std::cout << "Photons before = " << PhotonsBefore << "\n"
	      << "bulk Photons = " << bulkPhotons << "\n"
	      << "eb Photons = " << ebPhotons << "\n"
	      << "dom Photons = " << domPhotons << "\n"
	      << "Photons after = " << domPhotons+ebPhotons+bulkPhotons << "\n" << std::endl;
  }
#endif


}

void McPhoto::advancePhotonsTransient(ParticleContainer<Photon>& a_bulkPhotons,
					 ParticleContainer<Photon>& a_ebPhotons,
					 ParticleContainer<Photon>& a_domainPhotons,
					 ParticleContainer<Photon>& a_photons,
					 const Real                  a_dt){
  CH_TIME("McPhoto::advancePhotonsTransient");
  if(m_verbosity > 5){
    pout() << m_name + "::advancePhotonsTransient" << endl;
  }

  // TLDR: This routine iterates over the levels and boxes and does the following
  //
  //       Forall Photons in a_photons: {
  //          1. Check new Photon position
  //          2. Check if the Photon was absorbed on the interval
  //          3. Check if path intersects boundary, either EB or domain
  //          4. Move the Photon to appropriate data holder:
  //                 Path crossed EB   => a_ebPhotons
  //                 Path cross domain => a_domainPhotons
  //                 Absorbed in bulk  => a_bulkPhotons
  //       }
  //
  //       Remap a_bulkPhotons, a_ebPhotons, a_domainPhotons, a_photons

  // Low and high corners
  const RealVect prob_lo = m_amr->getProbLo();
  const RealVect prob_hi = m_amr->getProbHi();

  // Safety factor
  const Real SAFETY = 1.E-8;

  // This is the implicit function used for intersection tests
  RefCountedPtr<BaseIF> impfunc;
  if(m_phase == phase::gas){
    impfunc = m_computationalGeometry->getGasImplicitFunction();
  }
  else{
    impfunc = m_computationalGeometry->getSolidImplicitFunction();
  }

#if MC_PHOTO_DEBUG // Debug
  int PhotonsBefore = this->countPhotons(a_photons.getParticles());
#endif

  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const Real dx = m_amr->getDx()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<Photon>& bulkPhotons = a_bulkPhotons[lvl][dit()].listItems();
      List<Photon>& ebPhotons   = a_ebPhotons[lvl][dit()].listItems();
      List<Photon>& domPhotons  = a_domainPhotons[lvl][dit()].listItems();
      List<Photon>& allPhotons  = a_photons[lvl][dit()].listItems();

      
      // These must be cleared
      bulkPhotons.clear();
      ebPhotons.clear();
      domPhotons.clear();
      

      // Iterate over the Photons that will be moved. 
      for (ListIterator<Photon> lit(allPhotons); lit.ok(); ++lit){
	Photon& p = lit();

	// Move the Photon
	const RealVect oldPos  = p.position();
	const RealVect v       = p.velocity();
	const RealVect unitV   = v/v.vectorLength();
	const RealVect newPos  = oldPos + v*a_dt;
	const RealVect path    = newPos - oldPos;
	const Real     pathLen = path.vectorLength();

	// Check if we should check of different types of boundary intersections. These are checp initial tests that allow
	// us to skip more expensive intersection tests for some Photons.
	bool checkEB  = false;
	bool checkDom = false;

	if(!impfunc.isNull()){
	  checkEB = true;
	}
	for (int dir = 0; dir < SpaceDim; dir++){
	  if(newPos[dir] <= prob_lo[dir] || newPos[dir] >= prob_hi[dir]){
	    checkDom = true; 
	  }
	}

	// Handles for how long the Photons move
	bool absorbed_bulk   = false;
	bool absorbed_eb     = false;
	bool absorbed_domain = false;

	Real bulk_s = 1.E99;
	Real eb_s   = 1.E99;
	Real dom_s  = 1.E99;


	// Check absorption in the bulk
	const Real travelLen  = this->randomExponential(p.kappa());
	if(travelLen < pathLen){
	  absorbed_bulk = true;
	  bulk_s        = travelLen/pathLen;
	}

	// Check absorption on EBs and domain
	if(checkEB){
	  //	  absorbed_eb = ParticleOps::ebIntersectionBisect(impfunc, oldPos, newPos, pathLen, m_bisect_step, eb_s);
	  absorbed_eb = ParticleOps::ebIntersectionRaycast(impfunc, oldPos, newPos, 1.E-10*dx, eb_s);
	}
	if(checkDom){
	  absorbed_domain = ParticleOps::domainIntersection(oldPos, newPos, path, prob_lo, prob_hi, dom_s);
	  dom_s = (absorbed_domain) ? Max(0.0, dom_s - SAFETY) : dom_s;
	}

	const bool absorb = absorbed_bulk || absorbed_eb || absorbed_domain;
	if(!absorb){
	  p.position() = newPos;
	}
	else{
	  const Real s = Min(bulk_s, Min(eb_s, dom_s));
	  p.position() = oldPos + s*path;

	  // Now check where it was actually absorbed
	  if(absorbed_bulk && bulk_s < Min(eb_s, dom_s)){ 
	    bulkPhotons.transfer(lit);
	  }
	  else if(absorbed_eb && eb_s < Min(bulk_s, dom_s)){
	    ebPhotons.transfer(lit);
	  }
	  else if(absorbed_domain && dom_s < Min(bulk_s, eb_s)){
	    domPhotons.transfer(lit);
	  }
	  else{
	    MayDay::Abort("McPhoto::advancePhotonsTransient - logic bust");
	  }
	}
      }
    }
  }

  // Remap and clear. 
  a_bulkPhotons.remap();
  a_ebPhotons.remap();
  a_domainPhotons.remap();
  a_photons.remap();

#if MC_PHOTO_DEBUG // Debug
  int bulkPhotons = this->countPhotons(a_bulkPhotons.getParticles());
  int ebPhotons = this->countPhotons(a_ebPhotons.getParticles());
  int domPhotons = this->countPhotons(a_domainPhotons.getParticles());
  int afterPhotons = this->countPhotons(a_photons.getParticles());

  if(procID() == 0){
    std::cout << "Photons before = " << PhotonsBefore << "\n"
	      << "Photons after = " << afterPhotons << "\n"
	      << "bulk Photons = " << bulkPhotons << "\n"
	      << "eb Photons = " << ebPhotons << "\n"
	      << "dom Photons = " << domPhotons << "\n"
	      << "total = " << domPhotons+ebPhotons+bulkPhotons+afterPhotons << "\n" << std::endl;
  }


#endif
  a_bulkPhotons.remap();
}

void McPhoto::remap(){
  CH_TIME("McPhoto::remap()");
  if(m_verbosity > 5){
    pout() << m_name + "::remap()" << endl;
  }

  this->remap(m_photons);
}

void McPhoto::remap(ParticleContainer<Photon>& a_photons){
  CH_TIME("McPhoto::remap(Photons)");
  if(m_verbosity > 5){
    pout() << m_name + "::remap(Photons)" << endl;
  }

  a_photons.remap();
}

int McPhoto::countPhotons(const AMRParticles<Photon>& a_photons) const {
  CH_TIME("McPhoto::countPhotons");
  if(m_verbosity > 5){
    pout() << m_name + "::countPhotons" << endl;
  }

  int num_photons = 0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    num_photons += a_photons[lvl]->numValid();
  }

  return num_photons;
}

int McPhoto::countOutcast(const AMRParticles<Photon>& a_photons) const {
  CH_TIME("McPhoto::countOutcast");
  if(m_verbosity > 5){
    pout() << m_name + "::countOutcast" << endl;
  }

  int num_outcast = 0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    num_outcast += a_photons[lvl]->numOutcast();
  }

  return num_outcast;
}

void McPhoto::writePlotData(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("McPhoto::writePlotData");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotData" << endl;
  }

  if(m_plotPhi) {
    this->writeData(a_output, a_comp, m_phi,  false);
  }
  if(m_plotSource) {
    this->writeData(a_output, a_comp, m_source, false);
  }
  if(m_plot_phot){
    this->depositPhotons(m_scratch, m_photons.getParticles(), m_plot_deposition);
    this->writeData(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_bulk_phot){
    this->depositPhotons(m_scratch, m_bulkPhotons.getParticles(), m_plot_deposition);
    this->writeData(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_eb_phot){
    this->depositPhotons(m_scratch, m_ebPhotons.getParticles(), m_plot_deposition);
    this->writeData(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_dom_phot){
    this->depositPhotons(m_scratch, m_domainPhotons.getParticles(), m_plot_deposition);
    this->writeData(a_output, a_comp, m_scratch,  false);
  }
  if(m_plot_source_phot){
    this->depositPhotons(m_scratch, m_sourcePhotons.getParticles(), m_plot_deposition);
    this->writeData(a_output, a_comp, m_scratch,  false);
  }
}

ParticleContainer<Photon>& McPhoto::getPhotons(){
  CH_TIME("McPhoto::getPhotons");
  if(m_verbosity > 5){
    pout() << m_name + "::getPhotons" << endl;
  }

  return m_photons;
}

ParticleContainer<Photon>& McPhoto::getBulkPhotons(){
  CH_TIME("McPhoto::getBulkPhotons");
  if(m_verbosity > 5){
    pout() << m_name + "::getBulkPhotons" << endl;
  }

  return m_bulkPhotons;
}

ParticleContainer<Photon>& McPhoto::getEbPhotons(){
  CH_TIME("McPhoto::getEbPhotons");
  if(m_verbosity > 5){
    pout() << m_name + "::getEbPhotons" << endl;
  }

  return m_ebPhotons;
}

ParticleContainer<Photon>& McPhoto::getDomainPhotons(){
  CH_TIME("McPhoto::getDomainPhotons");
  if(m_verbosity > 5){
    pout() << m_name + "::getDomainPhotons" << endl;
  }

  return m_domainPhotons;
}

ParticleContainer<Photon>& McPhoto::getSourcePhotons(){
  CH_TIME("McPhoto::getSourcePhotons");
  if(m_verbosity > 5){
    pout() << m_name + "::getSourcePhotons" << endl;
  }

  return m_sourcePhotons;
}

#include <CD_NamespaceFooter.H>
