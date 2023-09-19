/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_McPhoto.cpp
  @brief  Implementation of CD_McPhoto.H
  @author Robert Marskar
  @todo   Create new classes for handling domain BCs -- the old ones stink. 
*/

// Std includes
#include <time.h>
#include <chrono>

// Chombo includes
#include <ParmParse.H>
#include <ParticleIO.H>

// Our includes
#include <CD_ParticleManagement.H>
#include <CD_Location.H>
#include <CD_McPhoto.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_ParticleOps.H>
#include <CD_Random.H>
#include <CD_NamespaceHeader.H>

McPhoto::McPhoto()
{
  CH_TIME("McPhoto::McPhoto");

  m_name      = "McPhoto";
  m_className = "McPhoto";

  m_stationary = false;
}

McPhoto::~McPhoto() {}

bool
McPhoto::advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zerophi)
{
  CH_TIME("McPhoto::advance");
  if (m_verbosity > 5) {
    pout() << m_name + "::advance" << endl;
  }

  // Note: This routine does an on-the-fly generation of photons based on the contents in a_source. This routine is primarily
  //       written for fluid methods since all input/output parameters happen on the mesh. Most particle solvers will use a different
  //       approach and will probably fill m_sourcePhotons through a kinetic interface and use advancePhotonsInstantaneous or advancePhotonsTransient.
  //
  //       If you find yourself calling this routine with a pure particle method, you're may be doing something you shouldn't....
  //
  //       This routine is a bit convoluted because it permits sub-sampling of the photons generated in the grid cells. That is, instead of generating
  //       all photons at once and moving them, we divide them into packets and advance the packets independently of one another. This results in additional
  //       MPI calls, but can reduce memory when many photons are generated per cell.

  DataOps::setValue(a_phi, 0.0);

  EBAMRCellData phi;
  EBAMRCellData numPhysPhotonsTotal;
  EBAMRCellData numPhysPhotonsPacket;

  m_amr->allocate(phi, m_realm, m_phase, 1);
  m_amr->allocate(numPhysPhotonsTotal, m_realm, m_phase, 1);
  m_amr->allocate(numPhysPhotonsPacket, m_realm, m_phase, m_numSamplingPackets);

  // Recall: m_photons are 'traveling' photons, m_sourcePhotons are photons to be transferred into m_photons before transport over dt,
  //         and m_ebPhotons and m_domainPhotons are photons that strike the EB or domain boundaries.
  //
  //         If we're doing an instantaneous solver, do a safety cleanout of m_photons first (past photons have been absorbed on the mesh).
  if (m_instantaneous) {
    this->clear(m_photons);
  }

  this->clear(m_bulkPhotons);
  this->clear(m_sourcePhotons);
  this->clear(m_ebPhotons);
  this->clear(m_domainPhotons);

  this->computeNumPhysicalPhotons(numPhysPhotonsTotal, numPhysPhotonsPacket, a_source, a_dt);

  // Advance photons either instantaneously or transiently. For transient advances
  if (m_instantaneous) {
    const size_t maxPhotonsPerPacket = m_maxPhotonsGeneratedPerCell / m_numSamplingPackets;
    const size_t remainder           = m_maxPhotonsGeneratedPerCell % m_numSamplingPackets;

    ParticleContainer<Photon> scratchPhotons;
    m_amr->allocate(scratchPhotons, m_realm);

    for (int i = 0; i < m_numSamplingPackets; i++) {
      const size_t maxPhotonsPerCell = (i == 0) ? maxPhotonsPerPacket + remainder : maxPhotonsPerPacket;

      const EBAMRCellData& numPhysPhotons = m_amr->slice(numPhysPhotonsPacket, Interval(i, i));

      this->generateComputationalPhotons(m_photons, numPhysPhotons, maxPhotonsPerCell);
      this->advancePhotonsInstantaneous(scratchPhotons, m_ebPhotons, m_domainPhotons, m_photons);

      // Absorb the bulk photons on the mesh.
      this->depositPhotons(phi, scratchPhotons, m_deposition);
      DataOps::incr(a_phi, phi, 1.0);

      // Store the photons that were absorbed.
      m_bulkPhotons.transferParticles(scratchPhotons);
      scratchPhotons.clearParticles();
    }
  }
  else {
    this->generateComputationalPhotons(m_photons, numPhysPhotonsTotal, m_maxPhotonsGeneratedPerCell);
    this->advancePhotonsTransient(m_bulkPhotons, m_ebPhotons, m_domainPhotons, m_photons, a_dt);
    this->remap(m_photons);
    this->depositPhotons(a_phi, m_bulkPhotons, m_deposition);
  }

  m_amr->conservativeAverage(a_phi, m_realm, m_phase);
  m_amr->interpGhost(a_phi, m_realm, m_phase);

  return true;
}

bool
McPhoto::isInstantaneous()
{
  return m_instantaneous;
}

void
McPhoto::parseOptions()
{
  CH_TIME("McPhoto::parseOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseOptions" << endl;
  }

  // For registering the CIC deposition mask.
  m_haloBuffer = 1;

  this->parseVerbosity();
  this->parseTransparentBoundaries();
  this->parsePseudoPhotons();
  this->parsePhotoGeneration();
  this->parseSourceType();
  this->parseDeposition();
  this->parseIntersectionEB();
  this->parsePlotVariables();
  this->parseInstantaneous();
  this->parseDivergenceComputation();
}

void
McPhoto::parseRuntimeOptions()
{
  CH_TIME("McPhoto::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }

  this->parseVerbosity();
  this->parseTransparentBoundaries();
  this->parsePseudoPhotons();
  this->parsePhotoGeneration();
  this->parseSourceType();
  this->parseDeposition();
  this->parseIntersectionEB();
  this->parsePlotVariables();
  this->parseInstantaneous();
  this->parseDivergenceComputation();
}

void
McPhoto::parseTransparentBoundaries()
{
  CH_TIME("McPhoto::parseTransparentBoundaries");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseTransparentBoundaries" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("transparent_eb", m_transparentEB);
}

void
McPhoto::parseDivergenceComputation()
{
  CH_TIME("McPhoto::parseDivergenceComputation");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseDivergenceComputation" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("blend_conservation", m_blendConservation);
}

void
McPhoto::parseInstantaneous()
{
  CH_TIME("McPhoto::parseInstantaneous");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseInstantaneous" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("instantaneous", m_instantaneous);
}

void
McPhoto::parsePseudoPhotons()
{
  CH_TIME("McPhoto::parsePseudoPhotons");
  if (m_verbosity > 5) {
    pout() << m_name + "::parsePseudoPhotons" << endl;
  }

  ParmParse pp(m_className.c_str());

  int maxPhotons;

  pp.get("num_sampling_packets", m_numSamplingPackets);
  pp.get("max_photons_per_cell", maxPhotons);

  if (maxPhotons <= 0) {
    m_maxPhotonsGeneratedPerCell = std::numeric_limits<size_t>::max();
  }
  else {
    m_maxPhotonsGeneratedPerCell = (size_t)maxPhotons;
  }

  if (m_numSamplingPackets <= 0) {
    m_numSamplingPackets = 1;
  }
}

void
McPhoto::parsePhotoGeneration()
{
  CH_TIME("McPhoto::parsePhotoGeneration");
  if (m_verbosity > 5) {
    pout() << m_name + "::parsePhotoGeneration" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("photon_generation", str);

  if (str == "deterministic") {
    m_photoGenerationMethod = PhotonGeneration::Deterministic;
  }
  else if (str == "stochastic") {
    m_photoGenerationMethod = PhotonGeneration::Stochastic;
  }
  else {
    MayDay::Error("McPhoto::set_PhotonGeneration - unknown photon generation type requested");
  }
}

void
McPhoto::parseSourceType()
{
  CH_TIME("McPhoto::parseSourceType");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseSourceType" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("source_type", str);

  if (str == "number") {
    m_sourceType = SourceType::Number;
  }
  else if (str == "volume") {
    m_sourceType = SourceType::PerVol;
  }
  else if (str == "volume_rate") {
    m_sourceType = SourceType::PerVolSecond;
  }
  else if (str == "rate") {
    m_sourceType = SourceType::PerSecond;
  }
  else {
    MayDay::Error("McPhoto::setSourceType - unknown source type requested");
  }
}

void
McPhoto::parseIntersectionEB()
{
  CH_TIME("McPhoto::parseIntersectionEB");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseIntersectionEB" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;

  pp.get("intersection_alg", str);
  pp.get("bisect_step", m_bisectStep);

  if (str == "raycast") {
    m_intersectionEB = IntersectionEB::Raycast;
  }
  else if (str == "bisection") {
    m_intersectionEB = IntersectionEB::Bisection;
  }
  else {
    MayDay::Error("McPhoto::parseIntersectionEB -- logic bust");
  }
}

void
McPhoto::parseDeposition()
{
  CH_TIME("McPhoto::parseDeposition");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseDeposition" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("deposition", str);

  m_depositNumber = false;
  if (str == "num") {
    m_deposition    = DepositionType::NGP;
    m_depositNumber = true;
  }
  else if (str == "ngp") {
    m_deposition = DepositionType::NGP;
  }
  else if (str == "cic") {
    m_deposition = DepositionType::CIC;
  }
  else {
    MayDay::Error("McPhoto::parseDeposition - unsupported deposition method requested");
  }

  pp.get("deposition_cf", str);
  if (str == "interp") {
    m_coarseFineDeposition = CoarseFineDeposition::Interp;
  }
  else if (str == "halo") {
    m_coarseFineDeposition = CoarseFineDeposition::Halo;
  }
  else if (str == "halo_ngp") {
    m_coarseFineDeposition = CoarseFineDeposition::HaloNGP;
  }
  else {
    MayDay::Error("McPhoto::parseDeposition - unknown coarse-fine deposition method requested.");
  }
}

void
McPhoto::parsePlotVariables()
{
  CH_TIME("McPhoto::parsePlotVariables");
  if (m_verbosity > 5) {
    pout() << m_name + "::parsePlotVariables" << endl;
  }

  m_plotPhi           = false;
  m_plotSource        = false;
  m_plotPhotons       = false;
  m_plotBulkPhotons   = false;
  m_plotEBPhotons     = false;
  m_plotDomainPhotons = false;
  m_plotSourcePhotons = false;

  ParmParse           pp(m_className.c_str());
  const int           num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++) {
    if (str[i] == "phi") {
      m_plotPhi = true;
    }
    else if (str[i] == "src") {
      m_plotSource = true;
    }
    else if (str[i] == "phot") {
      m_plotPhotons = true;
    }
    else if (str[i] == "bulk_phot") {
      m_plotBulkPhotons = true;
    }
    else if (str[i] == "eb_phot") {
      m_plotEBPhotons = true;
    }
    else if (str[i] == "dom_phot") {
      m_plotDomainPhotons = true;
    }
    else if (str[i] == "src_phot") {
      m_plotSourcePhotons = true;
    }
  }
}

void
McPhoto::clear(const WhichContainer& a_which)
{
  CH_TIME("McPhoto::clear(WhichContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::clear(WhichContainer)" << endl;
  }

  switch (a_which) {
  case WhichContainer::Photons: {
    m_photons.clearParticles();

    break;
  }
  case WhichContainer::Bulk: {
    m_bulkPhotons.clearParticles();

    break;
  }
  case WhichContainer::EB: {
    m_ebPhotons.clearParticles();

    break;
  }
  case WhichContainer::Domain: {
    m_domainPhotons.clearParticles();

    break;
  }
  case WhichContainer::Source: {
    m_sourcePhotons.clearParticles();

    break;
  }
  default: {
    MayDay::Error("McPhoto::sortPhotonsByCell -- logic bust");

    break;
  }
  }
}

void
McPhoto::clear()
{
  CH_TIME("McPhoto::clear()");
  if (m_verbosity > 5) {
    pout() << m_name + "::clear()" << endl;
  }

  this->clear(m_photons);
}

void
McPhoto::clear(ParticleContainer<Photon>& a_photons)
{
  CH_TIME("McPhoto::clear(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::clear(ParticleContainer)" << endl;
  }

  this->clear(a_photons.getParticles());
}

void
McPhoto::clear(AMRParticles<Photon>& a_photons)
{
  CH_TIME("McPhoto::clear(AMRParticles)");
  if (m_verbosity > 5) {
    pout() << m_name + "::clear(AMRParticles)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    a_photons[lvl]->clear();
  }
}

void
McPhoto::allocate()
{
  CH_TIME("McPhoto::allocate");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocate" << endl;
  }

  // Allocate mesh data
  m_amr->allocate(m_phi, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_source, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_scratch, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_depositionNC, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_massDiff, m_realm, m_phase, m_nComp);

  // Allocate particle data holders
  m_amr->allocate(m_photons, m_realm);
  m_amr->allocate(m_bulkPhotons, m_realm);
  m_amr->allocate(m_ebPhotons, m_realm);
  m_amr->allocate(m_domainPhotons, m_realm);
  m_amr->allocate(m_sourcePhotons, m_realm);
}

void
McPhoto::preRegrid(const int a_lmin, const int a_oldFinestLevel)
{
  CH_TIME("McPhoto::preRegrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::pre_grid" << endl;
  }

  m_photons.preRegrid(a_lmin);
  m_bulkPhotons.preRegrid(a_lmin);
  m_ebPhotons.preRegrid(a_lmin);
  m_domainPhotons.preRegrid(a_lmin);
  m_sourcePhotons.preRegrid(a_lmin);

  this->deallocate();
}

void
McPhoto::deallocate()
{
  CH_TIME("McPhoto::deallocate");
  if (m_verbosity > 5) {
    pout() << m_name + "::deallocate" << endl;
  }

  m_phi.clear();
  m_source.clear();
  m_scratch.clear();
  m_depositionNC.clear();
  m_massDiff.clear();
}

void
McPhoto::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("McPhoto::regrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::regrid" << endl;
  }

  // Mesh data regrids
  m_amr->allocate(m_phi, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_source, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_scratch, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_depositionNC, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_massDiff, m_realm, m_phase, m_nComp);

  m_amr->remapToNewGrids(m_photons, a_lmin, a_newFinestLevel);
  m_amr->remapToNewGrids(m_bulkPhotons, a_lmin, a_newFinestLevel);
  m_amr->remapToNewGrids(m_ebPhotons, a_lmin, a_newFinestLevel);
  m_amr->remapToNewGrids(m_domainPhotons, a_lmin, a_newFinestLevel);
  m_amr->remapToNewGrids(m_sourcePhotons, a_lmin, a_newFinestLevel);

  // Deposit
  this->depositPhotons();
}

void
McPhoto::sortPhotonsByCell(const WhichContainer& a_which)
{
  CH_TIME("McPhoto::sortPhotonsByCell(WhichContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::sortPhotonsByCell(WhichContainer)" << endl;
  }

  switch (a_which) {
  case WhichContainer::Photons: {
    m_photons.organizeParticlesByCell();

    break;
  }
  case WhichContainer::Bulk: {
    m_bulkPhotons.organizeParticlesByCell();

    break;
  }
  case WhichContainer::EB: {
    m_ebPhotons.organizeParticlesByCell();

    break;
  }
  case WhichContainer::Domain: {
    m_domainPhotons.organizeParticlesByCell();

    break;
  }
  case WhichContainer::Source: {
    m_sourcePhotons.organizeParticlesByCell();

    break;
  }
  default: {
    MayDay::Error("McPhoto::sortPhotonsByCell -- logic bust");

    break;
  }
  }
}

void
McPhoto::sortPhotonsByPatch(const WhichContainer& a_which)
{
  CH_TIME("McPhoto::sortPhotonsByPatch(WhichContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::sortPhotonsByPatch(WhichContainer)" << endl;
  }

  switch (a_which) {
  case WhichContainer::Photons: {
    m_photons.organizeParticlesByPatch();

    break;
  }
  case WhichContainer::Bulk: {
    m_bulkPhotons.organizeParticlesByPatch();

    break;
  }
  case WhichContainer::EB: {
    m_ebPhotons.organizeParticlesByPatch();

    break;
  }
  case WhichContainer::Domain: {
    m_domainPhotons.organizeParticlesByPatch();

    break;
  }
  case WhichContainer::Source: {
    m_sourcePhotons.organizeParticlesByPatch();

    break;
  }
  default: {
    MayDay::Error("McPhoto::sortPhotonsByPatch -- logic bust");

    break;
  }
  }
}

void
McPhoto::registerOperators()
{
  CH_TIME("McPhoto::registerOperators");
  if (m_verbosity > 5) {
    pout() << m_name + "::registerOperators" << endl;
  }

  if (m_amr.isNull()) {
    MayDay::Error("McPhoto::registerOperators - m_amr is null, need to set AmrMesh.");
  }
  else {
    m_amr->registerOperator(s_eb_coar_ave, m_realm, m_phase);
    m_amr->registerOperator(s_eb_fill_patch, m_realm, m_phase);
    m_amr->registerOperator(s_eb_redist, m_realm, m_phase);
    m_amr->registerOperator(s_particle_mesh, m_realm, m_phase);
    m_amr->registerOperator(s_noncons_div, m_realm, m_phase);

    // For CIC deposition
    m_amr->registerMask(s_particle_halo, m_haloBuffer, m_realm);
  }
}

void
McPhoto::computeBoundaryFlux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi)
{
  CH_TIME("McPhoto::computeBoundaryFlux");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeBoundaryFlux" << endl;
  }

  DataOps::setValue(a_ebFlux, 0.0);
}

void
McPhoto::computeDomainFlux(EBAMRIFData& a_domainFlux, const EBAMRCellData& a_phi)
{
  CH_TIME("McPhoto::computeDomainFlux");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDomainFlux" << endl;
  }

  DataOps::setValue(a_domainFlux, 0.0);
}

void
McPhoto::computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_phi)
{
  CH_TIME("McPhoto::computeFlux");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeFlux" << endl;
  }

  const std::string str =
    "McPhoto::computeFlux - Fluid flux can't be computed with discrete photons. Calling this is an error";

  MayDay::Error(str.c_str());
}

void
McPhoto::computeDensity(EBAMRCellData& a_isotropic, const EBAMRCellData& a_phi)
{
  CH_TIME("McPhoto::computeDensity");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDensity" << endl;
  }

  MayDay::Error(
    "McPhoto::computeDensity - can't compute density from RTE solution with discrete photons. Maybe you meant to deposit photons?");
}

void
McPhoto::writePlotFile()
{
  CH_TIME("McPhoto::writePlotFile");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotFile" << endl;
  }

  MayDay::Error("McPhoto::writePlotFile - not implemented for McPhoto (yet)");
}

#ifdef CH_USE_HDF5
void
McPhoto::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const
{
  CH_TIME("McPhoto::writeCheckpointLevel");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeCheckpointLevel" << endl;
  }

  // Write state vector
  write(a_handle, *m_phi[a_level], m_name);

  // Write particles. Must be implemented.
  std::string str = m_name + "_particles";
  writeParticlesToHDF(a_handle, m_photons[a_level], str);
}
#endif

#ifdef CH_USE_HDF5
void
McPhoto::readCheckpointLevel(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("McPhoto::readCheckpointLevel");
  if (m_verbosity > 5) {
    pout() << m_name + "::readCheckpointLevel" << endl;
  }

  // Read state vector
  read<EBCellFAB>(a_handle, *m_phi[a_level], m_name, m_amr->getGrids(m_realm)[a_level], Interval(0, 0), false);

  // Read particles. Should be implemented
  std::string str = m_name + "_particles";
  readParticlesFromHDF(a_handle, m_photons[a_level], str);
}
#endif

Vector<std::string>
McPhoto::getPlotVariableNames() const
{
  CH_TIME("McPhoto::getPlotVariableNames");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPlotVariableNames" << endl;
  }

  Vector<std::string> plotVarNames(0);

  if (m_plotPhi) {
    plotVarNames.push_back(m_name + " phi");
  }
  if (m_plotSource) {
    plotVarNames.push_back(m_name + " source");
  }
  if (m_plotPhotons) {
    plotVarNames.push_back(m_name + " active_photons");
  }
  if (m_plotBulkPhotons) {
    plotVarNames.push_back(m_name + " bulk_photons");
  }
  if (m_plotEBPhotons) {
    plotVarNames.push_back(m_name + " eb_photons");
  }
  if (m_plotDomainPhotons) {
    plotVarNames.push_back(m_name + " domain_photons");
  }
  if (m_plotSourcePhotons) {
    plotVarNames.push_back(m_name + " source_photons");
  }

  return plotVarNames;
}

int
McPhoto::getNumberOfPlotVariables() const
{
  CH_TIME("McPhoto::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }

  int numPlotVars = 0;

  if (m_plotPhi) {
    numPlotVars += 1;
  }
  if (m_plotSource) {
    numPlotVars += 1;
  }
  if (m_plotPhotons) {
    numPlotVars += 1;
  }
  if (m_plotBulkPhotons) {
    numPlotVars += 1;
  }
  if (m_plotEBPhotons) {
    numPlotVars += 1;
  }
  if (m_plotDomainPhotons) {
    numPlotVars += 1;
  }
  if (m_plotSourcePhotons) {
    numPlotVars += 1;
  }

  return numPlotVars;
}

int
McPhoto::domainBcMap(const int a_dir, const Side::LoHiSide a_side)
{
  const int iside = (a_side == Side::Lo) ? 0 : 1;

  return 2 * a_dir + iside;
}

Real
McPhoto::randomExponential(const Real a_mean)
{
  std::exponential_distribution<Real> dist(a_mean);
  return Random::get(dist);
}

void
McPhoto::computeNumPhysicalPhotons(EBAMRCellData&       a_numPhysPhotonsTotal,
                                   EBAMRCellData&       a_numPhysPhotonsPacket,
                                   const EBAMRCellData& a_source,
                                   const Real           a_dt) const noexcept
{
  CH_TIME("McPhoto::computeNumPhysicalPhotons");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeNumPhysicalPhotons" << endl;
  }

  CH_assert(a_numPhysPhotonsTotal[0]->nComp() == 1);
  CH_assert(a_numPhysPhotonsPacket[0]->nComp() == m_numSamplingPackets);
  CH_assert(a_source[0]->nComp() == 1);
  CH_assert(a_dt >= 0.0);

  DataOps::setValue(a_numPhysPhotonsTotal, 0.0);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];
    const Real               vol   = pow(dx, SpaceDim);

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const Box            cellBox    = dbl[dit()];
      const EBISBox&       ebisbox    = ebisl[dit()];
      const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_realm)[lvl])[dit()];

      const EBCellFAB& source    = (*a_source[lvl])[dit()];
      const FArrayBox& sourceReg = source.getFArrayBox();

      EBCellFAB& numPhysPhotonsTotal    = (*a_numPhysPhotonsTotal[lvl])[dit()];
      FArrayBox& numPhysPhotonsTotalReg = numPhysPhotonsTotal.getFArrayBox();

      EBCellFAB& numPhysPhotonsPacket    = (*a_numPhysPhotonsPacket[lvl])[dit()];
      FArrayBox& numPhysPhotonsPacketReg = numPhysPhotonsPacket.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        if (ebisbox.isRegular(iv) && validCells(iv)) {

          const size_t numPhysPhotons = this->drawPhotons(sourceReg(iv, 0), vol, a_dt);
          const size_t packetSize     = numPhysPhotons / m_numSamplingPackets;
          const size_t remainder      = numPhysPhotons % m_numSamplingPackets;

          numPhysPhotonsTotalReg(iv, 0) = Real(numPhysPhotons);

          for (int i = 0; i < m_numSamplingPackets; i++) {
            numPhysPhotonsPacketReg(iv, i) = Real(packetSize);
          }
          numPhysPhotonsPacketReg(iv, 0) += Real(remainder);
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();

        if (ebisbox.isIrregular(iv) && validCells(iv)) {
          const size_t numPhysPhotons = this->drawPhotons(source(vof, 0), vol, a_dt);
          const size_t packetSize     = numPhysPhotons / m_numSamplingPackets;
          const size_t remainder      = numPhysPhotons % m_numSamplingPackets;

          numPhysPhotonsTotal(vof, 0) = Real(numPhysPhotons);

          for (int i = 0; i < m_numSamplingPackets; i++) {
            numPhysPhotonsPacket(vof, i) = Real(packetSize);
          }
          numPhysPhotonsPacket(vof, 0) += Real(remainder);
        }
      };

      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

size_t
McPhoto::drawPhotons(const Real a_source, const Real a_volume, const Real a_dt) const noexcept
{
  CH_TIME("McPhoto::drawPhotons");
  if (m_verbosity > 5) {
    pout() << m_name + "::drawPhotons" << endl;
  }

  size_t numPhysicalPhotons = 0;

  // Check if we need any type of source term normalization
  Real factor;
  if (m_sourceType == SourceType::Number) {
    factor = 1.0;
  }
  else if (m_sourceType == SourceType::PerVol) {
    factor = a_volume;
  }
  else if (m_sourceType == SourceType::PerVolSecond) {
    factor = a_volume * a_dt;
  }
  else if (m_sourceType == SourceType::PerSecond) {
    factor = a_dt;
  }
  else {
    factor = 0.0;

    MayDay::Error("McPhoto::drawPhotons -- logic bust");
  }

  // Draw a number of Photons with the desired algorithm
  if (m_photoGenerationMethod == PhotonGeneration::Stochastic) {
    const Real mean    = a_source * factor;
    numPhysicalPhotons = Random::getPoisson<size_t>(mean);
  }
  else if (m_photoGenerationMethod == PhotonGeneration::Deterministic) {
    numPhysicalPhotons = round(a_source * factor);
  }
  else {
    MayDay::Error("mc::drawPhotons - unknown generation requested. Aborting...");
  }

  return numPhysicalPhotons;
}

void
McPhoto::generateComputationalPhotons(ParticleContainer<Photon>& a_photons,
                                      const EBAMRCellData&       a_numPhysPhotons,
                                      const size_t               a_maxPhotonsPerCell) const noexcept
{
  CH_TIME("McPhoto::generateComputationalPhotons");
  if (m_verbosity > 5) {
    pout() << m_name + "::generateComputationalPhotons" << endl;
  }

  CH_assert(a_numPhysPhotons[0]->nComp() == 1);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const RealVect           probLo = m_amr->getProbLo();
    const Real               dx     = m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box      cellBox = dbl[dit()];
      const EBISBox& ebisbox = ebisl[dit()];

      const EBCellFAB&     numPhysPhotons    = (*a_numPhysPhotons[lvl])[dit()];
      const FArrayBox&     numPhysPhotonsReg = numPhysPhotons.getFArrayBox();
      const BaseFab<bool>& validCells        = (*m_amr->getValidCells(m_realm)[lvl])[dit()];

      List<Photon>& photons = a_photons[lvl][dit()].listItems();

      // Regular cells. Note that we make superphotons if we have to. Also, only draw photons in valid cells,
      // grids that are covered by finer grids don't draw photons.
      auto regularKernel = [&](const IntVect& iv) -> void {
        if (ebisbox.isRegular(iv) && validCells(iv, 0)) {

          const size_t num = numPhysPhotonsReg(iv, 0);

          if (num > 0) {
            const std::vector<size_t> photonWeights = ParticleManagement::partitionParticleWeights(num,
                                                                                                   a_maxPhotonsPerCell);

            const RealVect lo = probLo + RealVect(iv) * dx;
            const RealVect hi = lo + RealVect::Unit * dx;

            for (size_t i = 0; i < photonWeights.size(); i++) {

              // Determine starting position within cell, propagation direction, absorption
              // length, and weight.
              const RealVect pos    = Random::randomPosition(lo, hi);
              const RealVect v      = Units::c * Random::getDirection();
              const Real     weight = (Real)photonWeights[i];
              const Real     kappa  = m_rtSpecies->getAbsorptionCoefficient(pos);

              photons.add(Photon(pos, v / v.vectorLength(), kappa, weight));
            }
          }
        }
      };

      // Irregular kernel. Same as the above really.
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();

        if (validCells(iv, 0)) {

          const size_t num = numPhysPhotons(vof, 0);

          if (num > 0) {
            const std::vector<size_t> photonWeights = ParticleManagement::partitionParticleWeights(num,
                                                                                                   a_maxPhotonsPerCell);

            // These are needed when drawing photon starting positions within cut-cells -- we compute the
            // minimum bounding box when we draw the position within the valid region of the cut-cell.
            const Real     volFrac       = ebisbox.volFrac(vof);
            const RealVect cellPos       = probLo + Location::position(Location::Cell::Center, vof, ebisbox, dx);
            const RealVect bndryCentroid = ebisbox.bndryCentroid(vof);
            const RealVect bndryNormal   = ebisbox.normal(vof);

            RealVect lo = -0.5 * RealVect::Unit;
            RealVect hi = 0.5 * RealVect::Unit;
            if (volFrac < 1.0) {
              DataOps::computeMinValidBox(lo, hi, bndryNormal, bndryCentroid);
            }

            for (size_t i = 0; i < photonWeights.size(); i++) {

              // Determine starting position within cell, propagation direction, absorption
              // length, and weight.
              const RealVect pos    = Random::randomPosition(cellPos, lo, hi, bndryCentroid, bndryNormal, dx, volFrac);
              const RealVect v      = Units::c * Random::getDirection();
              const Real     weight = (Real)photonWeights[i];

              photons.add(Photon(pos, v, m_rtSpecies->getAbsorptionCoefficient(pos), weight));
            }
          }
        }
      };

      // Run the kernels.
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

void
McPhoto::depositPhotons()
{
  CH_TIME("McPhoto::depositPhotons()");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositPhotons()" << endl;
  }

  this->depositPhotons(m_phi, m_photons, m_deposition);
}

void
McPhoto::depositPhotons(EBAMRCellData& a_phi, ParticleContainer<Photon>& a_photons, const DepositionType& a_deposition)
{
  CH_TIME("McPhoto::depositPhotons(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositPhotons(ParticleContainer)" << endl;
  }

  // a_phi contains only weights, i.e. not divided by kappa
  this->depositKappaConservative(a_phi, a_photons, a_deposition, m_coarseFineDeposition);

  // Compute m_depositionNC = sum(kappa*Wc)/sum(kappa)
  this->depositNonConservative(m_depositionNC, a_phi);

  // Compute hybrid deposition, including mass differnce
  this->depositHybrid(a_phi, m_massDiff, m_depositionNC);

  // Redistribute
  if (m_blendConservation) {
    Vector<RefCountedPtr<EBFluxRedistribution>>& redistOps = m_amr->getRedistributionOp(m_realm, m_phase);
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real     scale     = 1.0;
      const Interval variables = Interval(0, 0);
      const bool     hasCoar   = lvl > 0;
      const bool     hasFine   = lvl < m_amr->getFinestLevel();

      if (hasCoar) {
        redistOps[lvl]->redistributeCoar(*a_phi[lvl - 1], *m_massDiff[lvl], scale, variables);
      }

      redistOps[lvl]->redistributeLevel(*a_phi[lvl], *m_massDiff[lvl], scale, variables);

      if (hasFine) {
        redistOps[lvl]->redistributeFine(*a_phi[lvl + 1], *m_massDiff[lvl], scale, variables);
      }
    }
  }

  // Average down and interpolate
  m_amr->conservativeAverage(a_phi, m_realm, m_phase);
  m_amr->interpGhost(a_phi, m_realm, m_phase);
}

void
McPhoto::depositKappaConservative(EBAMRCellData&             a_phi,
                                  ParticleContainer<Photon>& a_particles,
                                  const DepositionType       a_deposition,
                                  const CoarseFineDeposition a_coarseFineDeposition)
{
  CH_TIME("McPhoto::depositKappaConservative");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositKappaConservative" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);

  switch (a_coarseFineDeposition) {
  case CoarseFineDeposition::Interp: {
    m_amr->depositParticles<Photon, &Photon::weight>(a_phi,
                                                     m_realm,
                                                     m_phase,
                                                     a_particles,
                                                     a_deposition,
                                                     CoarseFineDeposition::Interp,
                                                     false);

    break;
  }
  case CoarseFineDeposition::Halo: {

    // Copy particles living on the mask.
    const AMRMask& mask = m_amr->getMask(s_particle_halo, m_haloBuffer, m_realm);
    a_particles.copyMaskParticles(mask);

    m_amr->depositParticles<Photon, &Photon::weight>(a_phi,
                                                     m_realm,
                                                     m_phase,
                                                     a_particles,
                                                     a_deposition,
                                                     CoarseFineDeposition::Halo,
                                                     false);

    // Clear out the mask particles.
    a_particles.clearMaskParticles();

    break;
  }
  case CoarseFineDeposition::HaloNGP: {
    const AMRMask& mask = m_amr->getMask(s_particle_halo, m_haloBuffer, m_realm);

    // Transfer particles living on the mask.
    a_particles.transferMaskParticles(mask);

    m_amr->depositParticles<Photon, &Photon::weight>(a_phi,
                                                     m_realm,
                                                     m_phase,
                                                     a_particles,
                                                     a_deposition,
                                                     CoarseFineDeposition::HaloNGP,
                                                     false);

    // Transfer them back.
    a_particles.transferParticles(a_particles.getMaskParticles());

    break;
  }
  default: {
    MayDay::Error("McPhoto::depositKappaConservative -- logic bust due to unsupported coarse-fine deposition");
  }
  }
}

void
McPhoto::depositNonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC)
{
  CH_TIME("McPhoto::depositNonConservative");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositNonConservative" << endl;
  }

  if (m_blendConservation) {
    const IrregAmrStencil<NonConservativeDivergenceStencil>& stencils =
      m_amr->getNonConservativeDivergenceStencils(m_realm, m_phase);
    stencils.apply(a_depositionNC, a_depositionKappaC);
  }
  else {
    DataOps::setValue(a_depositionNC, 0.0);
  }
}

void
McPhoto::depositHybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_massDifference, const EBAMRIVData& a_depositionNC)
{
  CH_TIME("McPhoto::depositHybrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositHybrid" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      EBCellFAB&             divH   = (*a_depositionH[lvl])[dit()]; // On input, this contains kappa*depositionWeights
      BaseIVFAB<Real>&       deltaM = (*a_massDifference[lvl])[dit()];
      const BaseIVFAB<Real>& divNC  = (*a_depositionNC[lvl])[dit()];

      const Box      box     = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];

      // Iteration space for kernel.
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      // Kernel
      auto kernel = [&](const VolIndex& vof) -> void {
        const Real kappa = ebisbox.volFrac(vof);
        const Real dc    = divH(vof, m_comp);
        const Real dnc   = divNC(vof, m_comp);

        // Note that if dc - kappa*dnc can be negative, i.e. we may end up STEALING mass
        // from other cells. This is why there is a flag m_blendConservation which always
        // gives positive definite results.
        divH(vof, m_comp)   = dc + (1 - kappa) * dnc;           // On output, contains hybrid divergence
        deltaM(vof, m_comp) = (1 - kappa) * (dc - kappa * dnc); // Remember, dc already scaled by kappa.
      };

      BoxLoops::loop(vofit, kernel);
    }
  }
}

void
McPhoto::depositPhotonsNGP(LevelData<EBCellFAB>&            a_output,
                           const ParticleContainer<Photon>& a_photons,
                           const int                        a_level) const noexcept
{
  CH_TIME("McPhoto::depositPhotonsNGP");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositPhotonsNGP" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  const ProblemDomain&     domain = m_amr->getDomains()[a_level];
  const DisjointBoxLayout& dbl    = m_amr->getGrids(a_photons.getRealm())[a_level];
  const EBISLayout&        ebisl  = m_amr->getEBISLayout(a_photons.getRealm(), m_phase)[a_level];
  const Real               dx     = m_amr->getDx()[a_level];
  const RealVect           probLo = m_amr->getProbLo();

  CH_assert(a_output.disjointBoxLayout() == dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box      cellBox = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];

    EBParticleMesh particleMesh(domain, cellBox, ebisbox, dx * RealVect::Unit, probLo);

    EBCellFAB&          output  = a_output[dit()];
    const List<Photon>& photons = a_photons[a_level][dit()].listItems();

    particleMesh.deposit<Photon, &Photon::weight>(photons, output, DepositionType::NGP, true);
  }
}

void
McPhoto::advancePhotonsInstantaneous(ParticleContainer<Photon>& a_bulkPhotons,
                                     ParticleContainer<Photon>& a_ebPhotons,
                                     ParticleContainer<Photon>& a_domainPhotons,
                                     ParticleContainer<Photon>& a_photons)
{
  CH_TIME("McPhoto::advancePhotonsInstantaneous");
  if (m_verbosity > 5) {
    pout() << m_name + "::advancePhotonsInstantaneous" << endl;
  }

  // TLDR: This routine iterates over the levels and boxes and does the following
  //
  //       Forall photons in a_photons: {
  //          1. Draw random absorption position
  //          2. Determine if path intersected boundary, either EB or domain
  //          3. Move the photon to appropriate data holder:
  //                 Path crossed EB   => a_ebPhotons
  //                 Path cross domain => a_domainPhotons
  //                 Absorbed in bulk  => a_bulkPhotons
  //       }
  //
  //       Remap a_bulkPhotons, a_ebPhotons, a_domainPhotons

  // Low and high corners
  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  // Safety factor -- weird thing but if particles lie EXACTLY on the high side of the domain boundary they map to the outside
  // of the computational domain (for exactly the same reason that a photon on the interface between two cells map to the "upper one"). Adding
  // a negligible safety factor to prevent that from happening.
  constexpr Real SAFETY = 1.E-6;

  // This is the implicit function used for intersection tests
  const RefCountedPtr<BaseIF>& impFunc = m_computationalGeometry->getImplicitFunction(m_phase);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const Real               dx  = m_amr->getDx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      List<Photon>& bulkPhotons = a_bulkPhotons[lvl][dit()].listItems();
      List<Photon>& ebPhotons   = a_ebPhotons[lvl][dit()].listItems();
      List<Photon>& domPhotons  = a_domainPhotons[lvl][dit()].listItems();
      List<Photon>& allPhotons  = a_photons[lvl][dit()].listItems();

      // Iterate over the Photons that will be moved.
      for (ListIterator<Photon> lit(allPhotons); lit.ok(); ++lit) {
        Photon& p = lit();

        // Draw a new random absorption position
        const RealVect& oldPos    = p.position();
        const RealVect& direction = p.velocity() / (p.velocity().vectorLength());
        const RealVect  newPos    = oldPos + direction * this->randomExponential(p.kappa());

        // Check if we should check of different types of boundary intersections. These are cheap initial tests that allow
        // us to skip intersection tests for some photons.
        bool checkEB  = false;
        bool checkDom = false;

        if (!impFunc.isNull()) {
          checkEB = true;
        }
        for (int dir = 0; dir < SpaceDim; dir++) {
          if (newPos[dir] < probLo[dir] || newPos[dir] > probHi[dir]) {
            checkDom = true;
          }
        }

        if ((!checkEB && !checkDom) || m_transparentEB) {
          p.position() = newPos;
          bulkPhotons.add(p);
        }
        else {
          // Must do an intersection test (with either EB or domain). These tests work such that we parametrize the photon path as
          //
          // x(s) = x0 + s*(x1-x0), s = [0,1]
          //
          // where x0 is the starting position (oldPos) and x1 is the new position (newPos).
          //
          // We determine s for the domain boundaries and EBs. If the photon path cross both, smallest s takes precedence.

          Real sDom = std::numeric_limits<Real>::max();
          Real sEB  = std::numeric_limits<Real>::max();

          bool contactDomain = false;
          bool contactEB     = false;

          // Do intersection tests. These return true/false if the path crossed an object. If it returned true, the s-parameter
          // will have been defined as well.
          if (checkDom) {
            contactDomain = ParticleOps::domainIntersection(oldPos, newPos, probLo, probHi, sDom);
          }

          if (checkEB) {
            switch (m_intersectionEB) {
            case IntersectionEB::Raycast: {
              contactEB = ParticleOps::ebIntersectionRaycast(impFunc, oldPos, newPos, 1.E-3 * dx, sEB);

              break;
            }
            case IntersectionEB::Bisection: {
              contactEB = ParticleOps::ebIntersectionBisect(impFunc, oldPos, newPos, m_bisectStep, sEB);

              break;
            }
            default: {
              MayDay::Error("McPhoto::advancePhotonsInstantenous -- logic bust in eb intersection");

              break;
            }
            }
          }

          // Move the photon to the appropriate data holder
          if (!contactEB && !contactDomain) {
            p.position() = newPos;
            bulkPhotons.add(p);
          }
          else {
            const RealVect path = newPos - oldPos;

            if (sEB < sDom) {
              p.position() = oldPos + sEB * path;

              ebPhotons.add(p);
            }
            else {
              p.position() = oldPos + std::max((Real)0.0, sDom - SAFETY) * path;

              domPhotons.add(p);
            }
          }
        }
      }

      // Everything will have been absorbed, but we used ::add rather than ::transfer, so just clear out the starting photons.
      allPhotons.clear();
    }
  }

  // Need to remap because photons may/will have moved off the processor.
  a_bulkPhotons.remap();
  a_ebPhotons.remap();
  a_domainPhotons.remap();
}

void
McPhoto::advancePhotonsTransient(ParticleContainer<Photon>& a_bulkPhotons,
                                 ParticleContainer<Photon>& a_ebPhotons,
                                 ParticleContainer<Photon>& a_domainPhotons,
                                 ParticleContainer<Photon>& a_photons,
                                 const Real                 a_dt)
{
  CH_TIME("McPhoto::advancePhotonsTransient");
  if (m_verbosity > 5) {
    pout() << m_name + "::advancePhotonsTransient" << endl;
  }

  // TLDR: This routine iterates over the levels and boxes and does the following
  //
  //       Forall photons in a_photons: {
  //          1. Check new photon position
  //          2. Check if the photon was absorbed on the interval
  //          3. Check if path intersects boundary, either EB or domain
  //          4. Move the photon to appropriate data holder:
  //                 Path crossed EB   => a_ebPhotons
  //                 Path cross domain => a_domainPhotons
  //                 Absorbed in bulk  => a_bulkPhotons
  //       }
  //
  //       Remap a_bulkPhotons, a_ebPhotons, a_domainPhotons, a_photons

  // Low and high corners
  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  // Safety factor -- weird thing but if particles lie EXACTLY on the high side of the domain boundary they map to the outside
  // of the computational domain (for exactly the same reason that a photon on the interface between two cells map to the "upper one"). Adding
  // a negligible safety factor to prevent that from happening.
  constexpr Real SAFETY = 1.E-8;

  // This is the implicit function used for intersection tests
  const RefCountedPtr<BaseIF>& impFunc = m_computationalGeometry->getImplicitFunction(m_phase);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const Real               dx  = m_amr->getDx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      List<Photon>& bulkPhotons = a_bulkPhotons[lvl][dit()].listItems();
      List<Photon>& ebPhotons   = a_ebPhotons[lvl][dit()].listItems();
      List<Photon>& domPhotons  = a_domainPhotons[lvl][dit()].listItems();
      List<Photon>& allPhotons  = a_photons[lvl][dit()].listItems();

      // Iterate over the photons that will be moved.
      for (ListIterator<Photon> lit(allPhotons); lit.ok(); ++lit) {
        Photon& p = lit();

        // Move the Photon
        const RealVect oldPos  = p.position();
        const RealVect v       = p.velocity();
        const RealVect newPos  = oldPos + v * a_dt;
        const RealVect path    = newPos - oldPos;
        const Real     pathLen = path.vectorLength();

        // Check if we should check of different types of boundary intersections. These are checp initial tests that allow
        // us to skip more expensive intersection tests for some Photons.
        bool checkEB  = false;
        bool checkDom = false;

        if (!impFunc.isNull()) {
          checkEB = true;
        }
        for (int dir = 0; dir < SpaceDim; dir++) {
          if (newPos[dir] <= probLo[dir] || newPos[dir] >= probHi[dir]) {
            checkDom = true;
          }
        }

        // Return flags for particle intersection tests.
        bool absorbedBulk   = false;
        bool absorbedEB     = false;
        bool absorbedDomain = false;

        // Below, we do an intersection test (with either EB or domain). These tests work such that we parametrize the photon path as
        //
        // x(s) = x0 + s*(x1-x0), s = [0,1]
        //
        // where x0 is the starting position (oldPos) and x1 is the new position (newPos).
        //
        // We determine s for the domain boundaries and EBs. If the photon path cross both, smallest s takes precedence.

        // The s-parameters for various absorption scenarios.
        Real sBulk   = std::numeric_limits<Real>::max();
        Real sEB     = std::numeric_limits<Real>::max();
        Real sDomain = std::numeric_limits<Real>::max();

        // Check absorption in the bulk. We draw a propagation distance, if the photon propagates longer
        // than this distance the photon is absorbed.
        const Real travelLen = this->randomExponential(p.kappa());
        if (travelLen < pathLen) {
          absorbedBulk = true;
          sBulk        = travelLen / pathLen;
        }

        // Check absorption on EBs and domain
        if (checkDom) {
          absorbedDomain = ParticleOps::domainIntersection(oldPos, newPos, probLo, probHi, sDomain);
          sDomain        = (absorbedDomain) ? std::max((Real)0.0, sDomain - SAFETY) : sDomain;
        }

        if (checkEB) {
          switch (m_intersectionEB) {
          case IntersectionEB::Raycast: {
            absorbedEB = ParticleOps::ebIntersectionRaycast(impFunc, oldPos, newPos, 1.E-3 * dx, sEB);

            break;
          }
          case IntersectionEB::Bisection: {
            absorbedEB = ParticleOps::ebIntersectionBisect(impFunc, oldPos, newPos, m_bisectStep, sEB);

            break;
          }
          default: {
            MayDay::Error("McPhoto::advancePhotonsTransient -- logic bust in eb intersection");

            break;
          }
          }
        }

        // True if the photon was absorbed (on anything) and false otherwise.
        const bool absorb = absorbedBulk || absorbedEB || absorbedDomain;

        if (!absorb) {
          p.position() = newPos;
        }
        else { // Determine scenarios -- if the photon was absorbed the smallest s-parameter takes precedence.
          const Real s = std::min(sBulk, std::min(sEB, sDomain));

          p.position() = oldPos + s * path;

          // Now check where it was actually absorbed
          if (absorbedBulk && sBulk < std::min(sEB, sDomain)) {
            bulkPhotons.transfer(lit);
          }
          else if (absorbedEB && sEB < std::min(sBulk, sDomain)) {
            ebPhotons.transfer(lit);
          }
          else if (absorbedDomain && sDomain < std::min(sBulk, sEB)) {
            domPhotons.transfer(lit);
          }
          else {
            MayDay::Error("McPhoto::advancePhotonsTransient - logic bust");
          }
        }
      }
    }
  }

  // Photons will have moved off the MPI ranks so we need to remap.
  a_bulkPhotons.remap();
  a_ebPhotons.remap();
  a_domainPhotons.remap();
  a_photons.remap();
}

void
McPhoto::remap()
{
  CH_TIME("McPhoto::remap()");
  if (m_verbosity > 5) {
    pout() << m_name + "::remap()" << endl;
  }

  this->remap(m_photons);
}

void
McPhoto::remap(ParticleContainer<Photon>& a_photons)
{
  CH_TIME("McPhoto::remap(ParticleContainer<Photon>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::remap(ParticleContainer<Photon>)" << endl;
  }

  a_photons.remap();
}

int
McPhoto::countPhotons(const AMRParticles<Photon>& a_photons) const
{
  CH_TIME("McPhoto::countPhotons");
  if (m_verbosity > 5) {
    pout() << m_name + "::countPhotons" << endl;
  }

  int numPhotons = 0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    numPhotons += a_photons[lvl]->numValid();
  }

  return numPhotons;
}

int
McPhoto::countOutcast(const AMRParticles<Photon>& a_photons) const
{
  CH_TIME("McPhoto::countOutcast");
  if (m_verbosity > 5) {
    pout() << m_name + "::countOutcast" << endl;
  }

  int num_outcast = 0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    num_outcast += a_photons[lvl]->numOutcast();
  }

  return num_outcast;
}

void
McPhoto::writePlotData(LevelData<EBCellFAB>& a_output,
                       int&                  a_comp,
                       const std::string     a_outputRealm,
                       const int             a_level) const noexcept
{
  CH_TIMERS("McPhoto::writePlotData");
  CH_TIMER("McPhoto::writePlotData::mesh_data", t1);
  CH_TIMER("McPhoto::writePlotData::particle_data", t2);
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  CH_START(t1);
  if (m_plotPhi) {
    this->writeData(a_output, a_comp, m_phi, a_outputRealm, a_level, false, true);
  }
  if (m_plotSource) {
    this->writeData(a_output, a_comp, m_source, a_outputRealm, a_level, false, true);
  }
  CH_STOP(t1);

  CH_START(t2);
  if (m_plotPhotons) {
    this->depositPhotonsNGP(*m_scratch[a_level], m_photons, a_level);

    const Interval srcInterv = Interval(0, 0);
    const Interval dstInterv = Interval(a_comp, a_comp);

    m_amr->copyData(a_output,
                    *m_scratch[a_level],
                    a_level,
                    a_outputRealm,
                    m_realm,
                    dstInterv,
                    srcInterv,
                    CopyStrategy::ValidGhost,
                    CopyStrategy::ValidGhost);

    a_comp++;
  }
  if (m_plotBulkPhotons) {
    this->depositPhotonsNGP(*m_scratch[a_level], m_bulkPhotons, a_level);

    const Interval srcInterv = Interval(0, 0);
    const Interval dstInterv = Interval(a_comp, a_comp);

    m_amr->copyData(a_output,
                    *m_scratch[a_level],
                    a_level,
                    a_outputRealm,
                    m_realm,
                    dstInterv,
                    srcInterv,
                    CopyStrategy::ValidGhost,
                    CopyStrategy::ValidGhost);

    a_comp++;
  }
  if (m_plotEBPhotons) {
    this->depositPhotonsNGP(*m_scratch[a_level], m_ebPhotons, a_level);

    const Interval srcInterv = Interval(0, 0);
    const Interval dstInterv = Interval(a_comp, a_comp);

    m_amr->copyData(a_output,
                    *m_scratch[a_level],
                    a_level,
                    a_outputRealm,
                    m_realm,
                    dstInterv,
                    srcInterv,
                    CopyStrategy::ValidGhost,
                    CopyStrategy::ValidGhost);

    a_comp++;
  }
  if (m_plotDomainPhotons) {
    this->depositPhotonsNGP(*m_scratch[a_level], m_domainPhotons, a_level);

    const Interval srcInterv = Interval(0, 0);
    const Interval dstInterv = Interval(a_comp, a_comp);

    m_amr->copyData(a_output,
                    *m_scratch[a_level],
                    a_level,
                    a_outputRealm,
                    m_realm,
                    dstInterv,
                    srcInterv,
                    CopyStrategy::ValidGhost,
                    CopyStrategy::ValidGhost);

    a_comp++;
  }
  if (m_plotSourcePhotons) {
    this->depositPhotonsNGP(*m_scratch[a_level], m_sourcePhotons, a_level);

    const Interval srcInterv = Interval(0, 0);
    const Interval dstInterv = Interval(a_comp, a_comp);

    m_amr->copyData(a_output,
                    *m_scratch[a_level],
                    a_level,
                    a_outputRealm,
                    m_realm,
                    dstInterv,
                    srcInterv,
                    CopyStrategy::ValidGhost,
                    CopyStrategy::ValidGhost);

    a_comp++;
  }
  CH_STOP(t2);
}

ParticleContainer<Photon>&
McPhoto::getPhotons()
{
  CH_TIME("McPhoto::getPhotons");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPhotons" << endl;
  }

  return m_photons;
}

ParticleContainer<Photon>&
McPhoto::getBulkPhotons()
{
  CH_TIME("McPhoto::getBulkPhotons");
  if (m_verbosity > 5) {
    pout() << m_name + "::getBulkPhotons" << endl;
  }

  return m_bulkPhotons;
}

ParticleContainer<Photon>&
McPhoto::getEbPhotons()
{
  CH_TIME("McPhoto::getEbPhotons");
  if (m_verbosity > 5) {
    pout() << m_name + "::getEbPhotons" << endl;
  }

  return m_ebPhotons;
}

ParticleContainer<Photon>&
McPhoto::getDomainPhotons()
{
  CH_TIME("McPhoto::getDomainPhotons");
  if (m_verbosity > 5) {
    pout() << m_name + "::getDomainPhotons" << endl;
  }

  return m_domainPhotons;
}

ParticleContainer<Photon>&
McPhoto::getSourcePhotons()
{
  CH_TIME("McPhoto::getSourcePhotons");
  if (m_verbosity > 5) {
    pout() << m_name + "::getSourcePhotons" << endl;
  }

  return m_sourcePhotons;
}

int
McPhoto::getMaxPhotonsPerCell() const noexcept
{
  CH_TIME("McPhoto::getMaxPhotonsPerCell");
  if (m_verbosity > 5) {
    pout() << m_name + "::getMaxPhotonsPerCell" << endl;
  }

  return (int)m_maxPhotonsGeneratedPerCell;
}

int
McPhoto::getNumSamplingPackets() const noexcept
{
  CH_TIME("McPhoto::getNumSamplingPackets");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumSamplingPackets" << endl;
  }

  return m_numSamplingPackets;
}

void
McPhoto::computeLoads(Vector<long long>& a_loads, const DisjointBoxLayout& a_dbl, const int a_level) const noexcept
{
  CH_TIME("McPhoto::computeLoads");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeLoads" << endl;
  }

  a_loads.resize(a_dbl.size(), 1LL);

  for (DataIterator dit(a_dbl); dit.ok(); ++dit) {
    a_loads[dit().intCode()] += m_photons[a_level][dit()].listItems().length();
    a_loads[dit().intCode()] += m_bulkPhotons[a_level][dit()].listItems().length();
  }
}

#include <CD_NamespaceFooter.H>
