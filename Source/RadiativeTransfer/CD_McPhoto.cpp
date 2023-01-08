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

#define MC_PHOTO_DEBUG 0

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

  // If we're doing an instantaneous solver, do a safety cleanout first of m_photons first (past photons have been absorbed on the mesh)
  if (m_instantaneous) {
    this->clear(m_photons);
  }

  // Generate new photons and add them to m_sourcePhotons.
  this->clear(m_sourcePhotons);
  this->generatePhotons(m_sourcePhotons, a_source, a_dt);
  m_photons.addParticles(m_sourcePhotons);

  // Advance photons either instantaneously or transiently.
  if (m_instantaneous) {
    this->advancePhotonsInstantaneous(m_bulkPhotons, m_ebPhotons, m_domainPhotons, m_photons);
  }
  else {
    this->advancePhotonsTransient(m_bulkPhotons, m_ebPhotons, m_domainPhotons, m_photons, a_dt);
    this->remap(m_photons);
  }

  // Deposit photons on the mesh.
  this->depositPhotons(a_phi, m_bulkPhotons, m_deposition);

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

  Real maxPhotons;
  pp.get("max_photons", maxPhotons);

  // = -1 => no restriction
  if (maxPhotons <= 0.0) {
    m_maxPhotonsGeneratedPerCell = std::numeric_limits<size_t>::max();
  }
  else {
    m_maxPhotonsGeneratedPerCell = round(maxPhotons);
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

  pp.get("plot_deposition", str);
  m_plotNumbers = false;
  if (str == "num") {
    m_plotDeposition = DepositionType::NGP;
    m_plotNumbers    = true;
  }
  else if (str == "ngp") {
    m_plotDeposition = DepositionType::NGP;
  }
  else if (str == "cic") {
    m_plotDeposition = DepositionType::CIC;
  }
  else {
    MayDay::Error("McPhoto::parseDeposition - unsupported interpolant requested");
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
McPhoto::allocateInternals()
{
  CH_TIME("McPhoto::allocateInternals");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocateInternals" << endl;
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
}

void
McPhoto::deallocateInternals()
{
  CH_TIME("McPhoto::deallocateInternals");
  if (m_verbosity > 5) {
    pout() << m_name + "::deallocateInternals" << endl;
  }
}

void
McPhoto::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("McPhoto::regrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::regrid" << endl;
  }

  // Mesh data regrids
  m_amr->reallocate(m_phi, m_phase, a_lmin);
  m_amr->reallocate(m_source, m_phase, a_lmin);
  m_amr->reallocate(m_scratch, m_phase, a_lmin);
  m_amr->reallocate(m_depositionNC, m_phase, a_lmin);
  m_amr->reallocate(m_massDiff, m_phase, a_lmin);

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
    plotVarNames.push_back(m_name + " photons");
  }
  if (m_plotBulkPhotons) {
    plotVarNames.push_back(m_name + " bulkPhotons");
  }
  if (m_plotEBPhotons) {
    plotVarNames.push_back(m_name + " ebPhotons");
  }
  if (m_plotDomainPhotons) {
    plotVarNames.push_back(m_name + " domainPhotons");
  }
  if (m_plotSourcePhotons) {
    plotVarNames.push_back(m_name + " sourcePhotons");
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
McPhoto::generatePhotons(ParticleContainer<Photon>& a_photons, const EBAMRCellData& a_source, const Real a_dt)
{
  CH_TIME("McPhoto::generatePhotons");
  if (m_verbosity > 5) {
    pout() << m_name + "::generatePhotons" << endl;
  }

  CH_assert(a_source[0]->nComp() == 1);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const RealVect           probLo = m_amr->getProbLo();
    const Real               dx     = m_amr->getDx()[lvl];
    const Real               vol    = pow(dx, SpaceDim);

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box      cellBox = dbl.get(dit());
      const EBISBox& ebisbox = (*a_source[lvl])[dit()].getEBISBox();

      const EBCellFAB&     source     = (*a_source[lvl])[dit()];
      const FArrayBox&     srcFAB     = source.getFArrayBox();
      const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_realm)[lvl])[dit()];

      List<Photon>& photons = a_photons[lvl][dit()].listItems();
      photons.clear();

      // Regular cells. Note that we make superphotons if we have to. Also, only draw photons in valid cells,
      // grids that are covered by finer grids don't draw photons.
      auto regularKernel = [&](const IntVect& iv) -> void {
        if (ebisbox.isRegular(iv) && validCells(iv, 0)) {

          const size_t numPhysPhotons = this->drawPhotons(srcFAB(iv, 0), vol, a_dt);

          if (numPhysPhotons > 0) {
            const std::vector<size_t> photonWeights =
              ParticleManagement::partitionParticleWeights(numPhysPhotons, m_maxPhotonsGeneratedPerCell);

            const RealVect lo = probLo + RealVect(iv) * dx;
            const RealVect hi = lo + RealVect::Unit * dx;

            for (size_t i = 0; i < photonWeights.size(); i++) {

              // Determine starting position within cell, propagation direction, absorption
              // length, and weight.
              const RealVect pos    = Random::randomPosition(lo, hi);
              const RealVect v      = Units::c * Random::getDirection();
              const Real     weight = (Real)photonWeights[i];
              const Real     kappa  = m_rtSpecies->getAbsorptionCoefficient(pos);

              photons.add(Photon(pos, v, kappa, weight));
            }
          }
        }
      };

      // Irregular kernel. Same as the above really.
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();

        if (validCells(iv, 0)) {

          const size_t numPhysPhotons = this->drawPhotons(source(vof, 0), vol, a_dt);

          if (numPhysPhotons > 0) {
            const std::vector<size_t> photonWeights =
              ParticleManagement::partitionParticleWeights(numPhysPhotons, m_maxPhotonsGeneratedPerCell);

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
              const Real     kappa  = m_rtSpecies->getAbsorptionCoefficient(pos);
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

size_t
McPhoto::drawPhotons(const Real a_source, const Real a_volume, const Real a_dt)
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

  // Increment level redistribution register
  this->incrementRedist(m_massDiff);

  // Do the redistribution magic
  this->coarseFineIncrement(m_massDiff); // Compute C2F, F2C, and C2C mass transfers
  this->levelRedist(a_phi);              // Level redistribution. Weights is a dummy parameter
  this->coarseFineRedistribution(a_phi); // Do the coarse-fine redistribution

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
McPhoto::incrementRedist(const EBAMRIVData& a_massDifference)
{
  CH_TIME("McPhoto::incrementRedist");
  if (m_verbosity > 5) {
    pout() << m_name + "::incrementRedist" << endl;
  }

  const Interval interv(m_comp, m_comp);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    EBLevelRedist& levelRedist = *(m_amr->getLevelRedist(m_realm, m_phase)[lvl]);
    levelRedist.setToZero();

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      levelRedist.increment((*a_massDifference[lvl])[dit()], dit(), interv);
    }
  }
}

void
McPhoto::levelRedist(EBAMRCellData& a_phi)
{
  CH_TIME("McPhoto::levelRedist");
  if (m_verbosity > 5) {
    pout() << m_name + "::levelRedist" << endl;
  }

  const Interval interv(m_comp, m_comp);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    EBLevelRedist& levelRedist = *(m_amr->getLevelRedist(m_realm, m_phase)[lvl]);
    levelRedist.redistribute(*a_phi[lvl], interv);
    levelRedist.setToZero();
  }
}

void
McPhoto::coarseFineIncrement(const EBAMRIVData& a_massDifference)
{
  CH_TIME("McPhoto::coarseFineIncrement");
  if (m_verbosity > 5) {
    pout() << m_name + "::coarseFineIncrement" << endl;
  }

  const Interval interv(m_comp, m_comp);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    RefCountedPtr<EBFineToCoarRedist>& fine2coarRedist = m_amr->getFineToCoarRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToFineRedist>& coar2fineRedist = m_amr->getCoarToFineRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coarRedist = m_amr->getCoarToCoarRedist(m_realm, m_phase)[lvl];

    const bool hasCoar = lvl > 0;
    const bool hasFine = lvl < 0;

    if (hasCoar) {
      fine2coarRedist->setToZero();
    }
    if (hasFine) {
      coar2fineRedist->setToZero();
      coar2coarRedist->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      if (hasCoar) {
        fine2coarRedist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
      }

      if (hasFine) {
        coar2fineRedist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
        coar2coarRedist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
      }
    }
  }
}

void
McPhoto::coarseFineRedistribution(EBAMRCellData& a_phi)
{
  CH_TIME("McPhoto::coarseFineRedistribution");
  if (m_verbosity > 5) {
    pout() << m_name + "::coarseFineRedistribution" << endl;
  }

  const Interval interv(m_comp, m_comp);
  const int      finestLevel = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const bool hasCoar = lvl > 0;
    const bool hasFine = lvl < finestLevel;

    RefCountedPtr<EBCoarToFineRedist>& coar2fineRedist = m_amr->getCoarToFineRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coarRedist = m_amr->getCoarToCoarRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBFineToCoarRedist>& fine2coarRedist = m_amr->getFineToCoarRedist(m_realm, m_phase)[lvl];

    if (hasCoar) {
      fine2coarRedist->redistribute(*a_phi[lvl - 1], interv);
      fine2coarRedist->setToZero();
    }

    if (hasFine) {
      coar2fineRedist->redistribute(*a_phi[lvl + 1], interv);
      coar2coarRedist->redistribute(*a_phi[lvl], interv);

      coar2fineRedist->setToZero();
      coar2coarRedist->setToZero();
    }
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

#if MC_PHOTO_DEBUG // Debug
  const int photonsBefore = this->countPhotons(a_photons.getParticles());
#endif

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const Real               dx  = m_amr->getDx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      List<Photon>& bulkPhotons = a_bulkPhotons[lvl][dit()].listItems();
      List<Photon>& ebPhotons   = a_ebPhotons[lvl][dit()].listItems();
      List<Photon>& domPhotons  = a_domainPhotons[lvl][dit()].listItems();
      List<Photon>& allPhotons  = a_photons[lvl][dit()].listItems();

      // These must be cleared -- I assume that users will already have used any information
      // from a past time step already.
      bulkPhotons.clear();
      ebPhotons.clear();
      domPhotons.clear();

      // Iterate over the Photons that will be moved.
      for (ListIterator<Photon> lit(allPhotons); lit.ok(); ++lit) {
        Photon& p = lit();

        // Draw a new random absorption position
        const RealVect oldPos    = p.position();
        const RealVect direction = p.velocity() / (p.velocity().vectorLength());
        const RealVect newPos    = oldPos + direction * this->randomExponential(p.kappa());
        const RealVect path      = newPos - oldPos;

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

        if (!checkEB && !checkDom) { // No intersection test necessary, photons is guaranteed to end up on the mesh.
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

          // Move the Photon to the data holder where it belongs.
          if (!contactEB && !contactDomain) {
            p.position() = newPos;
            bulkPhotons.add(p);
          }
          else {
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

#if MC_PHOTO_DEBUG // Debug hook
  const int bulkPhotons = this->countPhotons(a_bulkPhotons.getParticles());
  const int ebPhotons   = this->countPhotons(a_ebPhotons.getParticles());
  const int domPhotons  = this->countPhotons(a_domainPhotons.getParticles());

  if (procID() == 0) {
    std::cout << "Photons before = " << photonsBefore << "\n"
              << "bulk Photons = " << bulkPhotons << "\n"
              << "eb Photons = " << ebPhotons << "\n"
              << "dom Photons = " << domPhotons << "\n"
              << "Photons after = " << domPhotons + ebPhotons + bulkPhotons << "\n"
              << std::endl;
  }
#endif
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

#if MC_PHOTO_DEBUG // Debug hook.
  const int photonsBefore = this->countPhotons(a_photons.getParticles());
#endif

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const Real               dx  = m_amr->getDx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      List<Photon>& bulkPhotons = a_bulkPhotons[lvl][dit()].listItems();
      List<Photon>& ebPhotons   = a_ebPhotons[lvl][dit()].listItems();
      List<Photon>& domPhotons  = a_domainPhotons[lvl][dit()].listItems();
      List<Photon>& allPhotons  = a_photons[lvl][dit()].listItems();

      // These must be cleared -- I assume that information from past time steps has already been processed.
      bulkPhotons.clear();
      ebPhotons.clear();
      domPhotons.clear();

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

        const bool absorb = absorbedBulk || absorbedEB ||
                            absorbedDomain; // True if the photon was absorbed (on anything) and false otherwise.

        if (!absorb) {
          p.position() = newPos;
        }
        else { // Determine scenarios -- if the photon was absorbed the smallest s-parameter takes precedence.
          const Real s = Min(sBulk, Min(sEB, sDomain));
          p.position() = oldPos + s * path;

          // Now check where it was actually absorbed
          if (absorbedBulk && sBulk < Min(sEB, sDomain)) {
            bulkPhotons.transfer(lit);
          }
          else if (absorbedEB && sEB < Min(sBulk, sDomain)) {
            ebPhotons.transfer(lit);
          }
          else if (absorbedDomain && sDomain < Min(sBulk, sEB)) {
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

#if MC_PHOTO_DEBUG // Debugging hook
  const int bulkPhotons  = this->countPhotons(a_bulkPhotons.getParticles());
  const int ebPhotons    = this->countPhotons(a_ebPhotons.getParticles());
  const int domPhotons   = this->countPhotons(a_domainPhotons.getParticles());
  const int afterPhotons = this->countPhotons(a_photons.getParticles());

  if (procID() == 0) {
    std::cout << "Photons before = " << photonsBefore << "\n"
              << "Photons after = " << afterPhotons << "\n"
              << "bulk Photons = " << bulkPhotons << "\n"
              << "eb Photons = " << ebPhotons << "\n"
              << "dom Photons = " << domPhotons << "\n"
              << "total = " << domPhotons + ebPhotons + bulkPhotons + afterPhotons << "\n"
              << std::endl;
  }
#endif
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
McPhoto::writePlotData(EBAMRCellData& a_output, int& a_comp)
{
  CH_TIME("McPhoto::writePlotData");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }

  if (m_plotPhi) {
    this->writeData(a_output, a_comp, m_phi, false);
  }
  if (m_plotSource) {
    this->writeData(a_output, a_comp, m_source, false);
  }
  if (m_plotPhotons) {
    this->depositPhotons(m_scratch, m_photons, m_plotDeposition);
    this->writeData(a_output, a_comp, m_scratch, false);
  }
  if (m_plotBulkPhotons) {
    this->depositPhotons(m_scratch, m_bulkPhotons, m_plotDeposition);
    this->writeData(a_output, a_comp, m_scratch, false);
  }
  if (m_plotEBPhotons) {
    this->depositPhotons(m_scratch, m_ebPhotons, m_plotDeposition);
    this->writeData(a_output, a_comp, m_scratch, false);
  }
  if (m_plotDomainPhotons) {
    this->depositPhotons(m_scratch, m_domainPhotons, m_plotDeposition);
    this->writeData(a_output, a_comp, m_scratch, false);
  }
  if (m_plotSourcePhotons) {
    this->depositPhotons(m_scratch, m_sourcePhotons, m_plotDeposition);
    this->writeData(a_output, a_comp, m_scratch, false);
  }
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

#include <CD_NamespaceFooter.H>
