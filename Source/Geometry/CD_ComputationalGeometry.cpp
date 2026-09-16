/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ComputationalGeometry.cpp
 * @brief  Implementation of CD_ComputationalGeometry.H
 * @author Robert Marskar
 */

// Std includes
#include <cmath>

// Chombo includes
#include <MFIndexSpace.H>
#include <IntersectionIF.H>
#include <UnionIF.H>
#include <AllRegularService.H>
#include <GeometryService.H>
#include <WrappedGShop.H>
#include <GeometryShop.H>
#include <ComplementIF.H>
#include <MayDay.H>

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NewIntersectionIF.H>
#include <CD_ScanShop.H>
#include <CD_PolyhedralGeometryShop.H>
#include <CD_MemoryReport.H>
#include <CD_NamespaceHeader.H>

ComputationalGeometry::ComputationalGeometry()
  : m_generator(Generator::GeometryShop),
    m_probLo(RealVect::Zero),
    m_eps0(1.0),
    m_refineAngle(0.0),
    m_maxGhostEB(0),
    m_startLevel(0),
    m_stopLevel(0),
    m_minBlockSize(0),
    m_maxBlockSize(0)
{
  CH_TIME("ComputationalGeometry::ComputationalGeometry()");

  // Default parameters.

  m_electrodes.resize(0);
  m_dielectrics.resize(0);

  m_scanDomain = ProblemDomain();

  m_multifluidIndexSpace = RefCountedPtr<MultiFluidIndexSpace>(new MultiFluidIndexSpace());
}

ComputationalGeometry::~ComputationalGeometry()
{
  CH_TIME("ComputationalGeometry::~ComputationalGeometry()");
}

void
ComputationalGeometry::useScanShop(const ProblemDomain& a_beginDomain)
{
  CH_TIME("ComputationalGeometry::useScanShop(ProblemDomain)");

  // TLDR: If you called this function you signal that ComputationalGeometry will use ScanShop for geometry generation.

  m_generator  = Generator::ScanShop;
  m_scanDomain = a_beginDomain;
}

void
ComputationalGeometry::usePolyhedralShop(const ProblemDomain& a_beginDomain)
{
  CH_TIME("ComputationalGeometry::usePolyhedralShop(ProblemDomain)");

  m_generator  = Generator::PolyhedralShop;
  m_scanDomain = a_beginDomain;
}

void
ComputationalGeometry::useChomboShop()
{
  CH_TIME("ComputationalGeometry::useChomboShop()");

  // TLDR: If you called this function you signal that ComputationalGeometry will use Chombo's GeometryShop for geometry
  // generation.
  m_generator  = Generator::GeometryShop;
  m_scanDomain = ProblemDomain();
}

const Vector<Dielectric>&
ComputationalGeometry::getDielectrics() const
{
  CH_TIME("ComputationalGeometry::getDielectrics()");

  return (m_dielectrics);
}

const Vector<Electrode>&
ComputationalGeometry::getElectrodes() const
{
  CH_TIME("ComputationalGeometry::getElectrodes()");

  return (m_electrodes);
}

const RefCountedPtr<BaseIF>&
ComputationalGeometry::getGasImplicitFunction() const
{
  CH_TIME("ComputationalGeometry::getGasImplicitFunction()");

  return (m_implicitFunctionGas);
}

const RefCountedPtr<BaseIF>&
ComputationalGeometry::getSolidImplicitFunction() const
{
  CH_TIME("ComputationalGeometry::getSolidImplicitFunction()");

  return (m_implicitFunctionSolid);
}

const RefCountedPtr<BaseIF>&
ComputationalGeometry::getImplicitFunction(const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::getImplicitFunction(phase::which_phase)");

  return (a_phase == phase::gas) ? m_implicitFunctionGas : m_implicitFunctionSolid;
}

Real
ComputationalGeometry::getGasPermittivity() const
{
  CH_TIME("ComputationalGeometry::getGasPermittivity()");

  return (m_eps0);
}

const RefCountedPtr<MultiFluidIndexSpace>&
ComputationalGeometry::getMfIndexSpace() const
{
  CH_TIME("ComputationalGeometry::getMfIndexSpace()");

  return (m_multifluidIndexSpace);
}

void
ComputationalGeometry::setDielectrics(const Vector<Dielectric>& a_dielectrics)
{
  CH_TIME("ComputationalGeometry::setDielectrics(Vector<Dielectric>)");

  m_dielectrics = a_dielectrics;
}

void
ComputationalGeometry::setElectrodes(const Vector<Electrode>& a_electrodes)
{
  CH_TIME("ComputationalGeometry::setElectrodes(Vector<Electrode>)");

  m_electrodes = a_electrodes;
}

void
ComputationalGeometry::setGasPermittivity(const Real a_eps0)
{
  CH_TIME("ComputationalGeometry::setGasPermittivity(Real)");

  m_eps0 = a_eps0;
}

void
ComputationalGeometry::buildGeometries(const ProblemDomain& a_finestDomain,
                                       const RealVect&      a_probLo,
                                       const Real           a_finestDx,
                                       const int            a_nCellMax,
                                       const int            a_maxGhostEB,
                                       const int            a_maxCoarsen)
{
  CH_TIME("ComputationalGeometry::buildGeometries(ProblemDomain, RealVect, Real, int, int, int)");

  // Set the default maximum number of EB ghosts that we will ever use. This is needed because ScanShop will look
  // through grown grid patches when it determines if a grid patch is irregular or not.
  m_maxGhostEB = a_maxGhostEB;

  // Build the composite implicit functions and then the GeometryService* objects which can be passed to Chombo.
  this->buildImplicitFunctions();

  Vector<GeometryService*> geoServices(2, nullptr);

  this->buildGasGeometry(geoServices[phase::gas], a_finestDomain, a_probLo, a_finestDx);
  this->buildSolidGeometry(geoServices[phase::solid], a_finestDomain, a_probLo, a_finestDx);

  // Define the multifluid index space.
  const bool useDistributedData = (m_generator != Generator::GeometryShop);

  m_multifluidIndexSpace->define(a_finestDomain.domainBox(), // Define MF
                                 a_probLo,
                                 a_finestDx,
                                 geoServices,
                                 useDistributedData,
                                 a_nCellMax,
                                 a_maxCoarsen);

  // Delete temps.
  for (int i = 0; i < 2; i++) {
    if (geoServices[i] != nullptr) {
      delete geoServices[i];
    }
  }
}

void
ComputationalGeometry::buildImplicitFunctions()
{
  CH_TIME("ComputationalGeometry::buildImplicitFunctions()");

  Vector<BaseIF*> dielectricParts;
  Vector<BaseIF*> electrodeParts;
  Vector<BaseIF*> allParts;

  for (int i = 0; i < m_dielectrics.size(); i++) {
    dielectricParts.push_back(&(*(m_dielectrics[i].getImplicitFunction())));
    allParts.push_back(&(*(m_dielectrics[i].getImplicitFunction())));
  }

  for (int i = 0; i < m_electrodes.size(); i++) {
    electrodeParts.push_back(&(*(m_electrodes[i].getImplicitFunction())));
    allParts.push_back(&(*(m_electrodes[i].getImplicitFunction())));
  }

  // The gas phase is the intersection of the region outside every object, so IntersectionIF is correct here.
  m_implicitFunctionGas = RefCountedPtr<BaseIF>(new NewIntersectionIF(allParts));

  // The solid phase is the region inside the dielectrics and outside the electrodes: the intersection of the
  // complement of "outside the dielectrics" with "outside the electrodes". There is none without dielectrics.
  if (dielectricParts.size() == 0) {
    m_implicitFunctionSolid = RefCountedPtr<BaseIF>();
  }
  else {
    Vector<BaseIF*> parts;

    RefCountedPtr<BaseIF> dielBaseIF = RefCountedPtr<BaseIF>(new NewIntersectionIF(dielectricParts));
    RefCountedPtr<BaseIF> elecBaseIF = RefCountedPtr<BaseIF>(new NewIntersectionIF(electrodeParts));
    RefCountedPtr<BaseIF> dielCompIF = RefCountedPtr<BaseIF>(new ComplementIF(*dielBaseIF));

    parts.push_back(&(*dielCompIF));
    parts.push_back(&(*elecBaseIF));

    m_implicitFunctionSolid = RefCountedPtr<BaseIF>(new IntersectionIF(parts));
  }
}

void
ComputationalGeometry::makeGrids(const ProblemDomain& a_startDomain,
                                 const ProblemDomain& a_stopDomain,
                                 const RealVect&      a_probLo,
                                 const Real           a_finestDx,
                                 const Real           a_refineAngle,
                                 const int            a_minBlockSize,
                                 const int            a_maxBlockSize,
                                 const int            a_maxGhostEB)
{
  CH_TIME("ComputationalGeometry::makeGrids");

  // Preconditions. Hard aborts rather than assertions: a violation here produces grids the index space would
  // serve silently wrong, and the cost of the check is nothing.
  if (a_startDomain.domainBox().isEmpty()) {
    MayDay::Error("ComputationalGeometry::makeGrids - the start domain is empty");
  }
  if (a_finestDx <= 0.0) {
    MayDay::Error("ComputationalGeometry::makeGrids - the finest grid spacing must be positive");
  }
  if (a_refineAngle < 0.0) {
    MayDay::Error("ComputationalGeometry::makeGrids - the refinement angle must not be negative");
  }
  if (a_maxGhostEB < 0) {
    MayDay::Error("ComputationalGeometry::makeGrids - the ghost width must not be negative");
  }
  if (a_minBlockSize <= 0) {
    MayDay::Error("ComputationalGeometry::makeGrids - the tile size must be positive");
  }
  if (a_maxBlockSize % a_minBlockSize != 0) {
    MayDay::Error("ComputationalGeometry::makeGrids - the super-tile size must be a multiple of the tile size");
  }

  // The one-tile nesting buffer between levels has to cover the ghost cells.
  if (a_minBlockSize < 2 * a_maxGhostEB) {
    MayDay::Error("ComputationalGeometry::makeGrids - the tile size must be at least twice the ghost width");
  }

  // The start domain and every level above it are tiled, so they must decompose into whole tiles. The levels
  // below the start domain are built whole and box by box, as ScanShop builds them, and may be smaller than a
  // tile; nothing is required of them.
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (a_startDomain.domainBox().size(dir) % a_minBlockSize != 0) {
      MayDay::Error("ComputationalGeometry::makeGrids - the start domain does not decompose into whole tiles");
    }
  }

  m_probLo       = a_probLo;
  m_minBlockSize = a_minBlockSize;
  m_maxBlockSize = a_maxBlockSize;
  m_maxGhostEB   = a_maxGhostEB;
  m_refineAngle  = a_refineAngle;

  // The levels are every factor-two coarsening of the start domain down to the coarsest domain that can still be
  // coarsened by two, the start domain itself, and every factor-two refinement of it up to the stop domain.
  // Level 0 is the coarsest domain.
  m_startLevel = 0;

  for (ProblemDomain coarDomain = a_startDomain; coarDomain.domainBox().coarsenable(2); coarDomain.coarsen(2)) {
    m_startLevel++;
  }

  m_stopLevel = m_startLevel;

  for (ProblemDomain fineDomain = a_startDomain; fineDomain.domainBox().size(0) < a_stopDomain.domainBox().size(0);
       fineDomain.refine(2)) {
    m_stopLevel++;
  }

  if (refine(a_startDomain, static_cast<int>(std::pow(2, m_stopLevel - m_startLevel))) != a_stopDomain) {
    MayDay::Error("ComputationalGeometry::makeGrids - the stop domain is not a refinement by two of the start domain");
  }

  const int numLevels = 1 + m_stopLevel;

  m_domains.resize(numLevels);
  m_dx.resize(numLevels);

  m_domains[m_stopLevel] = a_stopDomain;
  m_dx[m_stopLevel]      = a_finestDx;

  for (int lvl = m_stopLevel - 1; lvl >= 0; lvl--) {
    m_domains[lvl] = coarsen(m_domains[lvl + 1], 2);
    m_dx[lvl]      = 2.0 * m_dx[lvl + 1];
  }

  m_tiles.resize(numLevels);
  m_gasRegularBoxes.resize(numLevels);
  m_gasCoveredBoxes.resize(numLevels);
  m_gasIrregularBoxes.resize(numLevels);
  m_solidRegularBoxes.resize(numLevels);
  m_solidCoveredBoxes.resize(numLevels);
  m_solidIrregularBoxes.resize(numLevels);

  this->buildImplicitFunctions();

  // The algorithm, in the order it runs (the record of why it looks like this is PARTIAL_GRIDS.md):
  //
  //   0. Start level, per phase: domainSplit the whole domain and classify every box regular, covered or
  //      irregular (buildStartLevel, classifyBox).
  //   1. Upward, per phase, to the stop domain: a regular or covered box refines whole with its tag. An
  //      irregular box is split into classified pieces if the implicit function's normal turns by more than
  //      m_refineAngle between neighbouring cells near the surface; otherwise it is a leaf and nothing is built
  //      above it (buildFinerLevels, exceedsCurvature).
  //   2. Tiles, once, after both phases: the union of the two phases' irregular boxes on every level is tiled by
  //      TiledMeshRefine into a properly nested set common to both phases (makeTiles). This is the coverage the
  //      simulation regrids onto.
  //   3. Per phase: every tile lies inside an irregular box (it is irregular), inside a regular or covered box
  //      (it inherits, and the box is recorded as hit), or above a leaf of this phase (it is classified)
  //      (classifyTiles).
  //   4. Per phase: every hit regular or covered box is cut down to what the tiles left of it (decimateBoxes).
  //   5. The levels coarser than the start level, whole and classified box by box as ScanShop builds them; then,
  //      from the finest level down, a box containing a finer irregular box is irregular (buildCoarserLevels).
  //
  // Steps 0 and 1, per phase. A phase without an implicit function is one regular box on every level, and
  // takes no further part.
  const phase::which_phase phases[2] = {phase::gas, phase::solid};

  for (const phase::which_phase& curPhase : phases) {
    if (this->getImplicitFunction(curPhase).isNull()) {
      Vector<Vector<Box>>& regularBoxes = this->boxes(curPhase, GeometryService::Regular);

      for (int lvl = 0; lvl < numLevels; lvl++) {
        regularBoxes[lvl].push_back(m_domains[lvl].domainBox());
      }
    }
    else {
      this->buildStartLevel(curPhase);
      this->buildFinerLevels(curPhase);
    }
  }

  // Step 2, once.
  this->makeTiles();

  // Steps 3 and 4, per phase.
  for (const phase::which_phase& curPhase : phases) {
    Vector<Vector<GeometryService::InOut>> tileTypes(numLevels);
    Vector<Vector<int>>                    tileHosts(numLevels);

    this->classifyTiles(curPhase, tileTypes, tileHosts);
    this->decimateBoxes(curPhase, tileTypes, tileHosts);
  }

  // Step 5, once.
  this->buildCoarserLevels();
}

int
ComputationalGeometry::getNumGridLevels() const noexcept
{
  return m_domains.size();
}

int
ComputationalGeometry::getLevel(const ProblemDomain& a_domain) const noexcept
{
  int level = -1;

  for (int lvl = 0; lvl < m_domains.size(); lvl++) {
    if (m_domains[lvl].domainBox() == a_domain.domainBox()) {
      level = lvl;
    }
  }

  return level;
}

const Vector<Box>&
ComputationalGeometry::getBoxes(const phase::which_phase     a_phase,
                                const int                    a_level,
                                const GeometryService::InOut a_type) const noexcept
{
  return this->boxes(a_phase, a_type)[a_level];
}

Vector<Vector<Box>>&
ComputationalGeometry::boxes(const phase::which_phase a_phase, const GeometryService::InOut a_type) noexcept
{
  const auto& constThis = *this;

  return const_cast<Vector<Vector<Box>>&>(constThis.boxes(a_phase, a_type));
}

const Vector<Vector<Box>>&
ComputationalGeometry::boxes(const phase::which_phase a_phase, const GeometryService::InOut a_type) const noexcept
{
  const bool gas = (a_phase == phase::gas);

  switch (a_type) {
  case GeometryService::Regular: {
    return gas ? m_gasRegularBoxes : m_solidRegularBoxes;
  }
  case GeometryService::Covered: {
    return gas ? m_gasCoveredBoxes : m_solidCoveredBoxes;
  }
  default: {
    return gas ? m_gasIrregularBoxes : m_solidIrregularBoxes;
  }
  }
}

void
ComputationalGeometry::buildStartLevel(const phase::which_phase a_phase)
{
  CH_TIME("ComputationalGeometry::buildStartLevel");

  // Step 0. See the header for the pseudocode.
}

void
ComputationalGeometry::buildFinerLevels(const phase::which_phase a_phase)
{
  CH_TIME("ComputationalGeometry::buildFinerLevels");

  // Step 1. See the header for the pseudocode.
}

GeometryService::InOut
ComputationalGeometry::classifyBox(const Box& a_box, const int a_level, const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::classifyBox");

  // ScanShop::isRegular/isCovered on the implicit function. See the header for the pseudocode.
  return GeometryService::Irregular;
}

bool
ComputationalGeometry::exceedsCurvature(const Box& a_box, const int a_level, const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::exceedsCurvature");

  // The band test on finite-difference normals. See the header for the pseudocode.
  return false;
}

void
ComputationalGeometry::makeTiles()
{
  CH_TIME("ComputationalGeometry::makeTiles");

  // Step 2. See the header for the pseudocode.
}

void
ComputationalGeometry::classifyTiles(const phase::which_phase                a_phase,
                                     Vector<Vector<GeometryService::InOut>>& a_tileTypes,
                                     Vector<Vector<int>>&                    a_tileHosts) const
{
  CH_TIME("ComputationalGeometry::classifyTiles");

  // Step 3. See the header for the pseudocode.
}

void
ComputationalGeometry::decimateBoxes(const phase::which_phase                      a_phase,
                                     const Vector<Vector<GeometryService::InOut>>& a_tileTypes,
                                     const Vector<Vector<int>>&                    a_tileHosts)
{
  CH_TIME("ComputationalGeometry::decimateBoxes");

  // Step 4. See the header for the pseudocode.
}

void
ComputationalGeometry::buildCoarserLevels()
{
  CH_TIME("ComputationalGeometry::buildCoarserLevels");

  // Step 5. See the header for the pseudocode.
}

void
ComputationalGeometry::buildGasGeometry(GeometryService*&    a_geoserver,
                                        const ProblemDomain& a_finestDomain,
                                        const RealVect&      a_probLo,
                                        const Real           a_finestDx)
{
  CH_TIME("ComputationalGeometry::buildGasGeometry(GeometryService, ProblemDomain, RealVect, Real)");

  // Build the EBIS geometry. Use ScanShop, the polyhedral generator, or Chombo here.
  if (m_generator == Generator::PolyhedralShop) {
    auto* shop = new PolyhedralGeometryShop(*m_implicitFunctionGas,
                                            0,
                                            a_finestDx,
                                            a_probLo,
                                            a_finestDomain,
                                            m_scanDomain,
                                            m_maxGhostEB,
                                            s_thresh,
                                            s_strictGeometry);

    shop->setProfileFileName("PolyhedralShopReportGasPhase.dat");

    a_geoserver = static_cast<GeometryService*>(shop);
  }
  else if (m_generator == Generator::ScanShop) {
    auto* scanShop = new ScanShop(*m_implicitFunctionGas,
                                  0,
                                  a_finestDx,
                                  a_probLo,
                                  a_finestDomain,
                                  m_scanDomain,
                                  m_maxGhostEB,
                                  s_thresh);

    scanShop->setProfileFileName("ScanShopReportGasPhase.dat");

    a_geoserver = static_cast<GeometryService*>(scanShop);
  }
  else { // Chombo geometry generation
    a_geoserver = static_cast<GeometryService*>(
      new GeometryShop(*m_implicitFunctionGas, 0, a_finestDx * RealVect::Unit, s_thresh));
  }
}

void
ComputationalGeometry::buildSolidGeometry(GeometryService*&    a_geoserver,
                                          const ProblemDomain& a_finestDomain,
                                          const RealVect&      a_probLo,
                                          const Real           a_finestDx)
{
  CH_TIME("ComputationalGeometry::buildSolidGeometry(GeometryService, ProblemDomain, RealVect, Real)");

  // There is no solid phase without dielectrics.
  if (m_implicitFunctionSolid.isNull()) {
    a_geoserver = nullptr;
  }
  else {
    // Build the EBIS geometry. Use ScanShop, the polyhedral generator, or Chombo here.
    if (m_generator == Generator::PolyhedralShop) {
      auto* shop = new PolyhedralGeometryShop(*m_implicitFunctionSolid,
                                              0,
                                              a_finestDx,
                                              a_probLo,
                                              a_finestDomain,
                                              m_scanDomain,
                                              m_maxGhostEB,
                                              s_thresh,
                                              s_strictGeometry);

      shop->setProfileFileName("PolyhedralShopReportSolidPhase.dat");

      a_geoserver = static_cast<GeometryService*>(shop);
    }
    else if (m_generator == Generator::ScanShop) {
      auto* scanShop = new ScanShop(*m_implicitFunctionSolid,
                                    0,
                                    a_finestDx,
                                    a_probLo,
                                    a_finestDomain,
                                    m_scanDomain,
                                    m_maxGhostEB,
                                    s_thresh);

      scanShop->setProfileFileName("ScanShopReportSolidPhase.dat");

      a_geoserver = static_cast<GeometryService*>(scanShop);
    }
    else { // Chombo geometry generation
      a_geoserver = static_cast<GeometryService*>(
        new GeometryShop(*m_implicitFunctionSolid, 0, a_finestDx * RealVect::Unit, s_thresh));
    }
  }
}

#include <CD_NamespaceFooter.H>
