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

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NewIntersectionIF.H>
#include <CD_ScanShop.H>
#include <CD_PolyhedralGeometryShop.H>
#include <CD_MemoryReport.H>
#include <CD_NamespaceHeader.H>

ComputationalGeometry::ComputationalGeometry()
  : m_eps0(1.0),
    m_generator(Generator::GeometryShop),
    m_maxGhostEB(0),
    m_gridProbLo(RealVect::Zero),
    m_maxEbDepth(0),
    m_minBlockSize(0),
    m_maxBlockSize(0),
    m_refineAngle(0.0)
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
                                 const RealVect&      a_probLo,
                                 const Real           a_finestDx,
                                 const Real           a_refineAngle,
                                 const int            a_maxEbDepth,
                                 const int            a_minBlockSize,
                                 const int            a_maxBlockSize,
                                 const int            a_maxGhostEB)
{
  CH_TIME("ComputationalGeometry::makeGrids");

  // The design requires a_minBlockSize >= 2 * a_maxGhostEB (the one-tile nesting buffer must cover the ghost
  // cells) and a_maxBlockSize a multiple of a_minBlockSize. Neither is enforced yet, since nothing is built yet.
  CH_assert(a_maxEbDepth >= 0);

  m_gridProbLo   = a_probLo;
  m_maxEbDepth   = a_maxEbDepth;
  m_minBlockSize = a_minBlockSize;
  m_maxBlockSize = a_maxBlockSize;
  m_maxGhostEB   = a_maxGhostEB;
  m_refineAngle  = a_refineAngle;

  // Every level is a factor-two refinement of the start domain, and the finest one has the finest spacing.
  const int numLevels = 1 + m_maxEbDepth;

  m_gridDomains.resize(numLevels);
  m_gridDx.resize(numLevels);

  m_gridDomains[0] = a_startDomain;
  m_gridDx[0]      = a_finestDx * std::pow(2.0, m_maxEbDepth);

  for (int lvl = 1; lvl < numLevels; lvl++) {
    m_gridDomains[lvl] = refine(m_gridDomains[lvl - 1], 2);
    m_gridDx[lvl]      = 0.5 * m_gridDx[lvl - 1];
  }

  m_levelBoxes.resize(2);
  m_levelBoxTypes.resize(2);
  m_tileHost.resize(2);
  m_tileTypes.resize(2);
  m_grids.resize(2);
  m_gridTypes.resize(2);

  for (int p = 0; p < 2; p++) {
    m_levelBoxes[p].resize(numLevels);
    m_levelBoxTypes[p].resize(numLevels);
    m_tileHost[p].resize(numLevels);
    m_tileTypes[p].resize(numLevels);
    m_grids[p].resize(numLevels);
    m_gridTypes[p].resize(numLevels);
  }

  m_tiles.resize(numLevels);

  this->buildImplicitFunctions();

  // The algorithm, in the order it runs (the record of why it looks like this is PARTIAL_GRIDS.md):
  //
  //   0. Start level, per phase: domainSplit the whole domain and classify every box regular, covered or
  //      irregular (buildStartLevel, classifyBox).
  //   1. Upward, per phase, to m_maxEbDepth: a regular or covered box refines whole with its tag. An irregular
  //      box is split into classified pieces if the implicit function's normal turns by more than m_refineAngle
  //      between neighbouring cells near the surface; otherwise it is a leaf and nothing is built above it
  //      (buildFinerLevels, exceedsCurvature).
  //   2. Tiles, once, after both phases: the union of the two phases' irregular boxes on every level is tiled by
  //      TiledMeshRefine into a properly nested set common to both phases (makeTiles). This is the coverage the
  //      simulation regrids onto.
  //   3. Per phase: every tile lies inside an irregular box (it is irregular), inside a regular or covered box
  //      (it inherits, and the box is recorded as hit), or above a leaf of this phase (it is classified)
  //      (classifyTiles).
  //   4. Per phase: every hit regular or covered box is cut down to what the tiles left of it (decimateBoxes).
  //
  // Steps 0 and 1, per phase. A phase without an implicit function is one regular box on every level, and
  // takes no further part.
  for (int p = 0; p < 2; p++) {
    const phase::which_phase curPhase = static_cast<phase::which_phase>(p);

    if (this->getImplicitFunction(curPhase).isNull()) {
      for (int lvl = 0; lvl < numLevels; lvl++) {
        m_levelBoxes[p][lvl].push_back(m_gridDomains[lvl].domainBox());
        m_levelBoxTypes[p][lvl].push_back(GeometryService::Regular);
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
  for (int p = 0; p < 2; p++) {
    const phase::which_phase curPhase = static_cast<phase::which_phase>(p);

    this->classifyTiles(curPhase);
    this->decimateBoxes(curPhase);
  }
}

int
ComputationalGeometry::getNumGridLevels() const noexcept
{
  return m_gridDomains.size();
}

const Vector<Box>&
ComputationalGeometry::getGrids(const phase::which_phase a_phase, const int a_level) const noexcept
{
  return m_grids[a_phase][a_level];
}

const Vector<GeometryService::InOut>&
ComputationalGeometry::getGridTypes(const phase::which_phase a_phase, const int a_level) const noexcept
{
  return m_gridTypes[a_phase][a_level];
}

Vector<Box>
ComputationalGeometry::getBoxes(const phase::which_phase     a_phase,
                                const int                    a_level,
                                const GeometryService::InOut a_type) const noexcept
{
  Vector<Box> boxes;

  const Vector<Box>&                    levelBoxes = m_levelBoxes[a_phase][a_level];
  const Vector<GeometryService::InOut>& levelTypes = m_levelBoxTypes[a_phase][a_level];

  for (int i = 0; i < levelBoxes.size(); i++) {
    if (levelTypes[i] == a_type) {
      boxes.push_back(levelBoxes[i]);
    }
  }

  return boxes;
}

const Vector<Box>&
ComputationalGeometry::getTiles(const int a_level) const noexcept
{
  return m_tiles[a_level];
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
ComputationalGeometry::classifyTiles(const phase::which_phase a_phase)
{
  CH_TIME("ComputationalGeometry::classifyTiles");

  // Step 3. See the header for the pseudocode.
}

void
ComputationalGeometry::decimateBoxes(const phase::which_phase a_phase)
{
  CH_TIME("ComputationalGeometry::decimateBoxes");

  // Step 4. See the header for the pseudocode.
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
