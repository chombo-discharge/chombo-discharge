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
#include <BaseFab.H>
#include <BoxIterator.H>
#include <BRMeshRefine.H>
#include <IntVectSet.H>
#include <SPMD.H>
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
#include <CD_ParallelOps.H>
#include <CD_TiledMeshRefine.H>
#include <CD_Units.H>
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
                                 const Real           a_startDx,
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
  if (a_startDx <= 0.0) {
    MayDay::Error("ComputationalGeometry::makeGrids - the grid spacing must be positive");
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

  m_domains[m_startLevel] = a_startDomain;
  m_dx[m_startLevel]      = a_startDx;

  for (int lvl = m_startLevel - 1; lvl >= 0; lvl--) {
    m_domains[lvl] = coarsen(m_domains[lvl + 1], 2);
    m_dx[lvl]      = 2.0 * m_dx[lvl + 1];
  }

  for (int lvl = m_startLevel + 1; lvl <= m_stopLevel; lvl++) {
    m_domains[lvl] = refine(m_domains[lvl - 1], 2);
    m_dx[lvl]      = 0.5 * m_dx[lvl - 1];
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

  // The start domain decomposes into whole tiles (checked in makeGrids), so the block factor is the tile: every
  // box is a whole number of tiles, and at most a super-tile wide.
  Vector<Box> boxes;

  domainSplit(m_domains[m_startLevel], boxes, m_maxBlockSize, m_minBlockSize);

  const Vector<GeometryService::InOut> types = this->classifyBoxes(boxes, m_startLevel, a_phase);

  for (int i = 0; i < boxes.size(); i++) {
    this->boxes(a_phase, types[i])[m_startLevel].push_back(boxes[i]);
  }
}

void
ComputationalGeometry::buildFinerLevels(const phase::which_phase a_phase)
{
  CH_TIME("ComputationalGeometry::buildFinerLevels");

  Vector<Vector<Box>>& regularBoxes   = this->boxes(a_phase, GeometryService::Regular);
  Vector<Vector<Box>>& coveredBoxes   = this->boxes(a_phase, GeometryService::Covered);
  Vector<Vector<Box>>& irregularBoxes = this->boxes(a_phase, GeometryService::Irregular);

  for (int lvl = m_startLevel; lvl < m_stopLevel; lvl++) {

    // Regular and covered boxes refine whole, keeping their tag. They go in first, in their parents' order, so
    // that a parent's refinement sits at the parent's index on the next level.
    for (int i = 0; i < regularBoxes[lvl].size(); i++) {
      regularBoxes[lvl + 1].push_back(refine(regularBoxes[lvl][i], 2));
    }

    for (int i = 0; i < coveredBoxes[lvl].size(); i++) {
      coveredBoxes[lvl + 1].push_back(refine(coveredBoxes[lvl][i], 2));
    }

    // An irregular box splits if the surface inside it turns too sharply for this level; otherwise it is a leaf
    // and nothing is built above it. The pieces of every box that splits are classified together so that the
    // work is shared once per level rather than once per box.
    const Vector<int> flags = this->splitFlags(irregularBoxes[lvl], lvl, a_phase);

    Vector<Box> pieces;

    for (int i = 0; i < irregularBoxes[lvl].size(); i++) {
      if (flags[i] != 0) {
        Vector<Box> split;

        domainSplit(refine(irregularBoxes[lvl][i], 2), split, m_maxBlockSize, m_minBlockSize);

        pieces.append(split);
      }
    }

    const Vector<GeometryService::InOut> types = this->classifyBoxes(pieces, lvl + 1, a_phase);

    for (int i = 0; i < pieces.size(); i++) {
      this->boxes(a_phase, types[i])[lvl + 1].push_back(pieces[i]);
    }
  }
}

Vector<GeometryService::InOut>
ComputationalGeometry::classifyBoxes(const Vector<Box>&       a_boxes,
                                     const int                a_level,
                                     const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::classifyBoxes");

  // The classification travels as an integer so that one all-reduce assembles it: a rank writes only the
  // entries it owns and leaves the rest at zero, and the sum is the union.
  constexpr int regular   = 1;
  constexpr int covered   = 2;
  constexpr int irregular = 3;

  Vector<int> codes(a_boxes.size(), 0);

  for (int i = procID(); i < a_boxes.size(); i += numProc()) {
    switch (this->classifyBox(a_boxes[i], a_level, a_phase)) {
    case GeometryService::Regular: {
      codes[i] = regular;

      break;
    }
    case GeometryService::Covered: {
      codes[i] = covered;

      break;
    }
    default: {
      codes[i] = irregular;

      break;
    }
    }
  }

  ParallelOps::sum(codes);

  Vector<GeometryService::InOut> types(a_boxes.size(), GeometryService::Irregular);

  for (int i = 0; i < a_boxes.size(); i++) {
    switch (codes[i]) {
    case regular: {
      types[i] = GeometryService::Regular;

      break;
    }
    case covered: {
      types[i] = GeometryService::Covered;

      break;
    }
    case irregular: {
      types[i] = GeometryService::Irregular;

      break;
    }
    default: {
      MayDay::Error("ComputationalGeometry::classifyBoxes - a box was classified by no rank or by several");

      break;
    }
    }
  }

  return types;
}

Vector<int>
ComputationalGeometry::splitFlags(const Vector<Box>& a_boxes, const int a_level, const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::splitFlags");

  Vector<int> flags(a_boxes.size(), 0);

  for (int i = procID(); i < a_boxes.size(); i += numProc()) {
    flags[i] = this->exceedsCurvature(a_boxes[i], a_level, a_phase) ? 1 : 0;
  }

  ParallelOps::sum(flags);

  return flags;
}

GeometryService::InOut
ComputationalGeometry::classifyBox(const Box& a_box, const int a_level, const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::classifyBox");

  // ScanShop::isRegular/isCovered on the implicit function. A cell centre within half a cell diagonal of the
  // zero set may belong to a cut cell, and one such cell makes the box irregular. Otherwise every centre is on
  // one side, and that side is the classification; both sides without a cell in between is not possible for a
  // continuous function and is reported rather than resolved.
  const BaseIF& f = *(this->getImplicitFunction(a_phase));

  const Real dx           = m_dx[a_level];
  const Real halfDiagonal = 0.5 * dx * std::sqrt(static_cast<Real>(SpaceDim));
  const Box  grown        = grow(a_box, m_maxGhostEB) & m_domains[a_level].domainBox();

  bool anyFluid = false;
  bool anySolid = false;

  for (BoxIterator bit(grown); bit.ok(); ++bit) {
    const IntVect iv = bit();

    RealVect x = m_probLo;

    for (int dir = 0; dir < SpaceDim; dir++) {
      x[dir] += dx * (static_cast<Real>(iv[dir]) + 0.5);
    }

    const Real value = f.value(x);

    if (std::abs(value) <= halfDiagonal) {
      return GeometryService::Irregular;
    }

    anyFluid = anyFluid || (value < 0.0);
    anySolid = anySolid || (value > 0.0);
  }

  if (anyFluid && anySolid) {
    MayDay::Error("ComputationalGeometry::classifyBox - a box holds fluid and solid but no cell near the surface");
  }

  return anySolid ? GeometryService::Covered : GeometryService::Regular;
}

bool
ComputationalGeometry::exceedsCurvature(const Box& a_box, const int a_level, const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::exceedsCurvature");

  // Normals are taken on the cells near the zero set, one cell beyond the box so that a pair across the box
  // boundary is seen from both sides, by central differences of the implicit function. The first pair of
  // neighbouring band cells whose normals differ by more than the refinement angle decides.
  const BaseIF& f = *(this->getImplicitFunction(a_phase));

  const Real dx    = m_dx[a_level];
  const Real band  = dx * std::sqrt(static_cast<Real>(SpaceDim));
  const Real h     = 0.5 * dx;
  const Box  valid = a_box & m_domains[a_level].domainBox();
  const Box  grown = grow(a_box, 1) & m_domains[a_level].domainBox();

  BaseFab<Real> normal(grown, SpaceDim);
  BaseFab<int>  inBand(grown, 1);

  inBand.setVal(0);

  for (BoxIterator bit(grown); bit.ok(); ++bit) {
    const IntVect iv = bit();

    RealVect x = m_probLo;

    for (int dir = 0; dir < SpaceDim; dir++) {
      x[dir] += dx * (static_cast<Real>(iv[dir]) + 0.5);
    }

    if (std::abs(f.value(x)) <= band) {
      RealVect n = RealVect::Zero;

      for (int dir = 0; dir < SpaceDim; dir++) {
        RealVect xHi = x;
        RealVect xLo = x;

        xHi[dir] += h;
        xLo[dir] -= h;

        n[dir] = f.value(xHi) - f.value(xLo);
      }

      const Real length = n.vectorLength();

      if (length > 0.0) {
        n /= length;

        inBand(iv, 0) = 1;

        for (int dir = 0; dir < SpaceDim; dir++) {
          normal(iv, dir) = n[dir];
        }
      }
    }
  }

  const Real cosThreshold = std::cos(m_refineAngle * Units::pi / 180.0);

  for (BoxIterator bit(valid); bit.ok(); ++bit) {
    const IntVect iv = bit();

    if (inBand(iv, 0) == 0) {
      continue;
    }

    const Box neighbours = grow(Box(iv, iv), 1) & grown;

    for (BoxIterator nit(neighbours); nit.ok(); ++nit) {
      const IntVect jv = nit();

      if (jv == iv || inBand(jv, 0) == 0) {
        continue;
      }

      Real dot = 0.0;

      for (int dir = 0; dir < SpaceDim; dir++) {
        dot += normal(iv, dir) * normal(jv, dir);
      }

      if (dot < cosThreshold) {
        return true;
      }
    }
  }

  return false;
}

void
ComputationalGeometry::makeTiles()
{
  CH_TIME("ComputationalGeometry::makeTiles");

  const int numAbove = m_stopLevel - m_startLevel;

  if (numAbove == 0) {
    return;
  }

  // TiledMeshRefine tiles level k from tags on level k - 1 and takes the start domain as its level 0, so the
  // irregular boxes of builder level lvl enter, coarsened by two, as tags for tiler level lvl - m_startLevel.
  // Tags are rank-local and the tiler gathers them, so each rank tags only its share of the boxes.
  const phase::which_phase phases[2] = {phase::gas, phase::solid};

  Vector<IntVectSet> tags(numAbove);

  for (int lvl = m_startLevel + 1; lvl <= m_stopLevel; lvl++) {
    IntVectSet& levelTags = tags[lvl - m_startLevel - 1];

    for (const phase::which_phase& curPhase : phases) {
      const Vector<Box>& irregularBoxes = this->boxes(curPhase, GeometryService::Irregular)[lvl];

      for (int i = procID(); i < irregularBoxes.size(); i += numProc()) {
        levelTags |= coarsen(irregularBoxes[i], 2);
      }
    }
  }

  const Vector<int> refRatios(1 + numAbove, 2);

  TiledMeshRefine tiler(m_domains[m_startLevel],
                        refRatios,
                        m_minBlockSize * IntVect::Unit,
                        m_maxBlockSize * IntVect::Unit);

  Vector<Vector<Box>> tiles;

  const int finestTiled = tiler.regrid(tiles, tags);

  // Tiler level 0 is the start domain, whole and not tiled, and is discarded as AmrMesh discards it.
  for (int lvl = m_startLevel + 1; lvl <= m_startLevel + finestTiled; lvl++) {
    m_tiles[lvl] = tiles[lvl - m_startLevel];
  }
}

void
ComputationalGeometry::classifyTiles(const phase::which_phase                a_phase,
                                     Vector<Vector<GeometryService::InOut>>& a_tileTypes,
                                     Vector<Vector<int>>&                    a_tileHosts) const
{
  CH_TIME("ComputationalGeometry::classifyTiles");

  // Step 3. Waits for the containing-box lookup; see the header for the pseudocode.
}

void
ComputationalGeometry::decimateBoxes(const phase::which_phase                      a_phase,
                                     const Vector<Vector<GeometryService::InOut>>& a_tileTypes,
                                     const Vector<Vector<int>>&                    a_tileHosts)
{
  CH_TIME("ComputationalGeometry::decimateBoxes");

  // Step 4. Waits for step 3; see the header for the pseudocode.
}

void
ComputationalGeometry::buildCoarserLevels()
{
  CH_TIME("ComputationalGeometry::buildCoarserLevels");

  // The levels coarser than the start domain are built as ScanShop builds them: whole, each box classified on
  // its own. No block factor, since these levels are not tiled and the coarsest may be smaller than a tile. A
  // phase without an implicit function already has its one regular box on every level.
  const phase::which_phase phases[2] = {phase::gas, phase::solid};

  for (const phase::which_phase& curPhase : phases) {
    if (this->getImplicitFunction(curPhase).isNull()) {
      continue;
    }

    for (int lvl = m_startLevel - 1; lvl >= 0; lvl--) {
      Vector<Box> boxes;

      domainSplit(m_domains[lvl], boxes, m_maxBlockSize);

      const Vector<GeometryService::InOut> types = this->classifyBoxes(boxes, lvl, curPhase);

      for (int i = 0; i < boxes.size(); i++) {
        this->boxes(curPhase, types[i])[lvl].push_back(boxes[i]);
      }
    }
  }

  // The push-down -- a box containing a finer irregular box is irregular -- waits for the containing-box
  // lookup; see the header for the pseudocode.
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
