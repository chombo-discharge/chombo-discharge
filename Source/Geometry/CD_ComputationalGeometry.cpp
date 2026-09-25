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
#include <algorithm>
#include <limits>
#include <cmath>

// Chombo includes
#include <BaseFab.H>
#include <BoxIterator.H>
#include <BRMeshRefine.H>
#include <IntVectSet.H>
#include <TreeIntVectSet.H>
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
#include <ParmParse.H>

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NewIntersectionIF.H>
#include <CD_ParallelOps.H>
#include <CD_TiledMeshRefine.H>
#include <CD_Timer.H>
#include <CD_Units.H>
#include <CD_ScanShop.H>
#include <CD_PolyhedralGeometryShop.H>
#include <CD_CutCellBody.H>
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
    m_minBlockSize(8),
    m_maxBlockSize(8),
    m_profile(false),
    m_verbose(false)
{
  CH_TIME("ComputationalGeometry::ComputationalGeometry()");

  // Default parameters.

  ParmParse pp("ComputationalGeometry");

  pp.query("verbose", m_verbose);

  if (m_verbose) {
    pout() << "ComputationalGeometry::ComputationalGeometry()" << endl;
  }

  // The tile and super-tile the grids are built with are the geometry's own, separate from the simulation's
  // block sizes: curvature refinement follows the surface, and a small tile keeps the refined footprint close
  // to it.
  pp.query("min_block_size", m_minBlockSize);
  pp.query("max_block_size", m_maxBlockSize);
  pp.query("profile", m_profile);

  m_electrodes.resize(0);
  m_dielectrics.resize(0);

  m_scanDomain = ProblemDomain();

  m_multifluidIndexSpace = RefCountedPtr<MultiFluidIndexSpace>(new MultiFluidIndexSpace());
}

ComputationalGeometry::~ComputationalGeometry()
{
  CH_TIME("ComputationalGeometry::~ComputationalGeometry()");
  if (m_verbose) {
    pout() << "ComputationalGeometry::~ComputationalGeometry()" << endl;
  }
}

void
ComputationalGeometry::useScanShop(const ProblemDomain& a_beginDomain)
{
  CH_TIME("ComputationalGeometry::useScanShop(ProblemDomain)");
  if (m_verbose) {
    pout() << "ComputationalGeometry::useScanShop(ProblemDomain)" << endl;
  }

  // TLDR: If you called this function you signal that ComputationalGeometry will use ScanShop for geometry generation.

  m_generator  = Generator::ScanShop;
  m_scanDomain = a_beginDomain;
}

void
ComputationalGeometry::usePolyhedralShop(const ProblemDomain& a_beginDomain)
{
  CH_TIME("ComputationalGeometry::usePolyhedralShop(ProblemDomain)");
  if (m_verbose) {
    pout() << "ComputationalGeometry::usePolyhedralShop(ProblemDomain)" << endl;
  }

  m_generator  = Generator::PolyhedralShop;
  m_scanDomain = a_beginDomain;
}

void
ComputationalGeometry::useChomboShop()
{
  CH_TIME("ComputationalGeometry::useChomboShop()");
  if (m_verbose) {
    pout() << "ComputationalGeometry::useChomboShop()" << endl;
  }

  // TLDR: If you called this function you signal that ComputationalGeometry will use Chombo's GeometryShop for geometry
  // generation.
  m_generator  = Generator::GeometryShop;
  m_scanDomain = ProblemDomain();
}

const Vector<Dielectric>&
ComputationalGeometry::getDielectrics() const
{
  CH_TIME("ComputationalGeometry::getDielectrics()");
  if (m_verbose) {
    pout() << "ComputationalGeometry::getDielectrics()" << endl;
  }

  return (m_dielectrics);
}

const Vector<Electrode>&
ComputationalGeometry::getElectrodes() const
{
  CH_TIME("ComputationalGeometry::getElectrodes()");
  if (m_verbose) {
    pout() << "ComputationalGeometry::getElectrodes()" << endl;
  }

  return (m_electrodes);
}

const RefCountedPtr<BaseIF>&
ComputationalGeometry::getGasImplicitFunction() const
{
  CH_TIME("ComputationalGeometry::getGasImplicitFunction()");
  if (m_verbose) {
    pout() << "ComputationalGeometry::getGasImplicitFunction()" << endl;
  }

  return (m_implicitFunctionGas);
}

const RefCountedPtr<BaseIF>&
ComputationalGeometry::getSolidImplicitFunction() const
{
  CH_TIME("ComputationalGeometry::getSolidImplicitFunction()");
  if (m_verbose) {
    pout() << "ComputationalGeometry::getSolidImplicitFunction()" << endl;
  }

  return (m_implicitFunctionSolid);
}

const RefCountedPtr<BaseIF>&
ComputationalGeometry::getImplicitFunction(const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::getImplicitFunction(phase::which_phase)");
  if (m_verbose) {
    pout() << "ComputationalGeometry::getImplicitFunction(phase::which_phase)" << endl;
  }

  return (a_phase == phase::gas) ? m_implicitFunctionGas : m_implicitFunctionSolid;
}

Real
ComputationalGeometry::getGasPermittivity() const
{
  CH_TIME("ComputationalGeometry::getGasPermittivity()");
  if (m_verbose) {
    pout() << "ComputationalGeometry::getGasPermittivity()" << endl;
  }

  return (m_eps0);
}

const RefCountedPtr<MultiFluidIndexSpace>&
ComputationalGeometry::getMfIndexSpace() const
{
  CH_TIME("ComputationalGeometry::getMfIndexSpace()");
  if (m_verbose) {
    pout() << "ComputationalGeometry::getMfIndexSpace()" << endl;
  }

  return (m_multifluidIndexSpace);
}

void
ComputationalGeometry::setDielectrics(const Vector<Dielectric>& a_dielectrics)
{
  CH_TIME("ComputationalGeometry::setDielectrics(Vector<Dielectric>)");
  if (m_verbose) {
    pout() << "ComputationalGeometry::setDielectrics(Vector<Dielectric>)" << endl;
  }

  m_dielectrics = a_dielectrics;
}

void
ComputationalGeometry::setElectrodes(const Vector<Electrode>& a_electrodes)
{
  CH_TIME("ComputationalGeometry::setElectrodes(Vector<Electrode>)");
  if (m_verbose) {
    pout() << "ComputationalGeometry::setElectrodes(Vector<Electrode>)" << endl;
  }

  m_electrodes = a_electrodes;
}

void
ComputationalGeometry::setGasPermittivity(const Real a_eps0)
{
  CH_TIME("ComputationalGeometry::setGasPermittivity(Real)");
  if (m_verbose) {
    pout() << "ComputationalGeometry::setGasPermittivity(Real)" << endl;
  }

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
  if (m_verbose) {
    pout() << "ComputationalGeometry::buildGeometries(ProblemDomain, RealVect, Real, int, int, int)" << endl;
  }

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
  if (m_verbose) {
    pout() << "ComputationalGeometry::buildImplicitFunctions()" << endl;
  }

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
                                 const int            a_maxGhostEB)
{
  CH_TIME("ComputationalGeometry::makeGrids");
  if (m_verbose) {
    pout() << "ComputationalGeometry::makeGrids" << endl;
  }

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
  if (m_minBlockSize <= 0) {
    MayDay::Error("ComputationalGeometry::makeGrids - ComputationalGeometry.min_block_size must be positive");
  }
  if (m_maxBlockSize <= 0) {
    MayDay::Error("ComputationalGeometry::makeGrids - ComputationalGeometry.max_block_size must be positive");
  }
  if (m_maxBlockSize % m_minBlockSize != 0) {
    MayDay::Error(
      "ComputationalGeometry::makeGrids - ComputationalGeometry.max_block_size must be a multiple of min_block_size");
  }

  // The one-tile nesting buffer between levels has to cover the ghost cells.
  if (m_minBlockSize < 2 * a_maxGhostEB) {
    MayDay::Error("ComputationalGeometry::makeGrids - the tile size must be at least twice the ghost width");
  }

  // The start domain and every level above it are tiled, so they must decompose into whole tiles. The levels
  // below the start domain are built whole and box by box, as ScanShop builds them, and may be smaller than a
  // tile; nothing is required of them.
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (a_startDomain.domainBox().size(dir) % m_minBlockSize != 0) {
      MayDay::Error("ComputationalGeometry::makeGrids - the start domain does not decompose into whole tiles");
    }
  }

  m_probLo      = a_probLo;
  m_maxGhostEB  = a_maxGhostEB;
  m_refineAngle = a_refineAngle;

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

  m_cutTiles.resize(numLevels);
  m_boxes.resize(numLevels);
  m_splitCounts.resize(numLevels, Vector<int>(7, 0));
  m_splitBoxes.resize(numLevels);
  m_splitReasons.resize(numLevels);
  m_gasTypes.resize(numLevels);
  m_solidTypes.resize(numLevels);

  this->buildImplicitFunctions();

  // Without any implicit function every level is one regular box.
  if (m_implicitFunctionGas.isNull() && m_implicitFunctionSolid.isNull()) {
    for (int lvl = 0; lvl < numLevels; lvl++) {
      m_boxes[lvl].push_back(m_domains[lvl].domainBox());
      m_gasTypes[lvl].push_back(GeometryService::Regular);
      m_solidTypes[lvl].push_back(GeometryService::Regular);
    }

    this->buildBoxTrees();

    return;
  }

  // The algorithm, in the order it runs (the design record is kept with the roadmap, issue #733). One box
  // hierarchy serves both phases; every box is classified by both implicit functions.
  //
  //   0. Start level: domainSplit the whole domain and classify every box in both phases (buildStartLevel,
  //      classifyBox).
  //   1. Upward, to the stop level: a box regular or covered in both phases refines whole with its tags. A box
  //      irregular in some phase is split into classified pieces if, in any such phase, the implicit function's
  //      normal turns by more than m_refineAngle between neighbouring cells near the surface; otherwise it is a
  //      leaf and nothing is built above it (buildFinerLevels, exceedsCurvature).
  //   2. Tiles, once: the boxes irregular in either phase on every level are tiled by TiledMeshRefine into a
  //      properly nested set (makeTiles). This is the coverage the simulation regrids onto.
  //   3. Every tile lies inside a box or above a leaf; per phase it inherits a regular or covered box's tag and
  //      is otherwise classified by that phase's implicit function at its own level (classifyTiles).
  //   4. Every hit box is cut down to what the tiles left of it (decimateBoxes).
  //   5. The levels coarser than the start level, whole and classified box by box as ScanShop builds them; then,
  //      from the finest level down, a box containing a box irregular in a phase is irregular in that phase
  //      (buildCoarserLevels).
  Vector<Vector<int>> firstChild(numLevels);
  Vector<Vector<int>> numChildren(numLevels);

  Timer timer("ComputationalGeometry::makeGrids");

  timer.startEvent("Start level");
  this->buildStartLevel();
  this->buildBoxTree(m_startLevel);
  timer.stopEvent("Start level");

  timer.startEvent("Upward pass");
  this->buildFinerLevels(firstChild, numChildren);
  timer.stopEvent("Upward pass");

  timer.startEvent("Tiles");
  this->makeTiles();
  this->buildTileTrees();
  timer.stopEvent("Tiles");

  Vector<Vector<GeometryService::InOut>> gasTileTypes(numLevels);
  Vector<Vector<GeometryService::InOut>> solidTileTypes(numLevels);
  Vector<Vector<int>>                    tileHosts(numLevels);

  timer.startEvent("Classify tiles");
  this->classifyTiles(firstChild, numChildren, gasTileTypes, solidTileTypes, tileHosts);
  timer.stopEvent("Classify tiles");

  timer.startEvent("Decimate boxes");
  this->decimateBoxes(gasTileTypes, solidTileTypes, tileHosts);
  timer.stopEvent("Decimate boxes");

  timer.startEvent("Coarser levels");
  this->buildCoarserLevels();
  timer.stopEvent("Coarser levels");

  // The start level is whole and not tiled, but its boxes irregular in either phase are its cut tiles all the
  // same: super-tiles by construction, and where the level's cut cells are. Taken after the push-down, which can
  // make a start-level box irregular for what the tiles above it hold.
  for (int i = 0; i < m_boxes[m_startLevel].size(); i++) {
    if (m_gasTypes[m_startLevel][i] == GeometryService::Irregular ||
        m_solidTypes[m_startLevel][i] == GeometryService::Irregular) {
      m_cutTiles[m_startLevel].push_back(m_boxes[m_startLevel][i]);
    }
  }

  timer.startEvent("Index the boxes");
  this->buildBoxTrees();
  timer.stopEvent("Index the boxes");

  if (m_profile) {
    this->reportGrids();

    timer.eventReport(pout(), false);
  }
}

int
ComputationalGeometry::getNumGridLevels() const noexcept
{
  return m_domains.size();
}

GeometryService::InOut
ComputationalGeometry::classify(const Box& a_box, const int a_level, const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::classify");
  if (m_verbose) {
    pout() << "ComputationalGeometry::classify" << endl;
  }

  if (a_level < 0 || a_level >= m_boxes.size()) {
    MayDay::Error("ComputationalGeometry::classify - no such level");
  }
  if (a_box.isEmpty() || !m_domains[a_level].domainBox().contains(a_box)) {
    MayDay::Error("ComputationalGeometry::classify - the box is empty or not inside the level's domain");
  }

  const Vector<Box>&                    boxes = m_boxes[a_level];
  const Vector<GeometryService::InOut>& types = this->types(a_phase)[a_level];

  // A level with no boxes is not described here at all: nothing refined that far because everything below it is
  // a leaf, and the level below holds what there is to say -- including that the surface is in there, which an
  // empty list on its own would not say. A level with boxes and no index is one makeGrids has not finished with.
  if (boxes.size() == 0) {
    if (a_level == 0) {
      MayDay::Error("ComputationalGeometry::classify - the coarsest level has no boxes");
    }

    return this->classify(coarsen(a_box, 2), a_level - 1, a_phase);
  }

  if (a_level >= m_boxTrees.size() || !m_boxTrees[a_level]) {
    MayDay::Error("ComputationalGeometry::classify - the level's boxes are not indexed");
  }

  bool anyRegular = false;
  bool anyCovered = false;
  bool decided    = false;

  // The query box as a bounding volume, its cells taken as the unit cubes they are. A node is entered if its
  // bounds meet it; the exact box test then filters the candidates, which a touching node can produce.
  const BV query = this->boundingVolume(a_box);

  using Node = BoxTree::Node;

  EBGeometry::BVH::NodeKeyFactory<Node, bool> nodeKey = [&query](const Node& a_node) noexcept -> bool {
    return a_node.m_bv.intersects(query);
  };

  EBGeometry::BVH::PrunePredicate<Node, bool> prune = [&decided](const Node& /*a_node*/,
                                                                 const bool& a_meets) noexcept -> bool {
    return a_meets && !decided;
  };

  EBGeometry::BVH::PackedChildOrderer<bool, K> orderer =
    [](std::array<std::pair<uint32_t, bool>, K>& /*a_children*/) noexcept -> void {
  };

  EBGeometry::BVH::PackedLeafEvaluator<int, EBGeometry::BVH::ValueStorage<int>> evaluate =
    [&](const std::vector<int>& a_indices, size_t a_offset, size_t a_count) noexcept -> void {
    for (size_t i = a_offset; i < a_offset + a_count && !decided; i++) {
      const int index = a_indices[i];

      if (!boxes[index].intersectsNotEmpty(a_box)) {
        continue;
      }

      switch (types[index]) {
      case GeometryService::Regular: {
        anyRegular = true;

        break;
      }
      case GeometryService::Covered: {
        anyCovered = true;

        break;
      }
      default: {
        decided = true;

        break;
      }
      }

      // a box that meets both a regular and a covered box has the surface between them
      decided = decided || (anyRegular && anyCovered);
    }
  };

  m_boxTrees[a_level]->traverse(evaluate, prune, orderer, nodeKey);

  if (decided) {
    return GeometryService::Irregular;
  }

  return anyCovered ? GeometryService::Covered : GeometryService::Regular;
}

GeometryService::InOut
ComputationalGeometry::classify(const Box& a_box, const ProblemDomain& a_domain, const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::classify(ProblemDomain)");
  if (m_verbose) {
    pout() << "ComputationalGeometry::classify(ProblemDomain)" << endl;
  }

  const int level = this->getLevel(a_domain);

  if (level >= 0) {
    return this->classify(a_box, level, a_phase);
  }

  // Finer than every stored level: the finest one answers for the box coarsened onto it. Its domain must be a
  // refinement by two of the stored domain, or the two index spaces do not line up.
  const int finest = m_domains.size() - 1;

  const Box& fineBox = a_domain.domainBox();
  const Box& coarBox = m_domains[finest].domainBox();

  int ratio = 1;

  while (ratio * coarBox.size(0) < fineBox.size(0)) {
    ratio *= 2;
  }

  if (refine(m_domains[finest], ratio) != a_domain) {
    MayDay::Error("ComputationalGeometry::classify - the domain is not a refinement by two of a stored one");
  }

  return this->classify(coarsen(a_box, ratio), finest, a_phase);
}

int
ComputationalGeometry::getStartLevel() const noexcept
{
  return m_startLevel;
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

Vector<Box>
ComputationalGeometry::getBoxes(const phase::which_phase     a_phase,
                                const int                    a_level,
                                const GeometryService::InOut a_type) const noexcept
{
  Vector<Box> boxes;

  const Vector<Box>&                    levelBoxes = m_boxes[a_level];
  const Vector<GeometryService::InOut>& levelTypes = this->types(a_phase)[a_level];

  for (int i = 0; i < levelBoxes.size(); i++) {
    if (levelTypes[i] == a_type) {
      boxes.push_back(levelBoxes[i]);
    }
  }

  return boxes;
}

const Vector<Box>&
ComputationalGeometry::getBoxes(const int a_level) const noexcept
{
  return m_boxes[a_level];
}

const Vector<Box>&
ComputationalGeometry::getCutTiles(const int a_level) const noexcept
{
  return m_cutTiles[a_level];
}

const Vector<Box>&
ComputationalGeometry::getSplitBoxes(const int a_level, Vector<int>& a_reasons) const noexcept
{
  a_reasons = m_splitReasons[a_level];

  return m_splitBoxes[a_level];
}

const Vector<GeometryService::InOut>&
ComputationalGeometry::getTypes(const phase::which_phase a_phase, const int a_level) const noexcept
{
  return this->types(a_phase)[a_level];
}

const ProblemDomain&
ComputationalGeometry::getDomain(const int a_level) const noexcept
{
  return m_domains[a_level];
}

Real
ComputationalGeometry::getDx(const int a_level) const noexcept
{
  return m_dx[a_level];
}

Vector<Vector<GeometryService::InOut>>&
ComputationalGeometry::types(const phase::which_phase a_phase) noexcept
{
  return (a_phase == phase::gas) ? m_gasTypes : m_solidTypes;
}

const Vector<Vector<GeometryService::InOut>>&
ComputationalGeometry::types(const phase::which_phase a_phase) const noexcept
{
  return (a_phase == phase::gas) ? m_gasTypes : m_solidTypes;
}

void
ComputationalGeometry::buildStartLevel()
{
  CH_TIME("ComputationalGeometry::buildStartLevel");
  if (m_verbose) {
    pout() << "ComputationalGeometry::buildStartLevel" << endl;
  }

  // The start domain decomposes into whole tiles (checked in makeGrids), so the block factor is the tile: every
  // box is a whole number of tiles, and at most a super-tile wide.
  domainSplit(m_domains[m_startLevel], m_boxes[m_startLevel], m_maxBlockSize, m_minBlockSize);

  this->classifyBoxes(m_boxes[m_startLevel], m_startLevel, m_gasTypes[m_startLevel], m_solidTypes[m_startLevel]);
}

void
ComputationalGeometry::buildFinerLevels(Vector<Vector<int>>& a_firstChild, Vector<Vector<int>>& a_numChildren)
{
  CH_TIME("ComputationalGeometry::buildFinerLevels");
  if (m_verbose) {
    pout() << "ComputationalGeometry::buildFinerLevels" << endl;
  }

  auto isIrregular = [](const GeometryService::InOut a_type) -> bool {
    return a_type == GeometryService::Irregular;
  };

  for (int lvl = m_startLevel; lvl < m_stopLevel; lvl++) {
    const Vector<Box>&                    boxes      = m_boxes[lvl];
    const Vector<GeometryService::InOut>& gasTypes   = m_gasTypes[lvl];
    const Vector<GeometryService::InOut>& solidTypes = m_solidTypes[lvl];

    a_firstChild[lvl].resize(boxes.size(), -1);
    a_numChildren[lvl].resize(boxes.size(), 0);

    // A box regular or covered in both phases refines whole, keeping both tags.
    for (int i = 0; i < boxes.size(); i++) {
      if (!isIrregular(gasTypes[i]) && !isIrregular(solidTypes[i])) {
        a_firstChild[lvl][i]  = m_boxes[lvl + 1].size();
        a_numChildren[lvl][i] = 1;

        m_boxes[lvl + 1].push_back(refine(boxes[i], 2));
        m_gasTypes[lvl + 1].push_back(gasTypes[i]);
        m_solidTypes[lvl + 1].push_back(solidTypes[i]);
      }
    }

    // A box irregular in some phase splits if the surface inside it turns too sharply for this level in any
    // such phase; otherwise it is a leaf and nothing is built above it. The pieces of every box that splits are
    // classified together so that the work is shared once per level rather than once per box.
    const Vector<int> flags = this->splitFlags(boxes, lvl, gasTypes, solidTypes);

    for (int i = 0; i < boxes.size(); i++) {
      m_splitCounts[lvl][flags[i]]++;

      if (flags[i] != 0) {
        m_splitBoxes[lvl].push_back(boxes[i]);
        m_splitReasons[lvl].push_back(flags[i]);
      }
    }

    Vector<Box> pieces;

    for (int i = 0; i < boxes.size(); i++) {
      if (flags[i] != 0) {
        Vector<Box> split;

        domainSplit(refine(boxes[i], 2), split, m_maxBlockSize, m_minBlockSize);

        a_firstChild[lvl][i]  = m_boxes[lvl + 1].size() + pieces.size();
        a_numChildren[lvl][i] = split.size();

        pieces.append(split);
      }
    }

    Vector<GeometryService::InOut> pieceGasTypes;
    Vector<GeometryService::InOut> pieceSolidTypes;

    this->classifyBoxes(pieces, lvl + 1, pieceGasTypes, pieceSolidTypes);

    m_boxes[lvl + 1].append(pieces);
    m_gasTypes[lvl + 1].append(pieceGasTypes);
    m_solidTypes[lvl + 1].append(pieceSolidTypes);
  }
}

void
ComputationalGeometry::classifyBoxes(const Vector<Box>&              a_boxes,
                                     const int                       a_level,
                                     Vector<GeometryService::InOut>& a_gasTypes,
                                     Vector<GeometryService::InOut>& a_solidTypes) const
{
  CH_TIME("ComputationalGeometry::classifyBoxes");
  if (m_verbose) {
    pout() << "ComputationalGeometry::classifyBoxes" << endl;
  }

  // The classifications travel as integers so that one all-reduce assembles them: a rank writes only the
  // entries it owns and leaves the rest at zero, and the sum is the union. Two entries per box, gas then solid.
  constexpr int regular   = 1;
  constexpr int covered   = 2;
  constexpr int irregular = 3;

  auto encode = [&](const GeometryService::InOut a_type) -> int {
    switch (a_type) {
    case GeometryService::Regular: {
      return regular;
    }
    case GeometryService::Covered: {
      return covered;
    }
    default: {
      return irregular;
    }
    }
  };

  auto decode = [&](const int a_code) -> GeometryService::InOut {
    switch (a_code) {
    case regular: {
      return GeometryService::Regular;
    }
    case covered: {
      return GeometryService::Covered;
    }
    case irregular: {
      return GeometryService::Irregular;
    }
    default: {
      MayDay::Error("ComputationalGeometry::classifyBoxes - a box was classified by no rank or by several");

      return GeometryService::Irregular;
    }
    }
  };

  Vector<int> codes(2 * a_boxes.size(), 0);

  for (int i = procID(); i < a_boxes.size(); i += numProc()) {
    codes[2 * i]     = encode(this->classifyBox(a_boxes[i], a_level, phase::gas));
    codes[2 * i + 1] = encode(this->classifyBox(a_boxes[i], a_level, phase::solid));
  }

  ParallelOps::sum(codes);

  a_gasTypes.resize(a_boxes.size());
  a_solidTypes.resize(a_boxes.size());

  for (int i = 0; i < a_boxes.size(); i++) {
    a_gasTypes[i]   = decode(codes[2 * i]);
    a_solidTypes[i] = decode(codes[2 * i + 1]);
  }
}

Vector<int>
ComputationalGeometry::splitFlags(const Vector<Box>&                    a_boxes,
                                  const int                             a_level,
                                  const Vector<GeometryService::InOut>& a_gasTypes,
                                  const Vector<GeometryService::InOut>& a_solidTypes) const
{
  CH_TIME("ComputationalGeometry::splitFlags");
  if (m_verbose) {
    pout() << "ComputationalGeometry::splitFlags" << endl;
  }

  // The flag carries the reason, so the report can say why a level refined where it did.
  Vector<int> flags(a_boxes.size(), 0);

  for (int i = procID(); i < a_boxes.size(); i += numProc()) {
    SplitReason reason = SplitReason::None;

    if (a_gasTypes[i] == GeometryService::Irregular) {
      reason = this->exceedsCurvature(a_boxes[i], a_level, phase::gas);
    }

    if (reason == SplitReason::None && a_solidTypes[i] == GeometryService::Irregular) {
      reason = this->exceedsCurvature(a_boxes[i], a_level, phase::solid);
    }

    // A box the curvature leaves alone is a leaf, and a leaf is what the level above describes the other side
    // of. One that holds an edge crossed twice at the finer spacing cannot describe its own face there, so it
    // is refined until the crossing is resolved.
    if (reason == SplitReason::None) {
      const bool doubled = (a_gasTypes[i] == GeometryService::Irregular &&
                            this->doublyCrossedEdge(a_boxes[i], a_level, phase::gas)) ||
                           (a_solidTypes[i] == GeometryService::Irregular &&
                            this->doublyCrossedEdge(a_boxes[i], a_level, phase::solid));

      if (doubled) {
        reason = SplitReason::DoubleCrossing;
      }
    }

    flags[i] = static_cast<int>(reason);
  }

  ParallelOps::sum(flags);

  return flags;
}

GeometryService::InOut
ComputationalGeometry::classifyBox(const Box& a_box, const int a_level, const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::classifyBox");
  if (m_verbose) {
    pout() << "ComputationalGeometry::classifyBox" << endl;
  }

  const RefCountedPtr<BaseIF>& implicitFunction = this->getImplicitFunction(a_phase);

  if (implicitFunction.isNull()) {
    return GeometryService::Regular;
  }

  // The box is irregular if some cell of it, grown by the ghost width, is cut as the polyhedral shop reads cut:
  // its corner values disagree under isFluid. Two nodes of the grown region that disagree are joined by a path
  // of adjacent nodes along which some adjacent pair disagrees, and adjacent nodes are corners of one cell of
  // the region, so one node of each kind is the test. Node values are what the shop reconstructs from, so this
  // is exact for what it builds, whatever the function does away from its zero set.
  const BaseIF& f = *implicitFunction;

  const Real dx    = m_dx[a_level];
  const Box  grown = grow(a_box, m_maxGhostEB) & m_domains[a_level].domainBox();

  Box nodeBox = grown;
  nodeBox.surroundingNodes();

  bool anyFluid = false;
  bool anySolid = false;

  // How close the surface comes to the nodes, for the test below, which is reached only where every node agrees.
  Real closest = std::numeric_limits<Real>::max();

  for (BoxIterator bit(nodeBox); bit.ok(); ++bit) {
    const IntVect iv = bit();

    RealVect x = m_probLo;

    for (int dir = 0; dir < SpaceDim; dir++) {
      x[dir] += dx * static_cast<Real>(iv[dir]);
    }

    const Real value = PolyhedralGeometryShop::snappedValue(f, x, dx);

    closest = std::min(closest, std::abs(value));

    if (PolyhedralEB::isFluid(value)) {
      anyFluid = true;
    }
    else {
      anySolid = true;
    }

    if (anyFluid && anySolid) {
      return GeometryService::Irregular;
    }
  }

  // Every node agrees, and the cells are regular or covered as far as their corners can tell. The surface can
  // still pass through the box between the nodes, entering and leaving through one edge, and a box like that is
  // not regular: the level above it would see the crossings this level cannot, and a cell of it left on the
  // coarse side of a level boundary could not describe its own face. Such a box is irregular, so that it is
  // tagged, tiled and refined like any other, and the level above resolves what this one cannot represent.
  //
  // The finest level is left alone: nothing finer describes it, so the answer would cost tiles and buy nothing.
  // Away from the surface nothing is evaluated at all: a midpoint can only disagree with two agreeing ends if
  // the surface comes within half a cell of the edge, so a box whose nearest node value exceeds a cell width is
  // done -- the same reading of the implicit function as a distance that the scan-based pruning already makes.
  if (a_level >= m_stopLevel || closest > dx) {
    return anySolid ? GeometryService::Covered : GeometryService::Regular;
  }

  const Real fineDx = 0.5 * dx;

  for (BoxIterator bit(nodeBox); bit.ok(); ++bit) {
    const IntVect iv = bit();

    for (int dir = 0; dir < SpaceDim; dir++) {
      const IntVect jv = iv + BASISV(dir);

      if (!nodeBox.contains(jv)) {
        continue;
      }

      RealVect x = m_probLo;

      for (int d = 0; d < SpaceDim; d++) {
        x[d] += dx * static_cast<Real>(iv[d]);
      }

      x[dir] += fineDx;

      // the midpoint is a node of the level above, and is read at that level's spacing
      if (PolyhedralEB::isFluid(PolyhedralGeometryShop::snappedValue(f, x, fineDx)) != anyFluid) {
        return GeometryService::Irregular;
      }
    }
  }

  return anySolid ? GeometryService::Covered : GeometryService::Regular;
}

bool
ComputationalGeometry::doublyCrossedEdge(const Box& a_box, const int a_level, const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::doublyCrossedEdge");
  if (m_verbose) {
    pout() << "ComputationalGeometry::doublyCrossedEdge" << endl;
  }

  const RefCountedPtr<BaseIF>& implicitFunction = this->getImplicitFunction(a_phase);

  if (implicitFunction.isNull()) {
    return false;
  }

  const BaseIF& f = *implicitFunction;

  const Real dx     = m_dx[a_level];
  const Real fineDx = 0.5 * dx;
  const Box  valid  = a_box & m_domains[a_level].domainBox();

  // The node values of this level, as the cells of the box are built from, and the midpoint of every edge
  // between two of them, which is a node of the level above. The midpoint is snapped at the finer spacing,
  // since that is the spacing the level above would classify it at.
  BaseFab<Real> nodeValues;

  PolyhedralGeometryShop::fillNodeValues(f, nodeValues, valid, m_probLo, dx);

  const Box& nodeBox = nodeValues.box();

  for (BoxIterator bit(nodeBox); bit.ok(); ++bit) {
    const IntVect iv = bit();

    for (int dir = 0; dir < SpaceDim; dir++) {
      const IntVect jv = iv + BASISV(dir);

      if (!nodeBox.contains(jv)) {
        continue;
      }

      const Real loValue = nodeValues(iv, 0);
      const Real hiValue = nodeValues(jv, 0);

      const bool loZero = (loValue == 0.0);
      const bool hiZero = (hiValue == 0.0);

      const bool loFluid = PolyhedralEB::isFluid(loValue);
      const bool hiFluid = PolyhedralEB::isFluid(hiValue);

      // An edge whose ends disagree carries one crossing at this level and one at the next, in the half its own
      // crossing lies in. Only ends that agree can hide a pair -- but an end that is exactly zero is on the
      // surface, and which side of it that end belongs to is not decided at this spacing. The fluid rule breaks
      // the tie toward solid, and an edge that reads fluid-solid-zero would read fluid-solid-fluid had the tie
      // gone the other way, which is a pair. Such an end is therefore left open, and the ends count as agreeing
      // if any reading of it makes them.
      if (!loZero && !hiZero && loFluid != hiFluid) {
        continue;
      }

      RealVect x = m_probLo;

      for (int d = 0; d < SpaceDim; d++) {
        x[d] += dx * static_cast<Real>(iv[d]);
      }

      x[dir] += 0.5 * dx;

      const Real midValue = PolyhedralGeometryShop::snappedValue(f, x, fineDx);

      if (midValue == 0.0) {
        continue;
      }

      const bool midFluid = PolyhedralEB::isFluid(midValue);

      // The ends are read as the midpoint's opposite wherever they are free to be, since that is the reading
      // that hides a pair.
      const bool loAgainst = loZero ? !midFluid : loFluid;
      const bool hiAgainst = hiZero ? !midFluid : hiFluid;

      if (loAgainst == hiAgainst && midFluid != loAgainst) {
        return true;
      }
    }
  }

  return false;
}

ComputationalGeometry::SplitReason
ComputationalGeometry::exceedsCurvature(const Box& a_box, const int a_level, const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::exceedsCurvature");
  if (m_verbose) {
    pout() << "ComputationalGeometry::exceedsCurvature" << endl;
  }

  const RefCountedPtr<BaseIF>& implicitFunction = this->getImplicitFunction(a_phase);

  if (implicitFunction.isNull()) {
    return SplitReason::None;
  }

  // Normals are those of the interface polygons the polyhedral shop builds on this level, in the cut cells of the
  // box and one cell beyond it so that a pair across the box boundary is seen from both sides. They come from the
  // edge roots alone, so they depend only on the zero set: the implicit function is not a distance function
  // inside a body built by CSG, and the gradient there carries the kinks and medial shells of the construction,
  // which are not features of the surface. The node values are evaluated once per node and each crossed edge is
  // bisected once, shared by the cells around it, as the shop does when it builds the graph.
  const BaseIF& f = *implicitFunction;

  const Real dx    = m_dx[a_level];
  const Box  valid = a_box & m_domains[a_level].domainBox();
  const Box  grown = grow(a_box, 1) & m_domains[a_level].domainBox();

  BaseFab<Real> nodeValues;
  BaseFab<Real> intercept[SpaceDim];

  PolyhedralGeometryShop::fillNodeValues(f, nodeValues, grown, m_probLo, dx);
  PolyhedralGeometryShop::defineIntercepts(intercept, grown);

  BaseFab<Real> normal(grown, SpaceDim);
  BaseFab<int>  isCut(grown, 1);

  isCut.setVal(0);

  for (BoxIterator bit(grown); bit.ok(); ++bit) {
    const IntVect iv = bit();

    PolyhedralEB::CutCellSurface surface;

    PolyhedralGeometryShop::buildSurface(f, intercept, surface, nodeValues, iv, m_probLo, dx);

    if (PolyhedralEB::CutCellBody::classify(surface) != PolyhedralEB::CutCellBody::Kind::Cut) {
      continue;
    }

    // A body that does not close has no normal to compare; the shop reports such a cell itself when it builds
    // the level.
    PolyhedralEB::CutCellBody body;

    if (!body.define(surface)) {
      continue;
    }

    const RealVect n = body.normal();

    if (n.vectorLength() > 0.0) {
      isCut(iv, 0) = 1;

      for (int dir = 0; dir < SpaceDim; dir++) {
        normal(iv, dir) = n[dir];
      }
    }
  }

  const Real cosThreshold = std::cos(m_refineAngle * Units::pi / 180.0);

  SplitReason ring = SplitReason::None;

  for (BoxIterator bit(valid); bit.ok(); ++bit) {
    const IntVect iv = bit();

    if (isCut(iv, 0) == 0) {
      continue;
    }

    const Box neighbours = grow(Box(iv, iv), 1) & grown;

    for (BoxIterator nit(neighbours); nit.ok(); ++nit) {
      const IntVect jv = nit();

      if (jv == iv || isCut(jv, 0) == 0) {
        continue;
      }

      Real dot = 0.0;

      for (int dir = 0; dir < SpaceDim; dir++) {
        dot += normal(iv, dir) * normal(jv, dir);
      }

      // facing normals split whatever the threshold, so they are tested ahead of it
      if (dot < 0.0 || dot < cosThreshold) {
        // a pair inside the box is the box's own doing and decides at once; a pair reaching into the ring
        // is recorded and decides only if the whole box turns out to have no pair of its own, so that the
        // report attributes the split to the box itself whenever it can
        if (valid.contains(jv)) {
          return (dot < 0.0) ? SplitReason::Medial : SplitReason::Interior;
        }

        ring = (dot < 0.0) ? SplitReason::Medial : ((ring == SplitReason::None) ? SplitReason::Ring : ring);
      }
    }
  }

  // No pair of its own: a pair reaching into the ring decides.
  return ring;
}

void
ComputationalGeometry::makeTiles()
{
  CH_TIME("ComputationalGeometry::makeTiles");
  if (m_verbose) {
    pout() << "ComputationalGeometry::makeTiles" << endl;
  }

  const int numAbove = m_stopLevel - m_startLevel;

  if (numAbove == 0) {
    return;
  }

  // TiledMeshRefine tiles level k from tags on level k - 1 and takes the start domain as its level 0, so a box
  // irregular in either phase on builder level lvl enters, coarsened by two, as a tag for tiler level
  // lvl - m_startLevel. Tags are rank-local and the tiler gathers them, so each rank tags only its share.
  Vector<IntVectSet> tags(numAbove);

  for (int lvl = m_startLevel + 1; lvl <= m_stopLevel; lvl++) {
    IntVectSet& levelTags = tags[lvl - m_startLevel - 1];

    for (int i = procID(); i < m_boxes[lvl].size(); i += numProc()) {
      if (m_gasTypes[lvl][i] == GeometryService::Irregular || m_solidTypes[lvl][i] == GeometryService::Irregular) {
        levelTags |= coarsen(m_boxes[lvl][i], 2);
      }
    }
  }

  const Vector<int> refRatios(1 + numAbove, 2);

  TiledMeshRefine tiler(m_domains[m_startLevel],
                        refRatios,
                        m_minBlockSize * IntVect::Unit,
                        m_maxBlockSize * IntVect::Unit);

  // The tiles are built, then read back: a cell left on the coarse side of a level boundary whose edge the level
  // above crosses twice is tagged, and the tiles are built again. The pass has to come after the tiles exist
  // rather than before, because the tiled region is not a function of the box classification alone -- the tiler
  // nests, and nesting puts tiles above boxes that never refined, which is exactly where such a cell hides. Each
  // pass only adds tags, so the region grows and the loop ends; in practice one further pass finds nothing.
  for (int pass = 0; pass < s_maxTilePasses; pass++) {
    Vector<Vector<Box>> tiles;

    const int finestTiled = tiler.regrid(tiles, tags);

    for (int lvl = m_startLevel + 1; lvl <= m_stopLevel; lvl++) {
      m_cutTiles[lvl].clear();
    }

    // Tiler level 0 is the start domain, whole and not tiled, and is discarded as AmrMesh discards it.
    for (int lvl = m_startLevel + 1; lvl <= m_startLevel + finestTiled; lvl++) {
      m_cutTiles[lvl] = tiles[lvl - m_startLevel];
    }

    this->buildTileTrees();

    if (!this->tagUnresolvedSeams(tags)) {
      return;
    }

    if (m_verbose && procID() == 0) {
      pout() << "ComputationalGeometry::makeTiles - tiling again for the cells the level above would cross twice"
             << endl;
    }
  }

  MayDay::Error("ComputationalGeometry::makeTiles - the tiles never resolved every doubly crossed edge");
}

bool
ComputationalGeometry::tagUnresolvedSeams(Vector<IntVectSet>& a_tags) const
{
  CH_TIME("ComputationalGeometry::tagUnresolvedSeams");
  if (m_verbose) {
    pout() << "ComputationalGeometry::tagUnresolvedSeams" << endl;
  }

  const phase::which_phase phases[2] = {phase::gas, phase::solid};

  long long numTagged = 0;

  for (int lvl = m_startLevel; lvl < m_stopLevel; lvl++) {
    // On the start level every box is a candidate, since the level is whole; above it the tiles are what the
    // graph holds and the rest of the level is not described there at all.
    const Vector<Box>& coarse = (lvl == m_startLevel) ? m_boxes[lvl] : m_cutTiles[lvl];

    const Box& domainBox = m_domains[lvl].domainBox();

    IntVectSet& levelTags = a_tags[lvl - m_startLevel];

    for (int i = procID(); i < coarse.size(); i += numProc()) {
      const Box box   = coarse[i];
      const Box grown = grow(box, 1) & domainBox;

      // What the level above covers here, read off its tiles through their index.
      BaseFab<bool> refined(grown, 1);

      refined.setVal(false);

      const Vector<int> hits = this->tilesMeeting(lvl + 1, refine(grown, 2));

      for (int h = 0; h < hits.size(); h++) {
        const Box covered = coarsen(m_cutTiles[lvl + 1][hits[h]], 2) & grown;

        if (!covered.isEmpty()) {
          refined.setVal(true, covered, 0);
        }
      }

      for (BoxIterator bit(box); bit.ok(); ++bit) {
        const IntVect iv = bit();

        if (refined(iv, 0)) {
          continue;
        }

        // Only the coarse side of a level boundary: a cell whose neighbours are all at its own level describes
        // its faces with the same nodes they do. Face neighbours are the whole test, even though a doubly
        // crossed edge is shared by the four cells meeting along it and one of those can meet the refined
        // region along that edge alone. Of the other three, the two that share a face with the refined one are
        // tested here, carry the same edge, and are therefore tagged; once they are refined the fourth has a
        // refined face neighbour and is taken on the next pass, which is what the passes are for.
        bool onSeam = false;

        for (int dir = 0; dir < SpaceDim && !onSeam; dir++) {
          for (int side = 0; side < 2 && !onSeam; side++) {
            const IntVect jv = iv + (2 * side - 1) * BASISV(dir);

            onSeam = grown.contains(jv) && refined(jv, 0);
          }
        }

        if (!onSeam) {
          continue;
        }

        for (const phase::which_phase& curPhase : phases) {
          if (this->doublyCrossedEdge(Box(iv, iv), lvl, curPhase)) {
            levelTags |= iv;

            numTagged++;

            break;
          }
        }
      }
    }
  }

  return ParallelOps::sum(numTagged) > 0;
}

ComputationalGeometry::BV
ComputationalGeometry::boundingVolume(const Box& a_box) noexcept
{
  // A cell is the unit cube whose corners are its nodes, so a box spans [smallEnd, bigEnd + 1] in index space.
  // The bounding volumes are three-dimensional whatever SpaceDim is, and two of them overlap only where they do
  // so on every axis, so in two dimensions the third axis is given the unit thickness a cell has there too.
  Vec3 lo = Vec3(0.0, 0.0, 0.0);
  Vec3 hi = Vec3(1.0, 1.0, 1.0);

  for (int dir = 0; dir < SpaceDim; dir++) {
    lo[dir] = static_cast<Real>(a_box.smallEnd(dir));
    hi[dir] = static_cast<Real>(a_box.bigEnd(dir) + 1);
  }

  return BV(lo, hi);
}

void
ComputationalGeometry::buildBoxTree(const int a_level)
{
  CH_TIME("ComputationalGeometry::buildBoxTree");
  if (m_verbose) {
    pout() << "ComputationalGeometry::buildBoxTree" << endl;
  }

  m_boxTrees.resize(m_boxes.size());

  m_boxTrees[a_level] = this->buildTree(m_boxes[a_level]);
}

void
ComputationalGeometry::buildBoxTrees()
{
  CH_TIME("ComputationalGeometry::buildBoxTrees");
  if (m_verbose) {
    pout() << "ComputationalGeometry::buildBoxTrees" << endl;
  }

  m_boxTrees.resize(m_boxes.size());

  for (int lvl = 0; lvl < m_boxes.size(); lvl++) {
    m_boxTrees[lvl] = this->buildTree(m_boxes[lvl]);
  }
}

void
ComputationalGeometry::buildTileTrees()
{
  CH_TIME("ComputationalGeometry::buildTileTrees");
  if (m_verbose) {
    pout() << "ComputationalGeometry::buildTileTrees" << endl;
  }

  m_tileTrees.resize(m_cutTiles.size());

  for (int lvl = 0; lvl < m_cutTiles.size(); lvl++) {
    m_tileTrees[lvl] = this->buildTree(m_cutTiles[lvl]);
  }
}

std::shared_ptr<ComputationalGeometry::BoxTree>
ComputationalGeometry::buildTree(const Vector<Box>& a_boxes) const
{
  if (a_boxes.size() == 0) {
    return nullptr;
  }

  std::vector<std::pair<int, BV>> primitives;

  primitives.reserve(a_boxes.size());

  for (int i = 0; i < a_boxes.size(); i++) {
    primitives.emplace_back(i, ComputationalGeometry::boundingVolume(a_boxes[i]));
  }

  return std::make_shared<BoxTree>(std::move(primitives), s_treeLeafSize);
}

Vector<int>
ComputationalGeometry::meeting(const std::shared_ptr<BoxTree>& a_tree, const Vector<Box>& a_boxes, const Box& a_box)
{
  Vector<int> found;

  if (!a_tree) {
    return found;
  }

  using Node = BoxTree::Node;

  const BV query = ComputationalGeometry::boundingVolume(a_box);

  EBGeometry::BVH::NodeKeyFactory<Node, bool> nodeKey = [&query](const Node& a_node) noexcept -> bool {
    return a_node.m_bv.intersects(query);
  };

  EBGeometry::BVH::PrunePredicate<Node, bool> prune = [](const Node& /*a_node*/, const bool& a_meets) noexcept -> bool {
    return a_meets;
  };

  EBGeometry::BVH::PackedChildOrderer<bool, K> orderer =
    [](std::array<std::pair<uint32_t, bool>, K>& /*a_children*/) noexcept -> void {
  };

  EBGeometry::BVH::PackedLeafEvaluator<int, EBGeometry::BVH::ValueStorage<int>> evaluate =
    [&](const std::vector<int>& a_indices, size_t a_offset, size_t a_count) noexcept -> void {
    for (size_t i = a_offset; i < a_offset + a_count; i++) {
      if (a_boxes[a_indices[i]].intersectsNotEmpty(a_box)) {
        found.push_back(a_indices[i]);
      }
    }
  };

  a_tree->traverse(evaluate, prune, orderer, nodeKey);

  // The tree pads its K-ary nodes by repeating a child, so a leaf can be entered more than once and a box be
  // found more than once; the caller is answered with each box once, in increasing index order.
  std::vector<int>& hits = found.stdVector();

  std::sort(hits.begin(), hits.end());

  hits.erase(std::unique(hits.begin(), hits.end()), hits.end());

  return found;
}

Vector<int>
ComputationalGeometry::boxesMeeting(const int a_level, const Box& a_box) const
{
  // A level with no boxes has no tree and nothing to meet; a level with boxes and no tree was not indexed.
  if (a_level < 0 || a_level >= m_boxTrees.size() || (!m_boxTrees[a_level] && m_boxes[a_level].size() > 0)) {
    MayDay::Error("ComputationalGeometry::boxesMeeting - the level's boxes are not indexed");
  }

  return ComputationalGeometry::meeting(m_boxTrees[a_level], m_boxes[a_level], a_box);
}

Vector<int>
ComputationalGeometry::tilesMeeting(const int a_level, const Box& a_box) const
{
  // As boxesMeeting: no tiles is an empty answer, tiles without a tree is a level that was not indexed.
  if (a_level < 0 || a_level >= m_tileTrees.size() || (!m_tileTrees[a_level] && m_cutTiles[a_level].size() > 0)) {
    MayDay::Error("ComputationalGeometry::tilesMeeting - the level's tiles are not indexed");
  }

  return ComputationalGeometry::meeting(m_tileTrees[a_level], m_cutTiles[a_level], a_box);
}

int
ComputationalGeometry::containingBox(const int a_level, const IntVect& a_cell) const
{
  const Vector<int> hits = this->boxesMeeting(a_level, Box(a_cell, a_cell));

  // A whole level is a partition of its domain, so exactly one box holds the cell.
  if (hits.size() != 1) {
    MayDay::Error("ComputationalGeometry::containingBox - a cell of a whole level lies in no box or in several");
  }

  return hits[0];
}

void
ComputationalGeometry::classifyTiles(const Vector<Vector<int>>&              a_firstChild,
                                     const Vector<Vector<int>>&              a_numChildren,
                                     Vector<Vector<GeometryService::InOut>>& a_gasTileTypes,
                                     Vector<Vector<GeometryService::InOut>>& a_solidTileTypes,
                                     Vector<Vector<int>>&                    a_tileHosts) const
{
  CH_TIME("ComputationalGeometry::classifyTiles");
  if (m_verbose) {
    pout() << "ComputationalGeometry::classifyTiles" << endl;
  }

  constexpr int regular   = 1;
  constexpr int covered   = 2;
  constexpr int irregular = 3;

  auto encode = [&](const GeometryService::InOut a_type) -> int {
    return (a_type == GeometryService::Regular) ? regular : (a_type == GeometryService::Covered) ? covered : irregular;
  };

  auto decode = [&](const int a_code) -> GeometryService::InOut {
    if (a_code == regular) {
      return GeometryService::Regular;
    }
    else if (a_code == covered) {
      return GeometryService::Covered;
    }
    else if (a_code == irregular) {
      return GeometryService::Irregular;
    }

    MayDay::Error("ComputationalGeometry::classifyTiles - a tile was classified by no rank or by several");

    return GeometryService::Irregular;
  };

  for (int lvl = m_startLevel + 1; lvl <= m_stopLevel; lvl++) {
    const Vector<Box>& tiles = m_cutTiles[lvl];

    a_tileHosts[lvl].resize(tiles.size(), -1);

    // The box a tile lies in: the start-level box under it is the one holding its coarsening's first cell, and
    // from there the child links descend one level at a time, picking the child that contains the tile's
    // coarsening. Running out of children before the tile's level is the hole.
    for (int t = 0; t < tiles.size(); t++) {
      const IntVect startCell = coarsen(tiles[t], 1 << (lvl - m_startLevel)).smallEnd();

      int host = this->containingBox(m_startLevel, startCell);

      for (int k = m_startLevel; k < lvl && host >= 0; k++) {
        const int first = a_firstChild[k][host];
        const int num   = a_numChildren[k][host];

        const Box target = coarsen(tiles[t], 1 << (lvl - k - 1));

        host = -1;

        for (int c = first; c < first + num; c++) {
          if (m_boxes[k + 1][c].contains(target)) {
            host = c;

            break;
          }
        }
      }

      a_tileHosts[lvl][t] = host;
    }

    // Per phase, a tile inherits from a regular or covered host and is otherwise classified by the implicit
    // function at its own level: inside an irregular host because the host is conservative, in a hole because a
    // carried tile is generated at its level. Shared between the ranks by tile and assembled with one all-reduce.
    Vector<int> codes(2 * tiles.size(), 0);

    for (int t = procID(); t < tiles.size(); t += numProc()) {
      const int host = a_tileHosts[lvl][t];

      const phase::which_phase phases[2] = {phase::gas, phase::solid};

      for (int p = 0; p < 2; p++) {
        const GeometryService::InOut hostType = (host >= 0) ? this->types(phases[p])[lvl][host]
                                                            : GeometryService::Irregular;

        codes[2 * t + p] = (hostType == GeometryService::Irregular)
                             ? encode(this->classifyBox(tiles[t], lvl, phases[p]))
                             : encode(hostType);
      }
    }

    ParallelOps::sum(codes);

    a_gasTileTypes[lvl].resize(tiles.size());
    a_solidTileTypes[lvl].resize(tiles.size());

    for (int t = 0; t < tiles.size(); t++) {
      a_gasTileTypes[lvl][t]   = decode(codes[2 * t]);
      a_solidTileTypes[lvl][t] = decode(codes[2 * t + 1]);
    }
  }
}

void
ComputationalGeometry::decimateBoxes(const Vector<Vector<GeometryService::InOut>>& a_gasTileTypes,
                                     const Vector<Vector<GeometryService::InOut>>& a_solidTileTypes,
                                     const Vector<Vector<int>>&                    a_tileHosts)
{
  CH_TIME("ComputationalGeometry::decimateBoxes");
  if (m_verbose) {
    pout() << "ComputationalGeometry::decimateBoxes" << endl;
  }

  for (int lvl = m_startLevel + 1; lvl <= m_stopLevel; lvl++) {
    const Vector<Box>&                    oldBoxes      = m_boxes[lvl];
    const Vector<GeometryService::InOut>& oldGasTypes   = m_gasTypes[lvl];
    const Vector<GeometryService::InOut>& oldSolidTypes = m_solidTypes[lvl];
    const Vector<Box>&                    tiles         = m_cutTiles[lvl];

    // Which tiles each box hosts, as two flat arrays: hostedStart[i] .. hostedStart[i + 1] index into
    // hostedTiles for box i, filled by a counting sort over the tiles so that no box owns an allocation. For
    // the report, how many tiles are there because a box was irregular and how many because nesting put them
    // there.
    Vector<int> hostedStart(oldBoxes.size() + 1, 0);
    Vector<int> hostedTiles(tiles.size(), -1);

    for (int t = 0; t < tiles.size(); t++) {
      const int host = a_tileHosts[lvl][t];

      if (host >= 0) {
        hostedStart[host + 1]++;
      }

      const bool tagged = (host >= 0) && (oldGasTypes[host] == GeometryService::Irregular ||
                                          oldSolidTypes[host] == GeometryService::Irregular);

      m_splitCounts[lvl][tagged ? 5 : 6]++;
    }

    for (int i = 0; i < oldBoxes.size(); i++) {
      hostedStart[i + 1] += hostedStart[i];
    }

    Vector<int> hostedNext = hostedStart;

    for (int t = 0; t < tiles.size(); t++) {
      const int host = a_tileHosts[lvl][t];

      if (host >= 0) {
        hostedTiles[hostedNext[host]++] = t;
      }
    }

    // The tiles come first, with their own classifications.
    Vector<Box>                    newBoxes      = tiles;
    Vector<GeometryService::InOut> newGasTypes   = a_gasTileTypes[lvl];
    Vector<GeometryService::InOut> newSolidTypes = a_solidTileTypes[lvl];

    for (int i = 0; i < oldBoxes.size(); i++) {
      const bool irregular = (oldGasTypes[i] == GeometryService::Irregular) ||
                             (oldSolidTypes[i] == GeometryService::Irregular);

      const int firstHosted = hostedStart[i];
      const int numHosted   = hostedStart[i + 1] - hostedStart[i];

      if (numHosted == 0) {
        // Untouched: kept as it is.
        newBoxes.push_back(oldBoxes[i]);
        newGasTypes.push_back(oldGasTypes[i]);
        newSolidTypes.push_back(oldSolidTypes[i]);
      }
      else if (irregular) {
        // A box irregular in some phase is a super-tile the tiles cover in full; the tiles replace it.
        continue;
      }
      else {
        // A regular or covered box loses the tiles inside it. Done in tile coordinates, where the box and every
        // tile are whole cells, with TreeIntVectSet: one node for the box, a descent per tile, and one box per
        // remaining full node. Correct and octree-graded; not tight, and replaceable here.
        TreeIntVectSet remainder(coarsen(oldBoxes[i], m_minBlockSize));

        for (int j = firstHosted; j < firstHosted + numHosted; j++) {
          remainder -= coarsen(tiles[hostedTiles[j]], m_minBlockSize);
        }

        const Vector<Box> pieces = remainder.createBoxes();

        for (int j = 0; j < pieces.size(); j++) {
          newBoxes.push_back(refine(pieces[j], m_minBlockSize));
          newGasTypes.push_back(oldGasTypes[i]);
          newSolidTypes.push_back(oldSolidTypes[i]);
        }
      }
    }

    m_boxes[lvl]      = newBoxes;
    m_gasTypes[lvl]   = newGasTypes;
    m_solidTypes[lvl] = newSolidTypes;

    // The level's index held the boxes that were just replaced, and its primitives are indices into the list
    // that is gone. It is dropped rather than rebuilt here: makeGrids indexes every level once the lists are
    // final, and until then a query on this level meets the null check rather than the old list.
    if (lvl < m_boxTrees.size()) {
      m_boxTrees[lvl] = nullptr;
    }
  }
}

void
ComputationalGeometry::buildCoarserLevels()
{
  CH_TIME("ComputationalGeometry::buildCoarserLevels");
  if (m_verbose) {
    pout() << "ComputationalGeometry::buildCoarserLevels" << endl;
  }

  // The levels coarser than the start domain are built as ScanShop builds them: whole, each box classified on
  // its own, in both phases. No block factor, since these levels are not tiled and the coarsest may be smaller
  // than a tile.
  for (int lvl = m_startLevel - 1; lvl >= 0; lvl--) {
    domainSplit(m_domains[lvl], m_boxes[lvl], m_maxBlockSize);

    this->classifyBoxes(m_boxes[lvl], lvl, m_gasTypes[lvl], m_solidTypes[lvl]);

    this->buildBoxTree(lvl);
  }

  // The push-down: a box containing a box that is irregular in a phase is irregular in that phase, from the
  // finest level down. Above the start level the container of an irregular tile's coarsening is among the tiles
  // of the level below, by nesting; at and below the start level the levels are whole and the container is
  // lattice arithmetic. A coarsened tile may straddle two tiles, and then both are marked.
  const phase::which_phase phases[2] = {phase::gas, phase::solid};

  for (int lvl = m_stopLevel - 1; lvl >= 0; lvl--) {
    for (const phase::which_phase& curPhase : phases) {
      const Vector<GeometryService::InOut>& fineTypes = this->types(curPhase)[lvl + 1];
      Vector<GeometryService::InOut>&       coarTypes = this->types(curPhase)[lvl];

      for (int i = 0; i < m_boxes[lvl + 1].size(); i++) {
        if (fineTypes[i] != GeometryService::Irregular) {
          continue;
        }

        const Box coarsened = coarsen(m_boxes[lvl + 1][i], 2);

        if (lvl > m_startLevel) {
          const Vector<int> containers = this->tilesMeeting(lvl, coarsened);

          if (containers.size() == 0) {
            MayDay::Error("ComputationalGeometry::buildCoarserLevels - an irregular tile has no tile beneath it");
          }

          // The tiles come first in a decimated level's list, so a tile index is its index in m_boxes.
          for (int j = 0; j < containers.size(); j++) {
            coarTypes[containers[j]] = GeometryService::Irregular;
          }
        }
        else {
          coarTypes[this->containingBox(lvl, coarsened.smallEnd())] = GeometryService::Irregular;
        }
      }
    }
  }
}

void
ComputationalGeometry::reportGrids() const
{
  CH_TIME("ComputationalGeometry::reportGrids");
  if (m_verbose) {
    pout() << "ComputationalGeometry::reportGrids" << endl;
  }

  pout() << "ComputationalGeometry::makeGrids - levels " << m_domains.size() << ", start level " << m_startLevel
         << ", stop level " << m_stopLevel << endl;

  for (int lvl = 0; lvl < m_domains.size(); lvl++) {
    int gasCount[3]   = {0, 0, 0};
    int solidCount[3] = {0, 0, 0};

    for (int i = 0; i < m_boxes[lvl].size(); i++) {
      gasCount[static_cast<int>(m_gasTypes[lvl][i])]++;
      solidCount[static_cast<int>(m_solidTypes[lvl][i])]++;
    }

    pout() << "  level " << lvl << " domain " << m_domains[lvl].domainBox().size() << " dx " << m_dx[lvl] << ": boxes "
           << m_boxes[lvl].size() << ", tiles " << m_cutTiles[lvl].size() << "; gas regular/covered/irregular "
           << gasCount[GeometryService::Regular] << "/" << gasCount[GeometryService::Covered] << "/"
           << gasCount[GeometryService::Irregular] << "; solid " << solidCount[GeometryService::Regular] << "/"
           << solidCount[GeometryService::Covered] << "/" << solidCount[GeometryService::Irregular]
           << "; split interior/ring/medial/doubled " << m_splitCounts[lvl][1] << "/" << m_splitCounts[lvl][2] << "/"
           << m_splitCounts[lvl][3] << "/" << m_splitCounts[lvl][4] << " (leaves " << m_splitCounts[lvl][0]
           << "); tiles tagged/nesting " << m_splitCounts[lvl][5] << "/" << m_splitCounts[lvl][6] << endl;
  }
}

void
ComputationalGeometry::buildGasGeometry(GeometryService*&    a_geoserver,
                                        const ProblemDomain& a_finestDomain,
                                        const RealVect&      a_probLo,
                                        const Real           a_finestDx)
{
  CH_TIME("ComputationalGeometry::buildGasGeometry(GeometryService, ProblemDomain, RealVect, Real)");
  if (m_verbose) {
    pout() << "ComputationalGeometry::buildGasGeometry(GeometryService, ProblemDomain, RealVect, Real)" << endl;
  }

  // Build the EBIS geometry. Use ScanShop, the polyhedral generator, or Chombo here. The polyhedral generator
  // computes its moments over the full grids ScanShop builds; the grids makeGrids made are not yet what the
  // index space is generated over.
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
    shop->setGrids(*this, phase::gas);
    shop->buildGraphs();
    shop->verifySurface();

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
  if (m_verbose) {
    pout() << "ComputationalGeometry::buildSolidGeometry(GeometryService, ProblemDomain, RealVect, Real)" << endl;
  }

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
      shop->setGrids(*this, phase::solid);
      shop->buildGraphs();
      shop->verifySurface();

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
