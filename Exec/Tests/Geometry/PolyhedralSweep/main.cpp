/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Sweep a rotated cube and cylinder across the grid and check the polyhedral geometry on each pose.
 * @author Robert Marskar
 */

// Std includes
#include <string>
#include <vector>

// Chombo includes
#include <ParmParse.H>
#include <TransformIF.H>

// Our includes
#include <CD_BoxSdf.H>
#include <CD_CylinderSdf.H>
#include <CD_ComputationalGeometry.H>
#include <CD_Initialize.H>
#include <CD_ParallelOps.H>
#include <CD_Electrode.H>
#include <CD_PolyhedralEBGraph.H>
#include <CD_PolyhedralGeometryShop.H>
#include <CD_SphereSdf.H>
#include <CD_Timer.H>

using namespace ChomboDischarge;

/**
 * @brief One pose of one shape: which shape, how far it is shifted, and how far it is turned.
 */
struct Pose
{
  std::string shape;
  Real        shift;
  Real        angle;
};

/**
 * @brief Build the implicit function of a pose, centred on the domain and turned about its centre.
 * @details The cube is the shape whose faces and edges can align with the grid, and the cylinder the one whose
 * surface is curved in one direction and flat in the other; in two dimensions they are a square and a circle,
 * the circle being what a cylinder's cross-section is. The shift is in cells of the finest level and moves the
 * shape along the diagonal, so that a face, an edge and a corner each pass through nodes, edges and cell
 * interiors as the sweep runs.
 * @param[in] a_pose   The pose.
 * @param[in] a_centre Centre of the domain.
 * @param[in] a_size   Half-width of the shape.
 * @param[in] a_dx     Grid spacing of the finest level.
 * @return The implicit function, fluid outside the shape.
 */
RefCountedPtr<BaseIF>
buildShape(const Pose& a_pose, const RealVect& a_centre, const Real a_size, const Real a_dx)
{
  const RealVect centre = a_centre + a_pose.shift * a_dx * RealVect::Unit;

  BaseIF* shape = nullptr;

  if (a_pose.shape == "cube") {
    shape = new BoxSdf(centre - a_size * RealVect::Unit, centre + a_size * RealVect::Unit, false);
  }
  else if (a_pose.shape == "cylinder") {
#if CH_SPACEDIM == 3
    shape = new CylinderSdf(centre - a_size * RealVect(BASISV(2)),
                            centre + a_size * RealVect(BASISV(2)),
                            a_size,
                            false);
#else
    shape = new SphereSdf(centre, a_size, false);
#endif
  }
  else {
    MayDay::Error("PolyhedralSweep - unknown shape, expected cube or cylinder");
  }

  // The turn is about the domain's centre. A cylinder turned about its own axis is the cylinder again, so it is
  // turned about the first direction instead, which tilts the axis out of the grid.
  auto* turned = new TransformIF(*shape);

  turned->rotate(a_pose.angle * M_PI / 180.0,
                 centre,
                 RealVect(BASISV((a_pose.shape == "cylinder") ? 0 : SpaceDim - 1)));

  return RefCountedPtr<BaseIF>(static_cast<BaseIF*>(turned));
}

/**
 * @brief Build the geometry of one pose and check it.
 * @details Builds the grids, the polyhedral graph of every level from them, and runs the closure check over the
 * graphs, which is what fails if a cell's surface has a hole or if a face on a coarse-fine seam cannot be
 * described at the coarse level. Every failure stops the run from inside the generator, naming the level and the
 * cell; this reports the pose that was being built when it did.
 * @param[in] a_pose        The pose.
 * @param[in] a_nCells      Cells across the coarsest domain.
 * @param[in] a_depth       Levels above the coarsest domain.
 * @param[in] a_refineAngle Angle between neighbouring normals above which a box is split, in degrees.
 * @param[in] a_probLo      Lower-left corner of the domain.
 * @param[in] a_probHi      Upper-right corner of the domain.
 * @param[in] a_maxGhostEB  Ghost cells a box is grown by when it is classified.
 * @return Cut cells over every level of the gas phase, which is what the checks ran over.
 */
long long
checkPose(const Pose&    a_pose,
          const int      a_nCells,
          const int      a_depth,
          const Real     a_refineAngle,
          const RealVect a_probLo,
          const RealVect a_probHi,
          const int      a_maxGhostEB)
{
  const Box           coarseBox(IntVect::Zero, (a_nCells - 1) * IntVect::Unit);
  const ProblemDomain coarseDomain(coarseBox);
  const ProblemDomain fineDomain = refine(coarseDomain, static_cast<int>(std::pow(2, a_depth)));

  const Real coarseDx = (a_probHi[0] - a_probLo[0]) / static_cast<Real>(a_nCells);
  const Real fineDx   = coarseDx / std::pow(2, a_depth);

  const RealVect centre = 0.5 * (a_probLo + a_probHi);

  RefCountedPtr<BaseIF> shape = buildShape(a_pose, centre, 0.25 * (a_probHi[0] - a_probLo[0]), fineDx);

  ComputationalGeometry compGeom;

  compGeom.setElectrodes(Vector<Electrode>(1, Electrode(shape, true, 1.0)));

  compGeom.makeGrids(coarseDomain, fineDomain, a_probLo, coarseDx, a_refineAngle, a_maxGhostEB);

  PolyhedralGeometryShop
    shop(*compGeom.getGasImplicitFunction(), 0, fineDx, a_probLo, fineDomain, coarseDomain, a_maxGhostEB, 1.E-15, true);

  shop.setGrids(compGeom, phase::gas);
  shop.buildGraphs();
  shop.verifySurface();

  // Every level with cut tiles must have a graph, and the graphs must hold cut cells: a pose whose surface the
  // builder lost would otherwise pass every check by having nothing to check.
  long long numCutCells = 0;

  for (int lvl = compGeom.getStartLevel(); lvl < compGeom.getNumGridLevels(); lvl++) {
    const bool hasTiles = compGeom.getCutTiles(lvl).size() > 0;

    const PolyhedralEBGraph& graph = shop.getGraph(lvl);

    if (hasTiles != graph.isDefined()) {
      MayDay::Error("PolyhedralSweep - a level's cut tiles and its graph disagree");
    }

    if (!graph.isDefined()) {
      continue;
    }

    const DisjointBoxLayout& grids = graph.getGrids();

    for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box box = grids[dit()];

      const IntVectSet& cut = graph.getCutCells()[dit()];

      for (BoxIterator bit(box); bit.ok(); ++bit) {
        if (cut.contains(bit())) {
          numCutCells++;
        }
      }
    }
  }

  return ParallelOps::sum(numCutCells);
}

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  int  nCells      = 32;
  int  depth       = 2;
  int  maxGhostEB  = 2;
  Real refineAngle = 15.0;

  Vector<Real>        shifts;
  Vector<Real>        angles;
  Vector<std::string> shapes;

  Vector<Real> probLo(SpaceDim, 0.0);
  Vector<Real> probHi(SpaceDim, 1.0);

  ParmParse pp("PolyhedralSweep");

  pp.get("n_cells", nCells);
  pp.get("depth", depth);
  pp.get("refine_angle", refineAngle);
  pp.get("max_ghost_eb", maxGhostEB);
  pp.getarr("prob_lo", probLo, 0, SpaceDim);
  pp.getarr("prob_hi", probHi, 0, SpaceDim);
  pp.getarr("shifts", shifts, 0, pp.countval("shifts"));
  pp.getarr("angles", angles, 0, pp.countval("angles"));
  pp.getarr("shapes", shapes, 0, pp.countval("shapes"));

  RealVect lo = RealVect::Zero;
  RealVect hi = RealVect::Zero;

  for (int dir = 0; dir < SpaceDim; dir++) {
    lo[dir] = probLo[dir];
    hi[dir] = probHi[dir];
  }

  Timer timer("PolyhedralSweep");

  int numPoses = 0;

  for (int s = 0; s < shapes.size(); s++) {
    for (int a = 0; a < angles.size(); a++) {
      for (int t = 0; t < shifts.size(); t++) {
        const Pose pose{shapes[s], shifts[t], angles[a]};

        const std::string what = pose.shape + " turned " + std::to_string(pose.angle) + " degrees, shifted " +
                                 std::to_string(pose.shift) + " cells";

        // Announced before the pose is built, so that a generator that stops the run names the pose it stopped on.
        if (procID() == 0) {
          pout() << "PolyhedralSweep - " << what << endl;
        }

        timer.startEvent(pose.shape);
        const long long numCutCells = checkPose(pose, nCells, depth, refineAngle, lo, hi, maxGhostEB);
        timer.stopEvent(pose.shape);

        if (numCutCells == 0) {
          MayDay::Error("PolyhedralSweep - a pose produced no cut cells at all");
        }

        // The closure check is three-dimensional -- a two-dimensional interface is one chord per cell, which the
        // graph checks as it builds it -- so a two-dimensional pose is checked by the build alone.
        if (procID() == 0) {
#if CH_SPACEDIM == 3
          pout() << "PolyhedralSweep - " << what << ": " << numCutCells << " cut cells, surface closed" << endl;
#else
          pout() << "PolyhedralSweep - " << what << ": " << numCutCells << " cut cells, graph built" << endl;
#endif
        }

        numPoses++;
      }
    }
  }

  if (procID() == 0) {
    pout() << "PolyhedralSweep - " << numPoses << " poses built and checked" << endl;
  }

  timer.eventReport(pout(), false);

  ChomboDischarge::finalize();
}
