/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_AMRPolyhedralEBGraph.H
 * @brief  Declaration of the polyhedral graph as an AMR hierarchy of any layout holds it.
 * @author Robert Marskar
 */

#ifndef CD_AMRPOLYHEDRALEBGRAPH_H
#define CD_AMRPOLYHEDRALEBGRAPH_H

// Std includes
#include <vector>

// Chombo includes
#include <BaseFab.H>
#include <DisjointBoxLayout.H>
#include <IVSFAB.H>
#include <LayoutData.H>
#include <ProblemDomain.H>
#include <RefCountedPtr.H>

// Our includes
#include <CD_CutCellFaceOverrides.H>
#include <CD_CutCellSurface.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief The polyhedral graph as an AMR hierarchy of any layout holds it: every cell of every box, ghost cells
 * included, described by the geometry level that owns it.
 * @details The geometry's graph lives on its own levels and its own tiles. A consumer lives on its own hierarchy --
 * one DisjointBoxLayout per level, each box filled to a ghost width, less what the next finer level covers -- and a
 * cell of it may be owned by the geometry level at its own resolution, by a finer one, or by a coarser one. This holds,
 * for every box of every destination level and every geometry level that owns part of it, the data that level's cells
 * are rebuilt from: their states, their surfaces and their face overrides, as that level stores them and in its index
 * space. Nothing here is recomputed from the implicit function; PolyhedralGeometryShop::copyGraphs fills it from the
 * graphs alone, and a body is rebuilt from it by PolyhedralGeometryShop::defineBody exactly as on the graph itself.
 *
 * The data of a (destination level, geometry level) pair is kept per destination box, over the box of cells of the
 * geometry level that the destination box's pieces span. A cell there that no piece brought reads
 * PolyhedralGeometryShop::s_absent.
 */
class AMRPolyhedralEBGraph
{
public:
  /**
   * @brief Constructor. Holds nothing until define.
   */
  AMRPolyhedralEBGraph() noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~AMRPolyhedralEBGraph() noexcept;

  /**
   * @brief Set the destination hierarchy and allocate an empty set of containers for every pair of destination and
   * geometry level.
   * @param[in] a_grids           Destination layout of every level, coarsest first.
   * @param[in] a_domains         Domain of every destination level.
   * @param[in] a_geometryLevels  Geometry level of every destination level.
   * @param[in] a_numSourceLevels Number of geometry levels.
   * @param[in] a_ghost           Ghost width every destination box is filled to.
   */
  void
  define(const std::vector<DisjointBoxLayout>& a_grids,
         const std::vector<ProblemDomain>&     a_domains,
         const std::vector<int>&               a_geometryLevels,
         const int                             a_numSourceLevels,
         const int                             a_ghost);

  /**
   * @brief Whether define has run.
   * @return True after define.
   */
  bool
  isDefined() const noexcept;

  /**
   * @brief Number of destination levels.
   * @return The count.
   */
  int
  getNumLevels() const noexcept;

  /**
   * @brief Number of geometry levels.
   * @return The count.
   */
  int
  getNumSourceLevels() const noexcept;

  /**
   * @brief Ghost width every destination box is filled to.
   * @return The width.
   */
  int
  getGhost() const noexcept;

  /**
   * @brief Layout of a destination level.
   * @param[in] a_level Destination level.
   * @return The layout.
   */
  const DisjointBoxLayout&
  getGrids(const int a_level) const noexcept;

  /**
   * @brief Domain of a destination level.
   * @param[in] a_level Destination level.
   * @return The domain.
   */
  const ProblemDomain&
  getDomain(const int a_level) const noexcept;

  /**
   * @brief Geometry level a destination level lives on.
   * @param[in] a_level Destination level.
   * @return The geometry level.
   */
  int
  getGeometryLevel(const int a_level) const noexcept;

  /**
   * @brief Cell states of the cells a geometry level owns in every box of a destination level, in the geometry
   * level's index space; s_absent where it owns none.
   * @param[in] a_level       Destination level.
   * @param[in] a_sourceLevel Geometry level.
   * @return One fab per destination box, empty where the geometry level owns nothing of it.
   */
  LayoutData<BaseFab<signed char>>&
  getCellStates(const int a_level, const int a_sourceLevel) noexcept;

  /**
   * @brief Cell states, read only.
   * @param[in] a_level       Destination level.
   * @param[in] a_sourceLevel Geometry level.
   * @return One fab per destination box.
   */
  const LayoutData<BaseFab<signed char>>&
  getCellStates(const int a_level, const int a_sourceLevel) const noexcept;

  /**
   * @brief Surfaces of the cells that hold one, of those a geometry level owns in every box of a destination level.
   * @param[in] a_level       Destination level.
   * @param[in] a_sourceLevel Geometry level.
   * @return One set of surfaces per destination box.
   */
  LayoutData<IVSFAB<PolyhedralEB::CutCellSurface>>&
  getSurfaces(const int a_level, const int a_sourceLevel) noexcept;

  /**
   * @brief Surfaces, read only.
   * @param[in] a_level       Destination level.
   * @param[in] a_sourceLevel Geometry level.
   * @return One set of surfaces per destination box.
   */
  const LayoutData<IVSFAB<PolyhedralEB::CutCellSurface>>&
  getSurfaces(const int a_level, const int a_sourceLevel) const noexcept;

  /**
   * @brief Face overrides of the cells whose surfaces are overridden, of those a geometry level owns in every box of
   * a destination level.
   * @param[in] a_level       Destination level.
   * @param[in] a_sourceLevel Geometry level.
   * @return One set of overrides per destination box.
   */
  LayoutData<PolyhedralEB::CutCellFaceOverrides>&
  getFaceOverrides(const int a_level, const int a_sourceLevel) noexcept;

  /**
   * @brief Face overrides, read only.
   * @param[in] a_level       Destination level.
   * @param[in] a_sourceLevel Geometry level.
   * @return One set of overrides per destination box.
   */
  const LayoutData<PolyhedralEB::CutCellFaceOverrides>&
  getFaceOverrides(const int a_level, const int a_sourceLevel) const noexcept;

protected:
  /**
   * @brief Whether define has run.
   */
  bool m_isDefined;

  /**
   * @brief Ghost width every destination box is filled to.
   */
  int m_ghost;

  /**
   * @brief Number of geometry levels.
   */
  int m_numSourceLevels;

  /**
   * @brief Destination layout of every level.
   */
  std::vector<DisjointBoxLayout> m_grids;

  /**
   * @brief Domain of every destination level.
   */
  std::vector<ProblemDomain> m_domains;

  /**
   * @brief Geometry level of every destination level.
   */
  std::vector<int> m_geometryLevels;

  /**
   * @brief Cell states, per destination level and geometry level.
   */
  std::vector<std::vector<RefCountedPtr<LayoutData<BaseFab<signed char>>>>> m_cellStates;

  /**
   * @brief Surfaces, per destination level and geometry level.
   */
  std::vector<std::vector<RefCountedPtr<LayoutData<IVSFAB<PolyhedralEB::CutCellSurface>>>>> m_surfaces;

  /**
   * @brief Face overrides, per destination level and geometry level.
   */
  std::vector<std::vector<RefCountedPtr<LayoutData<PolyhedralEB::CutCellFaceOverrides>>>> m_faceOverrides;
};

#include <CD_NamespaceFooter.H>

#endif
