/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_AMRPolyhedralEBGraph.cpp
 * @brief  Implementation of CD_AMRPolyhedralEBGraph.H
 * @author Robert Marskar
 */

// Chombo includes
#include <CH_assert.H>
#include <CH_Timer.H>

// Our includes
#include <CD_AMRPolyhedralEBGraph.H>
#include <CD_NamespaceHeader.H>

AMRPolyhedralEBGraph::AMRPolyhedralEBGraph() noexcept
{
  m_isDefined       = false;
  m_ghost           = 0;
  m_numSourceLevels = 0;
}

AMRPolyhedralEBGraph::~AMRPolyhedralEBGraph() noexcept
{}

void
AMRPolyhedralEBGraph::define(const std::vector<DisjointBoxLayout>& a_grids,
                             const std::vector<ProblemDomain>&     a_domains,
                             const std::vector<int>&               a_geometryLevels,
                             const int                             a_numSourceLevels,
                             const int                             a_ghost)
{
  CH_TIME("AMRPolyhedralEBGraph::define");

  CH_assert(a_grids.size() == a_domains.size());
  CH_assert(a_grids.size() == a_geometryLevels.size());

  m_grids           = a_grids;
  m_domains         = a_domains;
  m_geometryLevels  = a_geometryLevels;
  m_numSourceLevels = a_numSourceLevels;
  m_ghost           = a_ghost;

  const int numLevels = static_cast<int>(a_grids.size());

  m_cellStates.assign(numLevels, {});
  m_surfaces.assign(numLevels, {});
  m_faceOverrides.assign(numLevels, {});

  for (int i = 0; i < numLevels; i++) {
    m_cellStates[i].resize(a_numSourceLevels);
    m_surfaces[i].resize(a_numSourceLevels);
    m_faceOverrides[i].resize(a_numSourceLevels);

    for (int lvl = 0; lvl < a_numSourceLevels; lvl++) {
      m_cellStates[i][lvl] = RefCountedPtr<LayoutData<BaseFab<signed char>>>(
        new LayoutData<BaseFab<signed char>>(a_grids[i]));
      m_surfaces[i][lvl] = RefCountedPtr<LayoutData<IVSFAB<PolyhedralEB::CutCellSurface>>>(
        new LayoutData<IVSFAB<PolyhedralEB::CutCellSurface>>(a_grids[i]));
      m_faceOverrides[i][lvl] = RefCountedPtr<LayoutData<PolyhedralEB::CutCellFaceOverrides>>(
        new LayoutData<PolyhedralEB::CutCellFaceOverrides>(a_grids[i]));
    }
  }

  m_isDefined = true;
}

bool
AMRPolyhedralEBGraph::isDefined() const noexcept
{
  return m_isDefined;
}

int
AMRPolyhedralEBGraph::getNumLevels() const noexcept
{
  return static_cast<int>(m_grids.size());
}

int
AMRPolyhedralEBGraph::getNumSourceLevels() const noexcept
{
  return m_numSourceLevels;
}

int
AMRPolyhedralEBGraph::getGhost() const noexcept
{
  return m_ghost;
}

const DisjointBoxLayout&
AMRPolyhedralEBGraph::getGrids(const int a_level) const noexcept
{
  return m_grids[a_level];
}

const ProblemDomain&
AMRPolyhedralEBGraph::getDomain(const int a_level) const noexcept
{
  return m_domains[a_level];
}

int
AMRPolyhedralEBGraph::getGeometryLevel(const int a_level) const noexcept
{
  return m_geometryLevels[a_level];
}

LayoutData<BaseFab<signed char>>&
AMRPolyhedralEBGraph::getCellStates(const int a_level, const int a_sourceLevel) noexcept
{
  return *m_cellStates[a_level][a_sourceLevel];
}

const LayoutData<BaseFab<signed char>>&
AMRPolyhedralEBGraph::getCellStates(const int a_level, const int a_sourceLevel) const noexcept
{
  return *m_cellStates[a_level][a_sourceLevel];
}

LayoutData<IVSFAB<PolyhedralEB::CutCellSurface>>&
AMRPolyhedralEBGraph::getSurfaces(const int a_level, const int a_sourceLevel) noexcept
{
  return *m_surfaces[a_level][a_sourceLevel];
}

const LayoutData<IVSFAB<PolyhedralEB::CutCellSurface>>&
AMRPolyhedralEBGraph::getSurfaces(const int a_level, const int a_sourceLevel) const noexcept
{
  return *m_surfaces[a_level][a_sourceLevel];
}

LayoutData<PolyhedralEB::CutCellFaceOverrides>&
AMRPolyhedralEBGraph::getFaceOverrides(const int a_level, const int a_sourceLevel) noexcept
{
  return *m_faceOverrides[a_level][a_sourceLevel];
}

const LayoutData<PolyhedralEB::CutCellFaceOverrides>&
AMRPolyhedralEBGraph::getFaceOverrides(const int a_level, const int a_sourceLevel) const noexcept
{
  return *m_faceOverrides[a_level][a_sourceLevel];
}

#include <CD_NamespaceFooter.H>
