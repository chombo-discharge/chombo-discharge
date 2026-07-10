/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MergeBlockGrid.cpp
  @brief  Implementation of CD_MergeBlockGrid.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <DataIterator.H>

// Our includes
#include <CD_MergeBlockGrid.H>
#include <CD_NamespaceHeader.H>

MergeBlockGrid::MergeBlockGrid() noexcept : m_isDefined(false), m_blockSize(0)
{
  CH_TIME("MergeBlockGrid::MergeBlockGrid(weak)");
}

MergeBlockGrid::MergeBlockGrid(const DisjointBoxLayout& a_dbl, const int a_blockSize) noexcept
{
  CH_TIME("MergeBlockGrid::MergeBlockGrid(full)");

  this->define(a_dbl, a_blockSize);
}

MergeBlockGrid::~MergeBlockGrid() noexcept
{
  CH_TIME("MergeBlockGrid::~MergeBlockGrid");
}

void
MergeBlockGrid::define(const DisjointBoxLayout& a_dbl, const int a_blockSize) noexcept
{
  CH_TIME("MergeBlockGrid::define");

  m_blockSize = a_blockSize;

  m_patchBlocks.define(a_dbl);

  for (DataIterator dit(a_dbl); dit.ok(); ++dit) {
    PatchBlocks& pb = m_patchBlocks[dit()];

    pb.cellBox      = a_dbl[dit()];
    pb.blockBox     = coarsen(pb.cellBox, a_blockSize);
    pb.haloBlockBox = grow(pb.blockBox, 1);
  }

  m_isDefined = true;
}

int
MergeBlockGrid::getBlockSize() const noexcept
{
  CH_assert(m_isDefined);

  return m_blockSize;
}

const MergeBlockGrid::PatchBlocks&
MergeBlockGrid::getPatchBlocks(const DataIndex& a_din) const noexcept
{
  CH_assert(m_isDefined);

  return m_patchBlocks[a_din];
}

IntVect
MergeBlockGrid::blockIndexOf(const IntVect& a_cell, const int a_blockSize) noexcept
{
  IntVect block;
  for (int dir = 0; dir < SpaceDim; dir++) {
    block[dir] = (a_cell[dir] >= 0) ? (a_cell[dir] / a_blockSize) : -(((-a_cell[dir] - 1) / a_blockSize) + 1);
  }
  return block;
}

void
MergeBlockGrid::blockAABB(const IntVect&  a_blockIndex,
                          const int       a_blockSize,
                          const RealVect& a_probLo,
                          const RealVect& a_dx,
                          RealVect&       a_lo,
                          RealVect&       a_hi) noexcept
{
  for (int dir = 0; dir < SpaceDim; dir++) {
    a_lo[dir] = a_probLo[dir] + a_blockIndex[dir] * a_blockSize * a_dx[dir];
    a_hi[dir] = a_lo[dir] + a_blockSize * a_dx[dir];
  }
}

#include <CD_NamespaceFooter.H>
