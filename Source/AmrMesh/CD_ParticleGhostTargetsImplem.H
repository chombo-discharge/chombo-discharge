/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleGhostTargetsImplem.H
  @brief  Implementation of CD_ParticleGhostTargets.H
  @author Robert Marskar
*/

#ifndef CD_PARTICLEGHOSTTARGETSIMPLEM_H
#define CD_PARTICLEGHOSTTARGETSIMPLEM_H

// Our includes
#include <CD_ParticleGhostTargets.H>
#include <CD_NamespaceHeader.H>

inline void
ParticleGhostTargets::define(const Box& a_box) noexcept
{
  m_box = a_box;

  m_count.resize(a_box, 1);
  m_count.setVal(0);
  m_start.resize(a_box, 1);
  m_start.setVal(-1);
  m_flat.clear();
}

inline void
ParticleGhostTargets::incrementCount(const IntVect& a_iv) noexcept
{
  m_count(a_iv, 0) += 1;
}

inline void
ParticleGhostTargets::allocate() noexcept
{
  int offset = 0;
  for (BoxIterator bit(m_box); bit.ok(); ++bit) {
    const int c = m_count(bit(), 0);
    if (c > 0) {
      m_start(bit(), 0) = offset;
      offset += c;
    }
  }
  m_flat.resize(offset);
}

inline void
ParticleGhostTargets::addTarget(const IntVect& a_iv, const LevelTiles::BoxIDs& a_target) noexcept
{
  CH_assert(m_start(a_iv, 0) >= 0);
  CH_assert(m_start(a_iv, 0) < static_cast<int>(m_flat.size()));

  // m_start doubles as the per-cell write head during the fill; finalize() rewinds it.
  m_flat[m_start(a_iv, 0)] = a_target;
  m_start(a_iv, 0) += 1;
}

inline void
ParticleGhostTargets::finalize() noexcept
{
  for (BoxIterator bit(m_box); bit.ok(); ++bit) {
    m_start(bit(), 0) -= m_count(bit(), 0);
  }
}

inline const Box&
ParticleGhostTargets::box() const noexcept
{
  return m_box;
}

inline int
ParticleGhostTargets::numTargets(const IntVect& a_iv) const noexcept
{
  return m_count(a_iv, 0);
}

inline const LevelTiles::BoxIDs&
ParticleGhostTargets::target(const IntVect& a_iv, const int a_i) const noexcept
{
  CH_assert(a_i >= 0 && a_i < m_count(a_iv, 0));

  return m_flat[m_start(a_iv, 0) + a_i];
}

inline unsigned int
ParticleGhostTargets::targetGridIndex(const IntVect& a_iv, const int a_i) const noexcept
{
  return this->target(a_iv, a_i).first;
}

inline unsigned int
ParticleGhostTargets::targetRank(const IntVect& a_iv, const int a_i) const noexcept
{
  return this->target(a_iv, a_i).second;
}

#include <CD_NamespaceFooter.H>

#endif
