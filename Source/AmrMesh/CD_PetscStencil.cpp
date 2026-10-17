/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscStencil.cpp
 * @brief  Implementation of CD_PetscStencil.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Std includes
#include <algorithm>

// Chombo includes
#include <CH_assert.H>

// Our includes
#include <CD_PetscStencil.H>
#include <CD_NamespaceHeader.H>

void
PetscStencil::clear()
{
  m_entries.resize(0);
}

int
PetscStencil::size() const noexcept
{
  return static_cast<int>(m_entries.size());
}

PetscInt
PetscStencil::row(const int a_entry) const noexcept
{
  CH_assert(a_entry >= 0 && a_entry < this->size());

  return m_entries[a_entry].first;
}

PetscScalar
PetscStencil::weight(const int a_entry) const noexcept
{
  CH_assert(a_entry >= 0 && a_entry < this->size());

  return m_entries[a_entry].second;
}

void
PetscStencil::add(const PetscInt& a_row, const PetscScalar& a_weight) noexcept
{
  m_entries.emplace_back(a_row, a_weight);
}

PetscStencil&
PetscStencil::operator+=(const PetscStencil& a_otherStencil) noexcept
{
  m_entries.insert(m_entries.end(), a_otherStencil.m_entries.begin(), a_otherStencil.m_entries.end());

  return *this;
}

void
PetscStencil::compress() noexcept
{
  // Stable, so duplicates are summed in insertion order and the result does not depend on the sort implementation.
  std::stable_sort(m_entries.begin(),
                   m_entries.end(),
                   [](const std::pair<PetscInt, PetscScalar>& a, const std::pair<PetscInt, PetscScalar>& b) -> bool {
                     return a.first < b.first;
                   });

  // Merge runs of equal rows into the first entry of each run.
  size_t numUnique = 0;

  for (size_t i = 0; i < m_entries.size(); i++) {
    if (numUnique > 0 && m_entries[numUnique - 1].first == m_entries[i].first) {
      m_entries[numUnique - 1].second += m_entries[i].second;
    }
    else {
      m_entries[numUnique] = m_entries[i];

      numUnique++;
    }
  }

  m_entries.resize(numUnique);
}

PetscStencil&
PetscStencil::operator*=(const PetscScalar& a_scalingFactor) noexcept
{
  for (auto& entry : m_entries) {
    entry.second *= a_scalingFactor;
  }

  return *this;
}

#include <CD_NamespaceFooter.H>

#endif
