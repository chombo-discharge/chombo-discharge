/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file CD_PetscAMRCell.cpp
 * @brief Implementation of CD_PetscAMRCell.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Std includes
#include <cstdint>
#include <cstring>
#include <type_traits>

// Our includes
#include <CD_PetscAMRCell.H>
#include <CD_NamespaceHeader.H>

namespace {
/**
 * @brief Copy a trivially-copyable value into a byte buffer and advance the buffer pointer.
 * @param[in,out] a_p Buffer pointer. Advanced by sizeof(T).
 * @param[in] a_x Value to write.
 */
template <class T>
void
pushValue(uint8_t*& a_p, const T& a_x) noexcept
{
  static_assert(std::is_trivially_copyable<T>::value, "pushValue expects trivially-copyable values");

  std::memcpy(a_p, &a_x, sizeof(T));

  a_p += sizeof(T);
}

/**
 * @brief Copy a trivially-copyable value out of a byte buffer and advance the buffer pointer.
 * @param[in,out] a_p Buffer pointer. Advanced by sizeof(T).
 * @param[out] a_x Value to read into.
 */
template <class T>
void
pullValue(const uint8_t*& a_p, T& a_x) noexcept
{
  static_assert(std::is_trivially_copyable<T>::value, "pullValue expects trivially-copyable values");

  std::memcpy(&a_x, a_p, sizeof(T));

  a_p += sizeof(T);
}
} // namespace

PetscAMRCell::PetscAMRCell() noexcept
{
  m_petscRows[0]         = -1;
  m_petscRows[1]         = -1;
  m_isCoveredByFinerGrid = false;
  m_isCoarCF             = false;
  m_isFineCF             = false;
  m_isGhostCF            = false;
  m_isDomainBoundaryCell = false;
}

PetscAMRCell::~PetscAMRCell() noexcept
{}

void
PetscAMRCell::define(const PetscInt a_petscRowPhase0,
                     const PetscInt a_petscRowPhase1,
                     const bool     a_isCoveredByFinerGrid,
                     const bool     a_isCoarCF,
                     const bool     a_isFineCF,
                     const bool     a_isGhostCF,
                     const bool     a_isDomainBoundaryCell) noexcept
{
  m_petscRows[0]         = a_petscRowPhase0;
  m_petscRows[1]         = a_petscRowPhase1;
  m_isCoveredByFinerGrid = a_isCoveredByFinerGrid;
  m_isCoarCF             = a_isCoarCF;
  m_isFineCF             = a_isFineCF;
  m_isGhostCF            = a_isGhostCF;
  m_isDomainBoundaryCell = a_isDomainBoundaryCell;
}

void
PetscAMRCell::setPetscRow(const int a_phase, const PetscInt a_row) noexcept
{
  CH_assert(a_phase == 0 || a_phase == 1);

  m_petscRows[a_phase] = a_row;
}

void
PetscAMRCell::setCoveredByFinerGrid(const bool a_coveredByFinerGrid) noexcept
{
  m_isCoveredByFinerGrid = a_coveredByFinerGrid;
}

void
PetscAMRCell::setCoarCF(const bool a_isCoarCF) noexcept
{
  m_isCoarCF = a_isCoarCF;
}

void
PetscAMRCell::setFineCF(const bool a_isFineCF) noexcept
{
  m_isFineCF = a_isFineCF;
}

void
PetscAMRCell::setGhostCF(const bool a_isGhostCF) noexcept
{
  m_isGhostCF = a_isGhostCF;
}

void
PetscAMRCell::setDomainBoundaryCell(const bool a_isDomainBoundaryCell) noexcept
{
  m_isDomainBoundaryCell = a_isDomainBoundaryCell;
}

PetscInt
PetscAMRCell::getPetscRow(const int a_iphase) const noexcept
{
  return m_petscRows[a_iphase];
}

bool
PetscAMRCell::isCoveredByFinerGrid() const noexcept
{
  return m_isCoveredByFinerGrid;
}

bool
PetscAMRCell::isCoarCF() const noexcept
{
  return m_isCoarCF;
}

bool
PetscAMRCell::isFineCF() const noexcept
{
  return m_isFineCF;
}

bool
PetscAMRCell::isGhostCF() const noexcept
{
  return m_isGhostCF;
}

bool
PetscAMRCell::isDomainBoundaryCell() const noexcept
{
  return m_isDomainBoundaryCell;
}

constexpr int
linearSize(const PetscAMRCell& a_amrCell)
{
  int size = 0;

  size += sizeof(PetscAMRCell::m_petscRows);
  size += sizeof(PetscAMRCell::m_isCoveredByFinerGrid);
  size += sizeof(PetscAMRCell::m_isCoarCF);
  size += sizeof(PetscAMRCell::m_isFineCF);
  size += sizeof(PetscAMRCell::m_isGhostCF);
  size += sizeof(PetscAMRCell::m_isDomainBoundaryCell);

  return size;
}

void
linearIn(PetscAMRCell& a_amrCell, const void* const a_buffer)
{
  const uint8_t* p = static_cast<const uint8_t*>(a_buffer);

  pullValue(p, a_amrCell.m_petscRows[0]);
  pullValue(p, a_amrCell.m_petscRows[1]);
  pullValue(p, a_amrCell.m_isCoveredByFinerGrid);
  pullValue(p, a_amrCell.m_isGhostCF);
  pullValue(p, a_amrCell.m_isCoarCF);
  pullValue(p, a_amrCell.m_isFineCF);
  pullValue(p, a_amrCell.m_isDomainBoundaryCell);
}

void
linearOut(void* const a_buffer, const PetscAMRCell& a_amrCell)
{
  uint8_t* p = static_cast<uint8_t*>(a_buffer);

  pushValue(p, a_amrCell.m_petscRows[0]);
  pushValue(p, a_amrCell.m_petscRows[1]);
  pushValue(p, a_amrCell.m_isCoveredByFinerGrid);
  pushValue(p, a_amrCell.m_isGhostCF);
  pushValue(p, a_amrCell.m_isCoarCF);
  pushValue(p, a_amrCell.m_isFineCF);
  pushValue(p, a_amrCell.m_isDomainBoundaryCell);
}

#include <CD_NamespaceFooter.H>

#endif
