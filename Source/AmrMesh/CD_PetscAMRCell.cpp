/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file CD_PetscAMRCell.cpp
  @brief Implementation of CD_PetscAMRCell.H
  @author Robert Marskar
*/

// Our includes
#include <CD_PetscAMRCell.H>
#include <CD_NamespaceHeader.H>

#warning \
  "I want new flags for cells that are on the fine side of the CF, and on the coarse side of the CF (in addition to the ghostCF flag)"

PetscAMRCell::PetscAMRCell() noexcept
{
  m_numPhases            = 0;
  m_isCoveredByFinerGrid = false;
  m_isCoarCF             = false;
  m_isFine               = false;
  m_isGhostCF            = false;
  m_isDomainBoundaryCell = false;
}

PetscAMRCell::~PetscAMRCell() noexcept
{}

void
PetscAMRCell::define(const unsigned int a_numPhases,
                     const bool         a_isCoveredByFinerGrid,
                     const bool         a_isCoarCF,
                     const bool         a_isFineCF,
                     const bool         a_isGhostCF,
                     const bool         a_isDomainBoundaryCell) noexcept
{
  m_numPhases            = a_numPhases;
  m_isCoveredByFinerGrid = a_isCoveredByFinerGrid;
  m_isCoarCF             = a_isCoarCF;
  m_isFineCF             = a_isFineCF;
  m_isGhostCF            = a_isGhostCF;
  m_isDomainBoundaryCell = a_isDomainBoundaryCell;
}

void
PetscAMRCell::setNumPhases(const unsigned int a_numPhases) noexcept
{
  m_numPhases = a_numPhases;
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

unsigned int
PetscAMRCell::getNumPhases() const noexcept
{
  return static_cast<unsigned int>(m_numPhases);
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

int
linearSize(const PetscAMRCell& a_amrCell)
{
  int size = 0;

  size += sizeof(PetscAMRCell::m_numPhases);
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
  const unsigned char* buffer = static_cast<const unsigned char*>(a_buffer);

  a_amrCell.m_numPhases = *buffer;
  buffer++;

  a_amrCell.m_isCoveredByFinerGrid = *buffer;
  buffer++;

  a_amrCell.m_isGhostCF = *buffer;
  buffer++;

  a_amrCell.m_isCoarCF = *buffer;
  buffer++;

  a_amrCell.m_isFineCF = *buffer;
  buffer++;

  a_amrCell.m_isDomainBoundaryCell = *buffer;
  buffer++;
}

void
linearOut(void* const a_buffer, const PetscAMRCell& a_amrCell)
{
  unsigned char* buffer = static_cast<unsigned char*>(a_buffer);

  *buffer = a_amrCell.m_numPhases;
  buffer++;

  *buffer = static_cast<unsigned char>(a_amrCell.m_isCoveredByFinerGrid);
  buffer++;

  *buffer = static_cast<unsigned char>(a_amrCell.m_isCoarCF);
  buffer++;

  *buffer = static_cast<unsigned char>(a_amrCell.m_isFineCF);
  buffer++;

  *buffer = static_cast<unsigned char>(a_amrCell.m_isGhostCF);
  buffer++;

  *buffer = static_cast<unsigned char>(a_amrCell.m_isDomainBoundaryCell);
  buffer++;
}

#include <CD_NamespaceFooter.H>
