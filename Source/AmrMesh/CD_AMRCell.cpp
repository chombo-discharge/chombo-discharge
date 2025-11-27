/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file CD_AMRCell.cpp
  @brief Implementation of CD_AMRCell.H
  @author Robert Marskar
*/

// Our includes
#include <CD_AMRCell.H>
#include <CD_NamespaceHeader.H>

AMRCell::AMRCell() noexcept
{
  m_numPhases            = 0;
  m_isCoveredByFinerGrid = false;
  m_isGhostCF            = false;
  m_isDomainCell         = false;
}

AMRCell::~AMRCell() noexcept
{}

void
AMRCell::define(const unsigned int a_numPhases,
                const bool         a_isCoveredByFinerGrid,
                const bool         a_isGhostCF,
                const bool         a_isDomainCell) noexcept
{
  m_numPhases            = a_numPhases;
  m_isCoveredByFinerGrid = a_isCoveredByFinerGrid;
  m_isGhostCF            = a_isGhostCF;
  m_isDomainCell         = a_isDomainCell;
}

void
AMRCell::setNumPhases(const unsigned int a_numPhases) noexcept
{
  m_numPhases = a_numPhases;
}

void
AMRCell::setCoveredByFinerGrid(const bool a_coveredByFinerGrid) noexcept
{
  m_isCoveredByFinerGrid = a_coveredByFinerGrid;
}

void
AMRCell::setGhostCF(const bool a_isGhostCF) noexcept
{
  m_isGhostCF = a_isGhostCF;
}

void
AMRCell::setDomainCell(const bool a_isDomainCell) noexcept
{
  m_isDomainCell = a_isDomainCell;
}

unsigned int
AMRCell::getNumPhases() const noexcept
{
  return static_cast<unsigned int>(m_numPhases);
}

bool
AMRCell::isCoveredByFinerGrid() const noexcept
{
  return m_isCoveredByFinerGrid;
}

bool
AMRCell::isGhostCF() const noexcept
{
  return m_isGhostCF;
}

bool
AMRCell::isDomainCell() const noexcept
{
  return m_isDomainCell;
}

int
linearSize(const AMRCell& a_amrCell)
{
  int size = 0;

  size += sizeof(AMRCell::m_numPhases);
  size += sizeof(AMRCell::m_isCoveredByFinerGrid);
  size += sizeof(AMRCell::m_isGhostCF);
  size += sizeof(AMRCell::m_isDomainCell);

  return size;
}

void
linearIn(AMRCell& a_amrCell, const void* const a_buffer)
{
  const unsigned char* buffer = static_cast<const unsigned char*>(a_buffer);

  a_amrCell.m_numPhases = *buffer;
  buffer++;

  a_amrCell.m_isCoveredByFinerGrid = *buffer;
  buffer++;

  a_amrCell.m_isGhostCF = *buffer;
  buffer++;

  a_amrCell.m_isDomainCell = *buffer;
  buffer++;
}

void
linearOut(void* const a_buffer, const AMRCell& a_amrCell)
{
  unsigned char* buffer = static_cast<unsigned char*>(a_buffer);

  *buffer = a_amrCell.m_numPhases;
  buffer++;

  *buffer = static_cast<unsigned char>(a_amrCell.m_isCoveredByFinerGrid);
  buffer++;

  *buffer = static_cast<unsigned char>(a_amrCell.m_isGhostCF);
  buffer++;

  *buffer = static_cast<unsigned char>(a_amrCell.m_isDomainCell);
  buffer++;
}

#include <CD_NamespaceFooter.H>
