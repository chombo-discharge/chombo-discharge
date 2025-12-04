/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file CD_PetscAMRCell.cpp
  @brief Implementation of CD_PetscAMRCell.H
  @author Robert Marskar
*/

#ifdef CH_USE_PETSC

// Our includes
#include <CD_PetscAMRCell.H>
#include <CD_NamespaceHeader.H>

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

int
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
  const PetscInt* bufferInt = static_cast<const PetscInt*>(a_buffer);

  a_amrCell.m_petscRows[0] = *bufferInt;
  bufferInt++;

  a_amrCell.m_petscRows[1] = *bufferInt;
  bufferInt++;

  const unsigned char* bufferChar = static_cast<const unsigned char*>((void*)bufferInt);

  a_amrCell.m_isCoveredByFinerGrid = *bufferChar;
  bufferChar++;

  a_amrCell.m_isGhostCF = *bufferChar;
  bufferChar++;

  a_amrCell.m_isCoarCF = *bufferChar;
  bufferChar++;

  a_amrCell.m_isFineCF = *bufferChar;
  bufferChar++;

  a_amrCell.m_isDomainBoundaryCell = *bufferChar;
  bufferChar++;
}

void
linearOut(void* const a_buffer, const PetscAMRCell& a_amrCell)
{
  PetscInt* bufferInt = static_cast<PetscInt*>(a_buffer);

  *bufferInt = a_amrCell.m_petscRows[0];
  bufferInt++;

  *bufferInt = a_amrCell.m_petscRows[1];
  bufferInt++;

  unsigned char* bufferChar = static_cast<unsigned char*>((void*)bufferInt);

  *bufferChar = static_cast<unsigned char>(a_amrCell.m_isCoveredByFinerGrid);
  bufferChar++;

  *bufferChar = static_cast<unsigned char>(a_amrCell.m_isGhostCF);
  bufferChar++;

  *bufferChar = static_cast<unsigned char>(a_amrCell.m_isCoarCF);
  bufferChar++;

  *bufferChar = static_cast<unsigned char>(a_amrCell.m_isFineCF);
  bufferChar++;

  *bufferChar = static_cast<unsigned char>(a_amrCell.m_isDomainBoundaryCell);
  bufferChar++;
}

#include <CD_NamespaceFooter.H>

#endif
