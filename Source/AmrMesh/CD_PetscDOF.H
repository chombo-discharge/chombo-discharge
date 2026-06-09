/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file CD_PetscDOF.H
  @brief Declaration of a simple mapping between a PETSc index and an AMR cell
  @author Robert Marskar
*/

#ifndef CD_PETSCDOF_H
#define CD_PETSCDOF_H

// Chombo includes
#include <IntVect.H>
#include <DataIndex.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Mapping between a PETSc local row index and an AMR grid cell.
*/
struct PetscDOF
{
  /** @brief Local PETSc row index. */
  PetscInt localRow;

  /** @brief AMR grid cell associated with this DOF. */
  IntVect gridCell;
};

#include <CD_NamespaceFooter.H>

#endif
