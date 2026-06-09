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

struct PetscDOF
{
  PetscInt localRow;
  IntVect  gridCell;
};

#include <CD_NamespaceFooter.H>

#endif
