/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file CD_PetscDOF.H
  @brief Declaration of a simple mapping between a PETSc index and an AMR cell
  @author Robert Marskar
*/

#ifndef CD_PetscDOF_H
#define CD_PetscDOF_H

// Chombo includes
#include <IntVect.H>
#include <DataIndex.H>

// Our includes
#include <CD_NamespaceHeader.H>

struct PetscDOF
{
  PetscInt  localRow;
  IntVect   gridCell;
};

#include <CD_NamespaceFooter.H>

#endif
