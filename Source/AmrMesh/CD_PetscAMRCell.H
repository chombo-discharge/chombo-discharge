/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file CD_PetscAMRCell.H
  @brief Declaration of a cell type for indicating if a cell is covered, ghost, etc.
  @author Robert Marskar
*/

#ifndef CD_PetscAMRCell_H
#define CD_PetscAMRCell_H

#ifdef CH_USE_PETSC

// Petsc includes
#include "petsc.h"

// Chombo includes
#include <SPMD.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for storing grid information and corresponding PETSc rows. Used for mapping
  between Chombo and PETSC

  @details This class is not very sophisticated. It is meant to be plugged into a BaseFab<PetscAMRCell>, and it will store the cell information
  of the corresponding cell (in an AMR context). Various cell flags are available, such as whether or not the cell is a ghost cell, if it is covered
  by another refined cell. It also stores the corresponding PETSc rows in the cell. Note that we permit two PETSc rows per cell since we're also
  dealing with jump conditions where we need to store the cell-centered potential on two phases. 
*/
class PetscAMRCell
{
public:
  /*!
    @brief Constructor - does nothing
  */
  PetscAMRCell() noexcept;

  /*!
    @brief Destructor - does nothing
  */
  ~PetscAMRCell() noexcept;

  /*!
    @brief Copy constructor -- can use default
  */
  PetscAMRCell(const PetscAMRCell&) = default;

  /*!
    @brief Move constructor -- can use default
  */
  PetscAMRCell(PetscAMRCell&&) = default;

  /*!
    @brief Copy assignment, can use default
  */
  PetscAMRCell&
  operator=(const PetscAMRCell&) = default;

  /*!
    @brief Move assignment, can use default
  */
  PetscAMRCell&
  operator=(PetscAMRCell&&) = default;

  /*!
    @brief Define operator. Sets all parameters.
    @param[in] a_petscRowPhase0 PETSc row corresponding to this grid cell on phase 0 (gas phase)
    @param[in] a_petscRowPhase1 PETSc row corresponding to this grid cell on phase 1 (solid phase)
    @param[in] a_isCoveredByFinerGrid Covered by finer grid or not.
    @param[in] a_isCoarCF If true, the cell lies on the coarse side of the refinement boundary
    @param[in] a_isFineCF If true, the cell lies on the fine side of the refinement boundary    
    @param[in] a_isGhostCF Cell lies outside the grid level.
    @param[in] a_isDomainBoundaryCell Cell abuts the domain boundary.
  */
  void
  define(const PetscInt a_petscRowPhase0,
         const PetscInt a_petscRowPhase1,
         const bool     a_isCoveredByFinerGrid,
         const bool     a_isCoarCF,
         const bool     a_isFineCF,
         const bool     a_isGhostCF,
         const bool     a_isDomainBoundaryCell) noexcept;

  /*!
    @brief Set corresponding row in PETSc for specific phase.
    @param[in] a_phase Phase
    @param[in] a_row PETSc row
  */
  void
  setPetscRow(const int a_phase, const PetscInt a_row) noexcept;

  /*!
    @brief Set covered by finer grid flag
    @param[in] a_coveredByFinerGrid If true, cell is covered by a finer grid.
  */
  void
  setCoveredByFinerGrid(const bool a_coveredByFinerGrid) noexcept;

  /*!
    @brief Set coarse-side CF flag
    @param[in] a_isCoarCF If true, cell abuts the refinement boundary from the coarse side
  */
  void
  setCoarCF(const bool a_isCoarCF) noexcept;

  /*!
    @brief Set fine-side CF flag
    @param[in] a_isFineCF If true, cell abuts the refinement boundary from the fine side
  */
  void
  setFineCF(const bool a_isFineCF) noexcept;

  /*!
    @brief Set ghosted coarse-fined flag
    @param[in] a_isGhostCF If true, cell is a ghost cell outside the grid level
  */
  void
  setGhostCF(const bool a_isGhostCF) noexcept;

  /*!
    @brief Set domain boundary flag
    @param[in] a_isDomainBoundaryCell If true, the cell abuts the domain boundary.
  */
  void
  setDomainBoundaryCell(const bool a_isDomainBoundaryCell) noexcept;

  /*!
    @brief Get the PETSc row corresponding to this cell and phase
  */
  PetscInt
  getPetscRow(const int a_iphase) const noexcept;

  /*!
    @brief Check if the cell is covered by a finer grid
  */
  bool
  isCoveredByFinerGrid() const noexcept;

  /*!
    @brief Check if the cell lies on the coarse-side of the refinement interface
  */
  bool
  isCoarCF() const noexcept;

  /*!
    @brief Check if the cell lies on the fine-side of the refinement interface
  */
  bool
  isFineCF() const noexcept;

  /*!
    @brief Check if the cell is a ghost cell outside this grid level
  */
  bool
  isGhostCF() const noexcept;

  /*!
    @brief Check if the cell is a domain cell
  */
  bool
  isDomainBoundaryCell() const noexcept;

  /*!
    @brief If >= 0, this contains the PETSc row corresponding to the grid cell on phase 0
  */
  PetscInt m_petscRows[2];

  /*!
    @brief If true, the cell is covered by a finer grid
  */
  bool m_isCoveredByFinerGrid;

  /*!
    @brief If true, the cell lies outside the grid level.
  */
  bool m_isGhostCF;

  /*!
    @brief If true, the cell lies on the coarse side of the coarse-fine boundary
  */
  bool m_isCoarCF;

  /*!
    @brief If true, the cell lies on the fine side of the coarse-fine boundary
  */
  bool m_isFineCF;

  /*!
    @brief If true, the cell is a domain boundary cell
   */
  bool m_isDomainBoundaryCell;
};

/*!
  @brief PetscAMRCell specialization of the linearization size function used by LevelData
*/
constexpr int
linearSize(const PetscAMRCell& a_amrCell);

/*!
  @brief PetscAMRCell specialization of the linearization function used by LevelData
  @param[out] a_amrCell Input object
  @param[in] a_buffer Buffer
*/
void
linearIn(PetscAMRCell& a_amrCell, const void* const a_buffer);

/*!
  @brief PetscAMRCell specialization of the linearization function used by LevelData
  @param[out] a_buffer Buffer to which we write the data
  @param[in] a_amrCell Object to linearize
*/
void
linearOut(void* const a_buffer, const PetscAMRCell& a_amrCell);

#include <CD_NamespaceFooter.H>

#endif

#endif
