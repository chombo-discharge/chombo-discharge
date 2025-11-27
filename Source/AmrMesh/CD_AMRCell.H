/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file CD_AMRCell.H
  @brief Declaration of a cell type for indicating if a cell is covered, ghost, etc. 
  @author Robert Marskar
*/

#ifndef CD_AMRCell_H
#define CD_AMRCell_H

// Chombo includes
#include <SPMD.H>

// Our includes
#include <CD_NamespaceHeader.H>

class AMRCell
{
public:
  /*!
    @brief Constructor - does nothing
  */
  AMRCell() noexcept;

  /*!
    @brief Destructor - does nothing
  */
  ~AMRCell() noexcept;

  /*!
    @brief Copy constructor -- can use default
  */
  AMRCell(const AMRCell&) = default;

  /*!
    @brief Move constructor -- can use default
  */
  AMRCell(AMRCell&&) = default;

  /*!
    @brief Copy assignment, can use default
  */
  AMRCell&
  operator=(const AMRCell&) = default;

  /*!
    @brief Move assignment, can use default
  */
  AMRCell&
  operator=(AMRCell&&) = default;

  /*!
    @brief Define operator. Sets all parameters.
    @param[in] a_numPhases Number of phases where this cell is a valid cell.
    @param[in] a_isCoveredByFinerGrid Covered by finer grid or not.
    @param[in] a_isGhostCF Cell lies outside the grid level.
    @param[in] a_isDomainCell Cell abuts the domain boundary.
  */
  void
  define(const unsigned int a_numPhases,
         const bool         a_isCoveredByFinerGrid,
         const bool         a_isGhostCF,
         const bool         a_isDomainCell) noexcept;

  /*!
    @brief Set number of phases
    @param[in] a_numPhases
  */
  void
  setNumPhases(const unsigned int a_numPhases) noexcept;

  /*!
    @brief Set covered by finer grid flag
    @param[in] a_coveredByFinerGrid If true, cell is covered by a finer grid.
  */
  void
  setCoveredByFinerGrid(const bool a_coveredByFinerGrid) noexcept;

  /*!
    @brief Set ghosted coarse-fined flag
    @param[in] a_isGhostCF If true, cell is a ghost cell outside the grid level
  */
  void
  setGhostCF(const bool a_isGhostCF) noexcept;

  /*!
    @brief Set domain boundary flag
    @param[in] a_isDomainCell If true, the cell abuts the domain boundary.
  */
  void
  setDomainCell(const bool a_isDomainCell) noexcept;

  /*!
    @brief Get number of phases.
  */
  unsigned int
  getNumPhases() const noexcept;

  /*!
    @brief Check if the cell is covered by a finer grid
  */
  bool
  isCoveredByFinerGrid() const noexcept;

  /*!
    @brief Check if the cell is a ghost cell outside this grid level
  */
  bool
  isGhostCF() const noexcept;

  /*!
    @brief Check if the cell is a domain cell
  */
  bool
  isDomainCell() const noexcept;

  /*!
    @brief If > 0, indicates the number of phases where this cell is a valid cell.
    @details If zero, it indicates the cell is an invalid cell (either entirely covered by both phases, or
    outside the domain.
  */
  unsigned char m_numPhases;

  /*!
    @brief If true, the cell is covered by a finer grid
  */
  bool m_isCoveredByFinerGrid;

  /*!
    @brief If > 0, the cell lies outside the grid level.
  */
  bool m_isGhostCF;

  /*!
    @brief If > 0, the cell is a domain boundary cell
   */
  bool m_isDomainCell;
};

/*!
  @brief AMRCell specialization of the linearization size function used by LevelData
*/
int
linearSize(const AMRCell& a_amrCell);

/*!
  @brief AMRCell specialization of the linearization function used by LevelData
  @param[out] a_amrCell Input object
  @param[in] a_buffer Buffer
*/
void
linearIn(AMRCell& a_amrCell, const void* const a_buffer);

/*!
  @brief AMRCell specialization of the linearization function used by LevelData
  @param[out] a_buffer Buffer to which we write the data
  @param[in] a_amrCell Object to linearize

*/
void
linearOut(void* const a_buffer, const AMRCell& a_amrCell);

#include <CD_NamespaceFooter.H>

#endif
