/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_EBCoarseFineParticleMesh.H
 * @brief  Declaration of a class for handling particle deposition arithmetic near refinement boundaries.
 * @author Robert Marskar
 */

#ifndef CD_EBCOARSEFINEPARTICLEMESH_H
#define CD_EBCOARSEFINEPARTICLEMESH_H

// Chombo includes
#include <DisjointBoxLayout.H>
#include <EBLevelGrid.H>
#include <EBCellFAB.H>
#include <ProblemDomain.H>

// Our includes
#include <CD_EBCoarAve.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief A class for moving deposited particle mass between two levels across a refinement boundary. This class does
 * not include any deposition functionality.
 * @details Both operations go through the coarsened fine grids:
 * 1. addFineGhostsToCoarse adds the part of the fine level's particle clouds that hangs over the refinement boundary
 *    (the fine ghost cells) to the coarse level.
 * 2. addInvalidCoarseToFine interpolates the coarse-level deposit under the fine grid to the fine level, for when the
 *    coarse-grid clouds deposit under the fine level but the mass should end up on the fine level.
 */
class EBCoarseFineParticleMesh
{
public:
  /**
   * @brief Default constructor. Leaves object in undefined state and subsequently requires a call to define.
   */
  EBCoarseFineParticleMesh() noexcept;

  /**
   * @brief Full constructor. Calls define.
   * @param[in] a_eblgCoar   Coarse grids
   * @param[in] a_eblgFine   Fine grids
   * @param[in] a_refRat     Refinement factory between levels
   * @param[in] a_ghost      Number of ghost cells in mesh data holders.
   */
  EBCoarseFineParticleMesh(const EBLevelGrid& a_eblgCoar,
                           const EBLevelGrid& a_eblgFine,
                           const int          a_refRat,
                           const IntVect      a_ghost) noexcept;

  /**
   * @brief Destructor (does nothing)
   */
  virtual ~EBCoarseFineParticleMesh() noexcept;

  /**
   * @brief Define function
   * @param[in] a_eblgCoar   Coarse grids
   * @param[in] a_eblgFine   Fine grids
   * @param[in] a_refRat     Refinement factory between levels
   * @param[in] a_ghost      Number of ghost cells in mesh data holders.
   */
  virtual void
  define(const EBLevelGrid& a_eblgCoar,
         const EBLevelGrid& a_eblgFine,
         const int          a_refRat,
         const IntVect      a_ghost) noexcept;

  /**
   * @brief Coarsen the fine-level ghost cells and add them to the coarse level.
   * @details On output, the contents in the ghost cells on the refinement boundary are added from the fine level to the
   * coarse level. This routine is used when the fine-level particles have particle clouds that deposit into the ghost
   * cells outside the fine level, i.e. over the refinement boundary. This routine coarsens that data and adds the mass
   * to the coarse level. The data is mass per cell volume, so each coarse VoF gets the sum of the fine ghost VoFs that
   * coarsen to it, divided by refRat^SpaceDim, which conserves mass.
   * @param[in,out] a_coarData Coarse mesh data.
   * @param[in]    a_fineData Fine mesh data.
   */
  virtual void
  addFineGhostsToCoarse(LevelData<EBCellFAB>& a_coarData, const LevelData<EBCellFAB>& a_fineData) const noexcept;

  /**
   * @brief Function which piecewise interpolates the data underneath the fine grid and adds it to the fine grid.
   * @details Every valid fine VoF gets the density of the coarse VoF it coarsens to: the data is mass per cell volume,
   * so in cut cells the coarse value is scaled by the ratio of the fine to the coarse volume fraction, which conserves
   * mass. The coarse valid data is copied onto the coarsened fine grids first, so no geometry finer than the fine
   * level's is needed.
   * @param[in,out] a_fineData Fine-level data.
   * @param[in]    a_coarData Coarse-level data.
   */
  virtual void
  addInvalidCoarseToFine(LevelData<EBCellFAB>& a_fineData, const LevelData<EBCellFAB>& a_coarData) const noexcept;

protected:
  /**
   * @brief Default component
   */
  static constexpr int m_comp = 0;

  /**
   * @brief Maximum number of components
   */
  static constexpr int m_nComp = 1;

  /**
   * @brief Is defined or not
   */
  bool m_isDefined;

  /**
   * @brief Verbose or not. Useful for debugging
   */
  bool m_verbose;

  /**
   * @brief Coarse EBLevelGrid
   */
  EBLevelGrid m_eblgCoar;

  /**
   * @brief Fine EBLevelGrid
   */
  EBLevelGrid m_eblgFine;

  /**
   * @brief Coarsening of m_eblgFine
   */
  EBLevelGrid m_eblgCoFi;

  /**
   * @brief Copier for copying data from the coarse level to the coarsened fine level.
   * @note This copier moves from valid in the source to valid in the destination.
   */
  Copier m_copierCoarToCoFi;

  /**
   * @brief Copier for copying data from the coarsened fine level to the coarse level.
   * @note This copier moves from valid+ghost in the source to valid+ghost in the destination.
   */
  Copier m_copierCoFiToCoarIncludeGhosts;

  /**
   * @brief VoFIterator over the fine ghost cells in the coarse cells of m_vofIterCoFiGhosts.
   */
  mutable LayoutData<VoFIterator> m_vofIterFineGhosts;

  /**
   * @brief VoFIterator over the coarse cells, on the coarsened fine grids, that are multi-valued or contain a
   * multi-valued fine ghost cell.
   * @details addFineGhostsToCoarse sums these per VoF; every other coarse cell takes the regular-cell sum.
   */
  mutable LayoutData<VoFIterator> m_vofIterCoFiGhosts;

  /**
   * @brief VoFIterator over the valid fine cells that lie in cut coarse cells.
   * @details addInvalidCoarseToFine interpolates these per VoF; every other fine cell takes the regular-cell value.
   */
  mutable LayoutData<VoFIterator> m_vofIterFineCut;

  /**
   * @brief Refinement ratio between the two levels
   */
  int m_refRat;

  /**
   * @brief Number of ghost cells in data holders.
   */
  IntVect m_ghost;

  /**
   * @brief Define the vof iterators
   */
  void
  defineVoFIterators() noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
