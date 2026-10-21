/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscGrid.H
 * @brief  File containing declaration of PetscGrid.
 * @author Robert Marskar
 */

#ifndef CD_PETSCGRID_H
#define CD_PETSCGRID_H

#ifdef CH_USE_PETSC

// Petsc includes
#include "petsc.h"
#include "petscis.h"

// Chombo includes
#include <DisjointBoxLayout.H>
#include <EBCellFAB.H>
#include <MFCellFAB.H>
#include <EBISLayout.H>
#include <VoFIterator.H>
#include <Stencils.H>

// Our includes
#include <CD_DischargeIO.H>
#include <CD_EBAMRData.H>
#include <CD_MFLevelGrid.H>
#include <CD_PetscDOF.H>
#include <CD_PetscAMRCell.H>
#include <CD_PetscStencil.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Class which maps a multifluid cut-cell AMR hierarchy to PETSc.
 * @details This class includes functions for putting chombo-discharge in PETSc and vice versa.
 * @note This class assumes that PETSc will only care about a single component per grid cell (two components in a
 * multifluid cell).
 */
class PetscGrid
{
public:
  /**
   * @brief Default constructor. Must subsequently call define
   */
  PetscGrid() noexcept;

  /**
   * @brief Destructor. Does nothing.
   */
  virtual ~PetscGrid() noexcept;

  /**
   * @brief Main define function. Puts object in usable state.
   * @details Note that a_levelGridsCoFi must be coarsened from finer grids. That is, entry lvl
   * in levelGridsCoFi must be a coarsened version of a_levelGrids[lvl+1]. Similarly, a_levelGridsFiCo
   * must be a refinement, i.e., a_levelGridsFiCo[lvl] must be a refinemetn of a_levelGrids[lvl].
   * @param[in] a_levelGrids chombo-discharge AMR grids
   * @param[in] a_levelGridsCoFi Coarsened grids.
   * @param[in] a_levelGridsFiCo Refined grids.
   * @param[in] a_validCells Data holder which is true for cells not covered by a finer grid and false otherwise.
   * @param[in] a_coarHaloCF Data holder which is true on cells that lie on the coarse side of a CF interface
   * @param[in] a_fineHaloCF Data holder which is true on cells that lie on the fine side of a CF interface
   * @param[in] a_refinementRatios Refinement ratios between grid level.
   * @param[in] a_dx Grid resolutions on each level.
   * @param[in] a_finestLevel Finest defined grid level.
   * @param[in] a_numGhost Number of ghost cells described in the grid.
   */
  virtual void
  define(const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGrids,
         const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGridsCoFi,
         const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGridsFiCo,
         const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_validCells,
         const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_coarHaloCF,
         const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_fineHaloCF,
         const Vector<int>&                                     a_refinementRatios,
         const Vector<Real>&                                    a_dx,
         const int                                              a_finestLevel,
         const int                                              a_numGhost) noexcept;

  /**
   * @brief Clear the object.
   * @details This deallocates all PETSc-relevant allocations and marks the grid as undefined, so calling it twice is
   * safe.
   */
  virtual void
  clear() noexcept;

  /**
   * @brief Is defined or not
   * @return True if the grid has been defined.
   */
  virtual bool
  isDefined() const noexcept;

  /**
   * @brief Get the number of phases.
   * @return Number of phases (1 or 2).
   */
  virtual int
  getNumPhases() const noexcept;

  /**
   * @brief Get the grids on each level.
   * @return m_levelGrids
   */
  virtual const Vector<RefCountedPtr<MFLevelGrid>>&
  getMFLevelGrids() const noexcept;

  /**
   * @brief Get the coarsened fine grids.
   * @details Entry lvl is level lvl+1 coarsened to the resolution of level lvl, and shares DataIndex with level lvl+1.
   * The finest level has no entry.
   * @return m_levelGridsCoFi
   */
  virtual const Vector<RefCountedPtr<MFLevelGrid>>&
  getMFLevelGridsCoFi() const noexcept;

  /**
   * @brief Get the refined coarse grids.
   * @details Entry lvl is level lvl-1 refined to the resolution of level lvl, and shares DataIndex with level lvl-1.
   * Level 0 has no entry.
   * @return m_levelGridsFiCo
   */
  virtual const Vector<RefCountedPtr<MFLevelGrid>>&
  getMFLevelGridsFiCo() const noexcept;

  /**
   * @brief Get the map from AMR cells to PETSc rows, on the grids of each level.
   * @details Valid cells and ghost cells both hold global PETSc rows.
   * @return m_amrToPetsc
   */
  virtual const Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>>&
  getAMRToPetsc() const noexcept;

  /**
   * @brief Get the map from AMR cells to PETSc rows, on the coarsened fine grids.
   * @details Entry lvl holds the rows of level lvl on the grids of getMFLevelGridsCoFi()[lvl].
   * @return m_amrToPetscCoFi
   */
  virtual const Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>>&
  getAMRToPetscCoFi() const noexcept;

  /**
   * @brief Get the map from AMR cells to PETSc rows, on the refined coarse grids.
   * @details Entry lvl holds the rows of level lvl on the grids of getMFLevelGridsFiCo()[lvl].
   * @return m_amrToPetscFiCo
   */
  virtual const Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>>&
  getAMRToPetscFiCo() const noexcept;

  /**
   * @brief Get iterators over the cut cells that are PETSc rows on a phase.
   * @details Non-const because a VoFIterator must be reset to be iterated over.
   * @param[in] a_iphase Phase index
   * @return m_vofIterators[a_iphase], one LayoutData per level.
   */
  virtual Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
  getVoFIterator(const int a_iphase) const noexcept;

  /**
   * @brief Rewrite a stencil over cells on one level as a stencil over global PETSc rows.
   * @details Each entry of a_stencil is a cell on level a_level, in the index space of grid a_din, possibly in its
   * ghost region. Every entry must be a PETSc row on a_phase, i.e. a valid cell of a_din or of a neighbouring box.
   * Anything else aborts with the reason: a cell covered by a finer level, a coarse-fine ghost cell, a cell covered by
   * the EB, or a cell outside the domain. Ghost cells and boundary conditions are the caller's to resolve first.
   * @param[in] a_stencil Stencil over cells on level a_level.
   * @param[in] a_level   Grid level.
   * @param[in] a_din     Grid index on level a_level.
   * @param[in] a_phase   Phase.
   * @return The same stencil, over global PETSc rows, compressed (unique, ascending rows).
   */
  virtual PetscStencil
  getPetscStencil(const VoFStencil& a_stencil,
                  const int         a_level,
                  const DataIndex&  a_din,
                  const int         a_phase) const noexcept;

  /**
   * @brief Rewrite a two-level stencil, over cells on a level and the next finer level, as a stencil over global PETSc
   * rows.
   * @details This is the form produced by LeastSquares::computeDualLevelStencils for a stencil on level a_level that
   * uses level a_level+1 cells. a_stencil resolves as in the single-level version. a_fineStencil is over cells on level
   * a_level+1, in the index space of grid a_din refined to level a_level+1 (getMFLevelGridsFiCo), and each of its
   * entries must be a PETSc row on level a_level+1.
   * @param[in] a_stencil     Stencil over cells on level a_level.
   * @param[in] a_fineStencil Stencil over cells on level a_level+1.
   * @param[in] a_level       Grid level of a_stencil. Must be below the finest level.
   * @param[in] a_din         Grid index on level a_level.
   * @param[in] a_phase       Phase.
   * @return The sum of both stencils, over global PETSc rows, compressed (unique, ascending rows).
   */
  virtual PetscStencil
  getPetscStencilWithFinerLevel(const VoFStencil& a_stencil,
                                const VoFStencil& a_fineStencil,
                                const int         a_level,
                                const DataIndex&  a_din,
                                const int         a_phase) const noexcept;

  /**
   * @brief Rewrite a two-level stencil, over cells on a level and the next coarser level, as a stencil over global
   * PETSc rows.
   * @details This is the form of the coarse-fine ghost cell stencils of EBMultigridInterpolator::getGhostStencil.
   * a_stencil resolves as in the single-level version. a_coarStencil is over cells on level a_level-1, in the index
   * space of the level a_level grid a_din coarsened to level a_level-1 (getMFLevelGridsCoFi), and each of its entries
   * must be a PETSc row on level a_level-1.
   * @param[in] a_stencil     Stencil over cells on level a_level.
   * @param[in] a_coarStencil Stencil over cells on level a_level-1.
   * @param[in] a_level       Grid level of a_stencil. Must be above level 0.
   * @param[in] a_din         Grid index on level a_level.
   * @param[in] a_phase       Phase.
   * @return The sum of both stencils, over global PETSc rows, compressed (unique, ascending rows).
   */
  virtual PetscStencil
  getPetscStencilWithCoarserLevel(const VoFStencil& a_stencil,
                                  const VoFStencil& a_coarStencil,
                                  const int         a_level,
                                  const DataIndex&  a_din,
                                  const int         a_phase) const noexcept;

  /**
   * @brief Create a PETSc vector.
   * @details This will allocate the PETSc vector. It is up to the user to destroy it when it is no longer needed.
   * @param[in,out] a_x PETSc vector
   */
  virtual void
  create(Vec& a_x) noexcept;

  /**
   * @brief Destroy a PETSc vector.
   * @details This will deallocate the PETSc vector.
   * @param[in,out] a_x PETSc vector
   */
  virtual void
  destroy(Vec& a_x) noexcept;

  /**
   * @brief Set values in a PETSc vector.
   * @details a_value Does not need to be the same across all ranks.
   * @param[in,out] a_x     Global PETSc vector
   * @param[in]     a_value Value to set
   */
  virtual void
  setValue(Vec& a_x, const PetscScalar a_value) const noexcept;

  /**
   * @brief Push data form a chombo-discharge data holder into PETSc
   * @param[out] a_x PETSc vector
   * @param[in] a_y chombo-discharge data holder
   */
  virtual void
  putChomboInPetsc(Vec& a_x, const MFAMRCellData& a_y) const noexcept;

  /**
   * @brief Pull data from a PETSc vector into a chombo-discharge data holder
   * @param[in,out] a_y chombo-discharge data holder
   * @param[in] a_x PETSc vector
   */
  virtual void
  putPetscInChombo(MFAMRCellData& a_y, const Vec& a_x) const noexcept;

  /**
   * @brief Write the PETSc grid to an output HDF5 file
   * @param[in] a_filename Output file name
   */
  virtual void
  dumpPetscGrid(const std::string& a_filename) const noexcept;

protected:
  /**
   * @brief Is defined or not
   */
  bool m_isDefined;

  /**
   * @brief Verbose or not
   */
  bool m_verbose;

  /**
   * @brief Run-time debug or not
   */
  bool m_debug;

  /**
   * @brief Run-time profile or not
   */
  bool m_profile;

  /**
   * @brief Finest grid level
   */
  int m_finestLevel;

  /**
   * @brief Number of ghost cells
   */
  int m_numGhost;

  /**
   * @brief Number of phases
   */
  int m_numPhases;

  /**
   * @brief Number of rows owned locally
   */
  PetscInt m_numLocalRows;

  /**
   * @brief Number of global rows
   */
  PetscInt m_numGlobalRows;

  /**
   * @brief Global row corresponding to the first local row
   */
  PetscInt m_localRowBegin;

  /**
   * @brief Local index to global index mapping for PETSc
   */
  ISLocalToGlobalMapping m_localToGlobalIS;

  /**
   * @brief Number of rows per MPI rank. Visible to all ranks.
   */
  Vector<PetscInt> m_numRowsPerRank;

  /**
   * @brief Grids on each level
   */
  Vector<RefCountedPtr<MFLevelGrid>> m_levelGrids;

  /**
   * @brief Coarsened grids
   */
  Vector<RefCountedPtr<MFLevelGrid>> m_levelGridsCoFi;

  /**
   * @brief Refined grids
   */
  Vector<RefCountedPtr<MFLevelGrid>> m_levelGridsFiCo;

  /**
   * @brief Data holder which is true for cells not covered by a finer grid.
   */
  Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> m_validCells;

  /**
   * @brief Coarse-side cells on refinement interface
   */
  Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> m_coarHaloCF;

  /**
   * @brief Fine-side cells on refinement interface
   */
  Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> m_fineHaloCF;

  /**
   * @brief Refinement ratios between levels.
   */
  Vector<int> m_refRat;

  /**
   * @brief Grid resolutions
   */
  Vector<Real> m_dx;

  /**
   * @brief Mapping structure for moving between AMR cells and PETSc rows.
   * @details PetscAMRCell contains extra information, such as whether or not a cell is a ghost cell, covered by another
   * cell, etc.
   * @note This holds global PETSc rows, also in ghost cells. Subtract m_localRowBegin to get the local row.
   */
  Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>> m_amrToPetsc;

  /**
   * @brief Mapping structure for moving between AMR cells and PETSc rows.
   * @details This is like m_amrToPetsc, except that it is defined on the coarsened grids. Index i thus contains the
   * coarsening of level i+1.
   * @note This holds global PETSc rows.
   */
  Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>> m_amrToPetscCoFi;

  /**
   * @brief Mapping structure for moving between AMR cells and PETSc rows.
   * @details This is like m_amrToPetsc, except that it is defined on the refined. Index i+1 thus contains the
   * refinement of level i.
   * @note This holds global PETSc rows.
   */
  Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>> m_amrToPetscFiCo;

  /**
   * @brief Mapping of PETSc local row index back to AMR
   * @details Locally owned. Maps a local row index to a level/grid/cell.
   */
  Vector<RefCountedPtr<LayoutData<Vector<PetscDOF>>>> m_petscToAMR[2];

  /**
   * @brief Iterators over the cut cells that are PETSc rows, per phase and level.
   * @details Mutable because iterating over a VoFIterator resets it.
   */
  mutable Vector<RefCountedPtr<LayoutData<VoFIterator>>> m_vofIterators[2];

  /**
   * @brief Define basic cell information for PetscGrid.
   * @details This defines basic cell properties, indicating whether or not a cell is covered by a finer grid, if it
   * abuts the coarse-fine interface, etc.
   */
  virtual void
  defineAMRCells() noexcept;

  /**
   * @brief Build the Chombo <--> Petsc index mapping for the valid cells.
   * @details This defines m_localRows and m_globalRows, m_numRowsPerRank, m_localToGlobalIS, and figures out the number
   * of degrees of freedom on each MPI rank. After this is called, it will be possible to index local grid cells into a
   * PETSc vector row.
   */
  virtual void
  definePetscRows() noexcept;

  /**
   * @brief Add the rows of a stencil's cells to a PETSc stencil.
   * @details Looks up each entry of a_stencil in a_map. Every entry must lie in a_map's box and be a PETSc row on
   * a_phase. Either failure aborts in every build, through abortInvalidEntry: an entry outside the map would read out
   * of bounds, and PETSc silently ignores negative indices.
   * @param[in,out] a_petscStencil PETSc stencil to add to.
   * @param[in]     a_stencil      Stencil over the cells of a_map.
   * @param[in]     a_map          Map from cells to PETSc rows, for one grid patch.
   * @param[in]     a_phase        Phase.
   * @param[in]     a_what         Which stencil this is, for the abort message.
   */
  virtual void
  addRows(PetscStencil&                a_petscStencil,
          const VoFStencil&            a_stencil,
          const BaseFab<PetscAMRCell>& a_map,
          const int                    a_phase,
          const std::string&           a_what) const noexcept;

  /**
   * @brief Abort because a stencil entry lies outside the map or is not a PETSc row, naming the cell and the reason.
   * @param[in] a_iv    Cell.
   * @param[in] a_map   Map from cells to PETSc rows, for one grid patch.
   * @param[in] a_phase Phase.
   * @param[in] a_what  Which stencil the entry belongs to.
   */
  virtual void
  abortInvalidEntry(const IntVect&               a_iv,
                    const BaseFab<PetscAMRCell>& a_map,
                    const int                    a_phase,
                    const std::string&           a_what) const noexcept;

  /**
   * @brief Define the coarsened/refined views of m_amrToPetsc.
   */
  virtual void
  defineCoFiBuffers() noexcept;
};

#include <CD_NamespaceFooter.H>

#endif

#endif
