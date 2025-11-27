/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PetscGrid.H
  @brief  File containing declaration of PetscGrid.
  @author Robert Marskar
*/

#ifndef CD_PetscGrid_H
#define CD_PetscGrid_H

#ifdef CH_USE_PETSC

#include "petsc.h"
#include "petscis.h"

// Chombo includes
#include <DisjointBoxLayout.H>
#include <EBCellFAB.H>
#include <MFCellFAB.H>
#include <EBISLayout.H>

// Our includes
#include <CD_EBAMRData.H>
#include <CD_MFLevelGrid.H>
#include <CD_PetscDOF.H>
#include <CD_AMRCell.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which maps a multifluid cut-cell AMR hierarchy to PETSc.
  @details This class includes functions for putting chombo-discharge in PETSc and vice versa.
  @note This class assumes that PETSc will only care about a single component per grid cell (two components in a multifluid cell).
*/
class PetscGrid
{
public:
  /*!
    @brief Default constructor. Must subsequently call define
  */
  PetscGrid() noexcept;

  /*!
    @brief Destructor. Does nothing.
  */
  virtual ~PetscGrid() noexcept;

  /*!
    @brief Main define function. Puts object in usable state.
    @details Note that a_levelGridsCoFi must be coarsened from finer grids. That is, entry lvl
    in levelGridsCoFi must be a coarsened version of a_levelGrids[lvl+1]. Similarly, a_levelGridsFiCo
    must be a refinement, i.e., a_levelGridsFiCo[lvl] must be a refinemetn of a_levelGrids[lvl].
    @param[in] a_levelGrids chombo-discharge AMR grids
    @param[in] a_levelGridsCoFi Coarsened grids.
    @param[in] a_levelGridsFiCo Refined grids.
    @param[in] a_amrCells Grid cell types
    @param[in] a_finestLevel Finest defined grid level.
    @param[in] a_numGhost Number of ghost cells described in the grid.
  */
  virtual void
  define(const Vector<RefCountedPtr<MFLevelGrid>>&                 a_levelGrids,
         const Vector<RefCountedPtr<MFLevelGrid>>&                 a_levelGridsCoFi,
         const Vector<RefCountedPtr<MFLevelGrid>>&                 a_levelGridsFiCo,
         const Vector<RefCountedPtr<LevelData<BaseFab<AMRCell>>>>& a_amrCells,
         const int                                                 a_finestLevel,
         const int                                                 a_numGhost) noexcept;

  /*!
    @brief Clear the object.
    @details This deallocates all PETSc-relevant allocations.
  */
  virtual void
  clear() noexcept;

  /*!
    @brief Create a PETSc vector.
    @details This will allocate the PETSc vector. It is up to the user to destroy it when it is no longer needed.
    param[in,out] x PETSc vector    
  */
  virtual void
  create(Vec& x) noexcept;

  /*!
    @brief Destroy a PETSc vector.
    @details This will deallocate the PETSc vector.
    param[in,out] x PETSc vector
  */
  virtual void
  destroy(Vec& x) noexcept;

  /*!
    @brief Set values in a PETSc vector.
    @details a_value Does not need to be the same across all ranks. 
    @param[in,out] x Global PETSc vector
  */
  virtual void
  setValue(Vec& a_x, const PetscScalar a_value) const noexcept;

  /*!
    @brief Push data form a chombo-discharge data holder into PETSc
    @param[out] a_x PETSc vector
    @param[in] a_y chombo-discharge data holder    
  */
  virtual void
  putChomboInPetsc(Vec& a_x, const MFAMRCellData& a_y) const noexcept;

  /*!
    @brief Pull data from a PETSc vector into a chombo-discharge data holder
    @param[in,out] a_y chombo-discharge data holder
    @param[in] a_x PETSc vector    
  */
  virtual void
  putPetscInChombo(MFAMRCellData& a_y, const Vec& a_x) const noexcept;

protected:
  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Verbose or not
  */
  bool m_verbose;

  /*!
    @brief Run-time debug or not
  */
  bool m_debug;

  /*!
    @brief Run-time profile or not
  */
  bool m_profile;

  /*!
    @brief Finest grid level
  */
  int m_finestLevel;

  /*!
    @brief Number of ghost cells
  */
  int m_numGhost;

  /*!
    @brief Number of phases
  */
  int m_numPhases;

  /*!
    @brief Number of rows owned locally
  */
  PetscInt m_numLocalRows;

  /*!
    @brief Number of global rows
  */
  PetscInt m_numGlobalRows;

  /*!
    @brief Global row corresponding to the first local row
  */
  PetscInt m_localRowBegin;

  /*!
    @brief Local index to global index mapping for PETSc
  */
  ISLocalToGlobalMapping m_localToGlobalIS;

  /*!
    @brief Number of rows per MPI rank. Visible to all ranks. 
  */
  Vector<PetscInt> m_numRowsPerRank;

  /*!
    @brief Valid grid cells
  */
  Vector<RefCountedPtr<LevelData<BaseFab<AMRCell>>>> m_amrCells;

  /*!
    @brief Grids on each level
  */
  Vector<RefCountedPtr<MFLevelGrid>> m_levelGrids;

  /*!
    @brief Coarsened grids
  */
  Vector<RefCountedPtr<MFLevelGrid>> m_levelGridsCoFi;

  /*!
    @brief Refined grids
  */
  Vector<RefCountedPtr<MFLevelGrid>> m_levelGridsFiCo;

  /*!
    @brief Mapping of AMR cells to local PETSc rows
    @details This contains, on each grid, the row number in the local PETSC vector corresponding to each valid grid cell
  */
  Vector<RefCountedPtr<LevelData<BaseFab<PetscInt>>>> m_localRows[2];

  /*!
    @brief Mapping of AMR cells to global PETSc rows
    @details This contains, on each grid, the row number in the global PETSc vector corresponding to each valid grid cell
  */
  Vector<RefCountedPtr<LevelData<BaseFab<PetscInt>>>> m_globalRows[2];

  /*!
    @brief Rows in the PETSc vector corresponding to a particular grid cell or face. If the value is negative it implies
    that the cell is either a ghost, invalid, covered by a finer grid, etc. 
    @details This contains the global row # in the PETsc vector corresponding to a grid cell (similar to m_globalRows). The
    only difference is that this contain also has ghost cells so it can access the global row number for cells that lie outside
    of the domain. Similarly, negative values are assigned to indicate if cells are ghost cells, covered by another grid, covered
    or covered by the geometry.
  */
  Vector<RefCountedPtr<LevelData<BaseFab<PetscInt>>>> m_rows[2];

  /*!
    @brief Rows in the PETSc vector corresponding to a particular grid cell or face. If the value is negative it implies
    that the cell is either a ghost, invalid, covered by a finer grid, etc.     
    @details This data holder is filled with row indices, including ghost cells around each patch. The grids on which
    the objects are defined are coarsenings from finer levels. I.e., entry on level lvl corresponds to a coarsening of level
    lvl+1. This means that if we are visiting level lvl+1 on m_rows, we can access the coarse grid by visiting level
    lvl in this data holder.
  */
  Vector<RefCountedPtr<LevelData<BaseFab<PetscInt>>>> m_rowsCoFi[2];

  /*!
    @brief Rows in the PETSc vector corresponding to a particular grid cell or face. If the value is negative it implies
    that the cell is either a ghost, invalid, covered by a finer grid, etc.         
    @details This data holder is filled with row indices, including ghost cells around each patch. The grids on which
    the objects are defined are refinements of coarser levels. I.e., the entry on level lvl+1 corresponds to a refinement
    of level lvl. This means that if we are visiting level lvl on m_rows, we can access the fine grids via this data
    holder (i.e., level lvl+1 in this data holder).
  */
  Vector<RefCountedPtr<LevelData<BaseFab<PetscInt>>>> m_rowsFiCo[2];

  /*!
    @brief Map over which cells belong to ghosts, covered, invalid, boundary, etc.
    @details This is defined over m_levelGrids, which is the fundamental grid structure on which solvers will operatore. 
  */
  Vector<RefCountedPtr<LevelData<BaseFab<int>>>> m_cellFlags;

  /*!
    @brief Map over which cells belong to ghosts, covered, invalid, boundary, etc.
    @details This is defined over coarsened versions m_levelGrids. I.e., the entry on level lvl is on a grid that is a coarsening
    from lvl+1
  */
  Vector<RefCountedPtr<LevelData<BaseFab<int>>>> m_cellFlagsCoFi;

  /*!
    @brief Map over which cells belong to ghosts, covered, invalid, boundary, etc.
    @details This is defined over refine4d versions m_levelGrids. I.e., the entry on level lvl+1 is on a grid that is a refinement
    from lvl    
  */
  Vector<RefCountedPtr<LevelData<BaseFab<int>>>> m_cellFlagsFiCo;

  /*!
    @brief Mapping of PETSc local row index back to AMR
    @details Locally owned. Maps a local row index to a level/grid/cell.
  */
  Vector<PetscDOF> m_petscToAMR;

  /*!
    @brief Build the Chombo <--> Petsc index mapping for the valid cells.
    @details This defines m_localRows and m_globalRows, m_numRowsPerRank, m_localToGlobalIS, and figures out the number
    of degrees of freedom on each MPI rank. After this is called, it will be possible to index local grid cells into a PETSc
    vector row.
  */
  virtual void
  definePetscRows() noexcept;

  /*!
    @brief Build the cell flag mappings
  */
  virtual void
  defineCellFlags() noexcept;

  /*!
    @brief Build the view of the global PETSc rows around each grid patch.
  */
  virtual void
  defineRowViews() noexcept;
};

#endif

#include <CD_NamespaceFooter.H>

#endif
