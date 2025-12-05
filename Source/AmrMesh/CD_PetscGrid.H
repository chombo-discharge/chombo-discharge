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

// Petsc includes
#include "petsc.h"
#include "petscis.h"

// Chombo includes
#include <DisjointBoxLayout.H>
#include <EBCellFAB.H>
#include <MFCellFAB.H>
#include <EBISLayout.H>

// Our includes
#include <CD_DischargeIO.H>
#include <CD_EBAMRData.H>
#include <CD_MFLevelGrid.H>
#include <CD_PetscDOF.H>
#include <CD_PetscAMRCell.H>
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
    @param[in] a_validCells Data holder which is true for cells not covered by a finer grid and false otherwise.
    @param[in] a_coarHaloCF Data holder which is true on cells that lie on the coarse side of a CF interface
    @param[in] a_fineHaloCF Data holder which is true on cells that lie on the fine side of a CF interface    
    @param[in] a_finestLevel Finest defined grid level.
    @param[in] a_numGhost Number of ghost cells described in the grid.
  */
  virtual void
  define(const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGrids,
         const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGridsCoFi,
         const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGridsFiCo,
         const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_validCells,
         const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_coarHaloCF,
         const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_fineHaloCF,
         const int                                              a_finestLevel,
         const int                                              a_numGhost) noexcept;

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

  /*!
    @brief Write the PETSc grid to an output HDF5 file
  */
  virtual void
  dumpPetscGrid(const std::string a_filename) const noexcept;

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
    @brief Data holder which is true for cells not covered by a finer grid.
  */
  Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> m_validCells;

  /*!
    @brief Coarse-side cells on refinement interface
  */
  Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> m_coarHaloCF;

  /*!
    @brief Fine-side cells on refinement interface
  */
  Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> m_fineHaloCF;

  /*!
    @brief Mapping structure for moving between AMR cells and PETSc rows.
    @details PetscAMRCell contains extra information, such as whether or not a cell is a ghost cell, covered by another cell, etc. 
  */
  Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>> m_amrToPetsc;

  /*!
    @brief Mapping structure for moving between AMR cells and PETSc rows.
    @details This is like m_amrToPetsc, except that it is defined on the coarsened grids. Index i thus contains the coarsening of
    level i+1.
  */
  Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>> m_amrToPetscCoFi;

  /*!
    @brief Mapping structure for moving between AMR cells and PETSc rows.
    @details This is like m_amrToPetsc, except that it is defined on the refined. Index i+1 thus contains the refinement of
    level i.
  */
  Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>> m_amrToPetscFiCo;      

  /*!
    @brief Mapping of PETSc local row index back to AMR
    @details Locally owned. Maps a local row index to a level/grid/cell.
  */  
  Vector<RefCountedPtr<LayoutData<Vector<PetscDOF>>>> m_petscToAMR[2];

  /*!
    @brief Define basic cell information for PetscGrid.
    @details This defines basic cell properties, indicating whether or not a cell is covered by a finer grid, if it abuts
    the coarse-fine interface, etc.
  */
  virtual void
  defineAMRCells() noexcept;

  /*!
    @brief Build the Chombo <--> Petsc index mapping for the valid cells.
    @details This defines m_localRows and m_globalRows, m_numRowsPerRank, m_localToGlobalIS, and figures out the number
    of degrees of freedom on each MPI rank. After this is called, it will be possible to index local grid cells into a PETSc
    vector row.
  */
  virtual void
  definePetscRows() noexcept;

  /*!
    @brief Define the coarsened/refined views of m_amrToPetsc.
  */
  virtual void
  defineCoFiBuffers() noexcept;
};

#endif

#include <CD_NamespaceFooter.H>

#endif
