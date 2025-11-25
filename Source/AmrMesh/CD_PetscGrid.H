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
    @brief A cell missing an identifier (usually signals a bug).
  */
  static const int InvalidCell;

  /*!
    @brief Genuine cell (i.e., cell not covered by geometry or a finer grid)
  */
  static const int GenuineCell;

  /*!
    @brief Ghost cell (i.e., a cell outside a fine-grid patch)
  */
  static const int GhostCell;

  /*!
    @brief Covered cell (i.e., a cell covered by a finer grid)
  */
  static const int CoveredCell;

  /*!
    @brief Exterior cell (i.e., a cell covered by the geometry, or outside the problem domain)
  */
  static const int ExteriorCell;    

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
    @param[in] a_levelGrids chombo-discharge AMR grids
    @param[in] a_levelGridsCoFi Refined grids
    @param[in] a_levelGridsFiCo Coarsened grids
    @param[in] a_validCells Valid grid cells
    @param[in] a_finestLevel Finest defined grid level.
    @param[in] a_numGhost Number of ghost cells described in the grid.
  */
  virtual void
  define(const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGrids,
         const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGridsCoFi,
         const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGridsFiCo,         
         const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_validCells,
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
    @brief Number of local DOFs
  */
  PetscInt m_numLocalDOFs;

  /*!
    @brief Number of global DOFs
  */
  PetscInt m_numGlobalDOFs;

  /*!
    @brief Global DOF corresponding to the first local DOF.
  */
  PetscInt m_localDOFBegin;

  /*!
    @brief Local index to global index mapping for PETSc
  */
  ISLocalToGlobalMapping m_localToGlobalIS;

  /*!
    @brief Number of DOFs per MPI rank. Visible to all ranks. 
  */
  Vector<PetscInt> m_numDOFsPerRank;

  /*!
    @brief Valid grid cells
  */
  Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> m_validCells;

  /*!
    @brief Encoded grid cells.
    @details These cells are encoded by whether or not they are a genuine cell, a ghost cell,
    covered by a finer grid, or outside the domain (completely inside EBs or outside the problem domain).
  */
  Vector<RefCountedPtr<LevelData<BaseFab<int>>>> m_cellTypes;

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
    @brief Mapping of AMR DOFs to local PETSc DOFs
    @details This contains, on each grid, the row number in the local PETSC vector corresponding to each valid DOF.
  */
  Vector<RefCountedPtr<LayoutData<BaseFab<PetscInt>>>> m_localIndices[2];
  
  /*!
    @brief Mapping of AMR DOFs to local PETSc DOFs
    @details This contains, on each grid, the row number in the global PETSc vector corresponding to each valid DOF.
  */
  Vector<RefCountedPtr<LayoutData<BaseFab<PetscInt>>>> m_globalIndices[2];  

  /*!
    @brief Mapping of PETSc local DOF index back to AMR
    @details Locally owned. Maps a local DOF index to a level/grid/cell DOF index.
  */
  Vector<PetscDOF> m_petscToAMR;

  /*!
    @brief Build the Chombo <--> Petsc index mapping for the valid cells.
    @details This defines m_localIndices and m_globalIndices, m_numDOFsPerRank, m_localToGlobalIS, and figures out the number
    of degrees of freedom on each MPI rank. After this is called, it will be possible to index local grid cells into a PETSc
    vector row.
  */
  virtual void 
  definePetscDOFs() noexcept;  

  /*!
    @brief Build the map over valid degrees of freedom, done so that each grid has a local view of the
    surrounding global DOFs. 
  */
  virtual void
  defineLocalCompositeView() noexcept;


};

#endif

#include <CD_NamespaceFooter.H>

#endif
