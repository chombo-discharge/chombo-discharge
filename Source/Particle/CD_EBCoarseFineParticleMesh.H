/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBCoarseFineParticleMesh.H
  @brief  Declaration of a class for handling particle deposition arithmetics near refinement boundaries.
  @author Robert Marskar
*/

#ifndef CD_EBCoarseFineParticleMesh_H
#define CD_EBCoarseFineParticleMesh_H

// Chombo includes
#include <DisjointBoxLayout.H>
#include <EBLevelGrid.H>
#include <EBCellFAB.H>
#include <ProblemDomain.H>

// Our includes
#include <CD_EBCoarAve.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief A class for handling particle deposition clouds that hang over refinement boundaries. This class does not include any deposition functionality but provides
  buffers for moving mass between two levels. 
  @details To use this class, the user can directly some functions:
     1. addFineGhostsToCoarse which adds the contribution of the fine level's particles to the coarse level. This
        will add the part of the particle clouds that hang over the refinement boundary to the coarse level. Likewise, the user can get a buffer (getFiCoBuffer) for
        getting a temporary data holder on which he can use special deposition procedures.
     2. addFiCoDataToFine where the user will have deposited coarse-level particles on a fine-grid buffer, and he/she wants to add the result back to the fine grid.
     3. addInvalidCoarseToFine where the coarse-grid deposition clouds are interpolated to the fine grid. This is the case when e.g. the coarse-grid clouds
        deposit underneath the fine level but the mass should end up on the fine level.
     4. exchangeAndAddFiCoData where deposited particles on the refined coarse grids have their masses synchronized (ghost cells are added across patches).
*/
class EBCoarseFineParticleMesh
{
public:
  /*!
    @brief Simple enum for switch between arithmetic and conservative averaging
  */
  enum Average
  {
    Arithmetic,
    Conservative
  };

  /*!
    @brief Default constructor. Leaves object in undefined state and subsequently requires a call to define. 
  */
  EBCoarseFineParticleMesh() noexcept;

  /*!
    @brief Full constructor. Calls define. 
    @param[in] a_eblgCoar   Coarse grids
    @param[in] a_eblgFine   Fine grids
    @param[in] a_refRat     Refinement factory between levels
    @param[in] a_ghost      Number of ghost cells in mesh data holders. 
  */
  EBCoarseFineParticleMesh(const EBLevelGrid& a_eblgCoar,
                           const EBLevelGrid& a_eblgFine,
                           const int          a_refRat,
                           const IntVect      a_ghost) noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBCoarseFineParticleMesh() noexcept;

  /*!
    @brief Define function
    @param[in] a_eblgCoar   Coarse grids
    @param[in] a_eblgFine   Fine grids
    @param[in] a_refRat     Refinement factory between levels
    @param[in] a_ghost      Number of ghost cells in mesh data holders. 
  */
  virtual void
  define(const EBLevelGrid& a_eblgCoar,
         const EBLevelGrid& a_eblgFine,
         const int          a_refRat,
         const IntVect      a_ghost) noexcept;

  /*!
    @brief Coarsen the fine-level ghost cells and add them to the coarse level. 
    @details On output, the contents in the ghost cells on the refinement boundary are added from the fine level to the coarse level. This routine
    is used when the fine-level particles have particle clouds that deposit into the ghost cells outside the fine level, i.e. over the refinement
    boundary. This routine coarsens that data and adds the mass to the coarse level. 
    @param[inout] a_coarData Coarse mesh data. 
    @param[in]    a_fineData Fine mesh data. 
  */
  virtual void
  addFineGhostsToCoarse(LevelData<EBCellFAB>& a_coarData, const LevelData<EBCellFAB>& a_fineData) const noexcept;

  /*!
    @brief Add the valid+ghost cells on the refined coarse level to the valid cells on the fine level.
    @details This routine is used when the user have deposited coarse-level particles into a_fiCoData (which is a refinement of the coarse level). This 
    routine then adds that data to the fine level. The motion plan defines the computation region as valid+ghost -> valid, which is what we want.
    @param[inout] a_fineData Fine mesh data
    @param[in]    a_fiCoData Refined coarse data buffer. Users will have deposited the mesh data onto this 
  */
  virtual void
  addFiCoDataToFine(LevelData<EBCellFAB>& a_fineData, const LevelData<EBCellFAB>& a_fiCoData) const noexcept;

  /*!
    @brief Function which piecewise interpolates the data underneath the fine grid and adds it to the fine grid. 
    @details The data-motion plan is to use the valid+ghost cells on the coarse level and add it into the valid region on the fine level. However, the 
    buffers' ghost regions are set to zero so this should only add non-ghosted data between the levels.
    @param[inout] a_fineData Fine-level data.
    @param[in]    a_coarData Coarse-level data.
  */
  virtual void
  addInvalidCoarseToFine(LevelData<EBCellFAB>& a_fineData, const LevelData<EBCellFAB>& a_coarData) const noexcept;

  /*!
    @brief Function which performs an exchange+add operation on the refined coarse data
  */
  virtual void
  exchangeAndAddFiCoData(LevelData<EBCellFAB>& a_fiCoData) const noexcept;

  /*!
    @brief Coarsen the data on the refined coarse grid to the coarse grid. The grids must be defined over the same layout
    @param[inout] a_coarData Coarse-grid data
    @param[in] a_fineData Data on the refined coarse grid
    @param[in] a_average Averaging method. Either arithmetic or conservative.
  */
  virtual void
  restrictAndAddFiCoDataToCoar(LevelData<EBCellFAB>&                   a_coarData,
                               const LevelData<EBCellFAB>&             a_fiCoData,
                               const EBCoarseFineParticleMesh::Average a_average) const noexcept;

  /*!
    @brief Get the refined version of the coarse eblevelgrid
  */
  const EBLevelGrid&
  getEblgFiCo() const;

protected:
  /*!
    @brief Default component
  */
  static constexpr int m_comp = 0;

  /*!
    @brief Maximum number of components
  */
  static constexpr int m_nComp = 1;

  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Verbose or not. Useful for debugging
  */
  bool m_verbose;

  /*!
    @brief Coarse EBLevelGrid
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Fine EBLevelGrid
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Refined coarse EBLevelGrid
  */
  EBLevelGrid m_eblgFiCo;

  /*!
    @brief Coarsening of m_eblgFine
  */
  EBLevelGrid m_eblgCoFi;

  /*!
    @brief Copier for copying data from the refined coarse level to the fine level.
    @note This copier moves from valid+ghost in the source to valid+ghost in the destination.
  */
  Copier m_copierFiCoToFineIncludeGhosts;

  /*!
    @brief Copier for copying data from the coarsened fine level to the coarse level.
    @note This copier moves from valid+ghost in the source to valid+ghost in the destination.
  */
  Copier m_copierCoFiToCoarIncludeGhosts;

  /*!
    @brief Copier for copying data from the refined coarse level to the fine level.
    @details This copier moves from valid+ghost in the source to valid in the destination (as opposed to m_copierCOFiToCoarIncludeGhosts which adds mass into the ghost cells).
  */
  Copier m_copierFiCoToFineNoGhosts;

  /*!
    @brief Copier for that moves data from valid+ghost -> valid on the refined coarse level.
    @details Needed for some of the deposition methods that deposit on the refined coarse level and then later needs to synchronize their data
  */
  Copier m_copierFiCoToFiCo;

  /*!
    @brief VoFIterator for fine-grid irregular ghost cells. 
  */
  mutable LayoutData<VoFIterator> m_vofIterFineGhosts;

  /*!
    @brief VoFIterator for coarse grid cells that are coarsenings of fine-grid ghost irregular cells.
    @details Essentially a coarsening of m_vofIterFineGhosts.
  */
  mutable LayoutData<VoFIterator> m_vofIterCoFiGhosts;

  /*!
    @brief VoFIterator for coarse cells
  */
  mutable LayoutData<VoFIterator> m_vofIterCoar;

  /*!
    @brief Stencils for arithmetic coarsening of volume data.
    @details This is defined on the coarse grid, and indexes into the refined coarse grid. 
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_cellArithmeticStencils;

  /*!
    @brief Stencils for conservative coarsening of volume data
    @details This is defined on the coarse grid, and indexes into the refined coarse grid.     
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_cellConservativeStencils;

  /*!
    @brief Refinement ratio between the two levels
  */
  int m_refRat;

  /*!
    @brief Number of ghost cells in data holders.
  */
  IntVect m_ghost;

  /*!
    @brief Define the vof iterators
  */
  void
  defineVoFIterators() noexcept;

  /*!
    @brief Define the vof iterators
  */
  void
  defineStencils() noexcept;

  /*!
    @brief Conservatively average and add data from the refined coarse data to the coarse data
    @param[inout] a_coarData Coarse-grid data
    @param[in] a_fineData Fine-grid data (must be on a refined layout of the coarse grid)
    @param[in] a_din Grid index
  */
  void
  conservativeAverageAndAdd(EBCellFAB& a_coarData, const EBCellFAB& a_fineData, const DataIndex& a_din) const noexcept;

  /*!
    @brief Arithmetically average and add data from the refined coarse data to the coarse data
    @param[inout] a_coarData Coarse-grid data
    @param[in] a_fineData Fine-grid data (must be on a refined layout of the coarse grid)
    @param[in] a_din Grid idnex
  */
  void
  arithmeticAverageAndAdd(EBCellFAB& a_coarData, const EBCellFAB& a_fineData, const DataIndex& a_din) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
