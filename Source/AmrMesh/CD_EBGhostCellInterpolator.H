/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBGhostCellInterpolator.H
  @brief  Declaration of a class which interpolates ghost cells across the coarse-fine interface. 
  @author Robert Marskar
*/

#ifndef CD_EBGhostCellInterpolator_H
#define CD_EBGhostCellInterpolator_H

// Std includes
#include <map>

// Chombo includes
#include <EBLevelGrid.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which can interpolate ghost cells across the coarse-fine interface. To select an interpolation method, use the type argument. 

  @details This class works by computing slopes on the coarse grid. We then interpolate within the coarse-grid cell to the fine-grid ghost cell, using
  the computed slopes. Note that all fine-grid ghost cells are interpolated from the coarse grid. Thus, this interpolator is useless for multigrid. 

  @todo The cut-cell slope calculation and interpolation uses on-the-fly coarsening of the ghost cell. This could be made more efficient but since there
  are few of these cells I don't think it is necessary (but only benchmark timings will answer this quantitatively).
*/
class EBGhostCellInterpolator
{
public:
  /*!
    @brief Type of interpolation slopes supported.
  */
  enum Type
  {
    PWC,
    MinMod,
    MonotonizedCentral,
    Superbee
  };

  /*!
    @brief Disallowed constructor. 
  */
  EBGhostCellInterpolator() noexcept;

  /*!
    @brief Disallowed copy constructor.
  */
  EBGhostCellInterpolator(const EBGhostCellInterpolator& a_other) = delete;

  /*!
    @brief Full constructor. Calls the define function. 
    @param[in] a_eblgFine     Fine grids
    @param[in] a_eblgCoFi     Coarsened fine grids
    @param[in] a_eblgCoar     Coarse grids
    @param[in] a_dataLocation Interpration of data centering. Either on the cell center or the cell centroid. 
    @param[in] a_ghostVector  (Minimum) number of ghost cells in input/output data. 
    @param[in] a_refRat       Refinement ratio between coarse and fine grids.
    @param[in] a_ghostCF      Width of the ghost region to be filled. Only relevant near the EBCF. 
  */
  EBGhostCellInterpolator(const EBLevelGrid& a_eblgFine,
                          const EBLevelGrid& a_eblgCoFi,
                          const EBLevelGrid& a_eblgCoar,
                          const IntVect&     a_ghostVector,
                          const int          a_refRat,
                          const int          a_ghostCF) noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBGhostCellInterpolator() noexcept;

  /*!
    @brief Define method. Puts object in usable state. 
    @param[in] a_eblgFine     Fine grids
    @param[in] a_eblgCoFi     Coarsened fine grids
    @param[in] a_eblgCoar     Coarse grids
    @param[in] a_dataLocation Interpration of data centering. Either on the cell center or the cell centroid. 
    @param[in] a_ghostVector  (Minimum) number of ghost cells in input/output data. 
    @param[in] a_refRat       Refinement ratio between coarse and fine grids.
    @param[in] a_ghostCF      Width of the ghost region to be filled. Only relevant near the EBCF. 
  */
  virtual void
  define(const EBLevelGrid& a_eblgFine,
         const EBLevelGrid& a_eblgCoFi,
         const EBLevelGrid& a_eblgCoar,
         const IntVect&     a_ghostVector,
         const int          a_refRat,
         const int          a_ghostCF) noexcept;

  /*!
    @brief Do inhomogeneous interpolation
    @param[inout] a_phiFine   Fine phi
    @param[in]    a_phiCoar   Coarse phi
    @param[in]    a_variables Variables to interpolate
    @param[in]    a_type      Interpolation type
  */
  virtual void
  interpolate(LevelData<EBCellFAB>&       a_phiFine,
              const LevelData<EBCellFAB>& a_phiCoar,
              const Interval              a_variables,
              const Type                  a_interpType) const noexcept;

protected:
  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Fine grids
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Coarse grids
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Coarsened fine grids. 
  */
  EBLevelGrid m_eblgCoFi;

  /*!
    @brief Minimum number of ghost cells in input data. 
  */
  IntVect m_ghostVector;

  /*!
    @brief Refinement factor between fine and coarse level. 
  */
  int m_refRat;

  /*!
    @brief Number of ghost cells to fill across coarse-fine interface. 
    @note Only applies the the cut-cells. 
  */
  int m_ghostCF;

  /*!
    @brief Regular ghost regions to be interpolated. 
  */
  LayoutData<std::map<std::pair<int, Side::LoHiSide>, Box>> m_regularGhostRegions;

  /*!
    @brief Fine-grid ghost cells on the coarse-fine interface that overlap an irregular coarse cell. 
  */
  mutable LayoutData<VoFIterator> m_fineIrregCells;

  /*!
    @brief Coarse-grid cells on the caorse-side of the interface that are irregular cells.
  */
  mutable LayoutData<VoFIterator> m_coarIrregCells;

  /*!
    @brief Buffer data defined on the coarsened fine grids (over m_coarIrregCells). Used for holding slopes in the
    coarse-grid cut-cells
  */
  mutable LayoutData<BaseIVFAB<Real>> m_coarIrregSlopes;

  /*!
    @brief Buffer data
  */
  mutable LevelData<EBCellFAB>* m_grownCoarData;

  /*!
    @brief Define ghost regions to be interpolated
  */
  virtual void
  defineGhostRegions() noexcept;

  /*!
    @brief Regular interpolation, as if the EB was not there. 
    @param[in] a_phiFine Fine-grid data
    @param[in] a_phiCoar Coarse-grid data
    @param[in] a_dit Grid index on the fine grid
    @param[in] a_fineVar Fine-grid variable
    @param[in] a_coarVar Coarse-grid variable
    @param[in] a_interType Slope type
  */
  virtual void
  interpolateRegular(FArrayBox&       a_phiFine,
                     const FArrayBox& a_phiCoar,
                     const DataIndex& a_dit,
                     const int        a_fineVar,
                     const int        a_coarVar,
                     const Type       a_interpType) const noexcept;

  /*!
    @brief Ghost cell interpolation near the EB, correcting slopes that broke in interpolateRegular. 
    @param[in] a_phiFine Fine-grid data
    @param[in] a_phiCoar Coarse-grid data
    @param[in] a_dit Grid index on the fine grid
    @param[in] a_fineVar Fine-grid variable
    @param[in] a_coarVar Coarse-grid variable
    @param[in] a_interType Slope type
  */
  virtual void
  interpolateIrregular(EBCellFAB&       a_phiFine,
                       const EBCellFAB& a_phiCoar,
                       const DataIndex& a_dit,
                       const int        a_fineVar,
                       const int        a_coarVar,
                       const Type       a_interpType) const noexcept;

  /*!
    @brief minmod slope function. 
    @param[in] dwl Left slope
    @param[in] dwr Right slope
  */
  inline Real
  minmod(const Real& dwl, const Real& dwr) const noexcept;

  /*!
    @brief Superbee slope limiter. 
    @param[in] dwl Left slope
    @param[in] dwr Right slope
  */
  inline Real
  superbee(const Real& dwl, const Real& dwr) const noexcept;

  /*!
    @brief Monotonized central difference slope limiter.
    @param[in] dwl Left slope
    @param[in] dwr Right slope
  */
  inline Real
  monotonizedCentral(const Real& dwl, const Real& dwr) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
