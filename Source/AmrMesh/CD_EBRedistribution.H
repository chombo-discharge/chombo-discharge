/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBRedistribution.H
  @brief  Declaration of a class which can perform redistribution in an EB-AMR context.
  @author Robert Marskar
*/

#ifndef CD_EBRedistribution_H
#define CD_EBRedistribution_H

// Chombo includes
#include <EBLevelGrid.H>
#include <LevelData.H>
#include <EBCellFAB.H>
#include <BaseIVFAB.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for redistribution ala Chombo's flux redistribution. 
  @details This uses the same algorithm as the core Chombo code, but includes much faster define functions and smaller memory 
  footprint. Because this class is used in a non-subcycling context, it has way fewer bells and whistles. Note that this class
  permits redistribution outside of the domain, in which case some mass will (correctly) be lost.

  The underlying stencils for in this class can be understood as follows: Assume that we have some quantity deltaM_i that we
  want to redistribute from cell i to nearby cells. We fetch all cells that are within a radius of 1 around this cell. If one
  of these cells lie on other side of the refinement boundary, we coarsen it and fetch the corresponding coarse-grid cell. 
  Likewise, if a neighborhood cell lies underneath the fine grid, we refine it and redistribute to the fine grid cells. In each
  cut-cell we thus compute the total volume of the surrounding neighborhood as

     totalVolume = sum_(level_cells) kappa_j + 1/NrefToFine^D * sum_(fine_cells) kappa_j + Nref*ToCoar^D * sum(coar_cells) kappa_j.

  Every cell in the neighborhood of cell i gets a weight 1/totalVolume. Note that proper division by the grid resolution is not
  made here; the user will have to ensure proper scaling of his/her variables when calling the actual redistribution functions. 
  The scaling factors depend on what is actually being redistributed, and the user will need to work this out on paper. 
*/
class EBRedistribution
{
public:
  /*!
    @brief Weak constructor. Need to call define afterwards
  */
  EBRedistribution() noexcept;

  /*!
    @brief Copy constructor is not allowed.
  */
  EBRedistribution(const EBRedistribution& a_other) = delete;

  /*!
    @brief Full constructor. Calls the define function and puts object in usable state. 
    @details If there is a finer level, a_eblgFine and a_eblgRefined must be defined. If there is a coarser grid level then
    a_eblgCoar and a_eblgCoarsened must be defined. 
    @param[in] a_eblgCoar             Coarser grid level
    @param[in] a_eblgCoarsened        Coarsened grid
    @param[in] a_eblg                 Grids on this level
    @param[in] a_eblgRefinedd         Refined grid
    @param[in] a_eblgFine             Fine grid
    @param[in] a_refToCoar            Refinement factor between this level and the coarse level
    @param[in] a_refToFine            Refinement factor between this level and the fine level
    @param[in] a_redistributeOutside  Turn on/off redistribution outside domain. 
  */
  EBRedistribution(const EBLevelGrid& a_eblgCoar,
                   const EBLevelGrid& a_eblgCoarsened,
                   const EBLevelGrid& a_eblg,
                   const EBLevelGrid& a_eblgRefined,
                   const EBLevelGrid& a_eblgFine,
                   const int          a_refToCoar,
                   const int          a_refToFine,
                   const bool         a_redistributeOutside) noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBRedistribution() noexcept;

  /*!
    @brief Define fucntion. Puts object in usable state. 
    @details If there is a finer level, a_eblgFine and a_eblgRefined must be defined. If there is a coarser grid level then
    a_eblgCoar and a_eblgCoarsened must be defined. 
    @param[in] a_eblgCoar             Coarser grid level
    @param[in] a_eblgCoarsened        Coarsened grids
    @param[in] a_eblg                 Grids on this level
    @param[in] a_eblgRefined          Refined grids
    @param[in] a_eblgFine             Fine grid
    @param[in] a_refToCoar            Refinement factor between this level and the coarse level
    @param[in] a_refToFine            Refinement factor between this level and the fine level
    @param[in] a_redistributeOutside  Turn on/off redistribution outside domain. 
  */
  virtual void
  define(const EBLevelGrid& a_eblgCoar,
         const EBLevelGrid& a_eblgCoarsened,
         const EBLevelGrid& a_eblg,
         const EBLevelGrid& a_eblgRefined,
         const EBLevelGrid& a_eblgFine,
         const int          a_refToCoar,
         const int          a_refToFine,
         const bool         a_redistributeOutside) noexcept;

  /*!
    @brief Redistribute the input quantity into the coarse level, this level, and the fine level.
    @details If calling this, pass in a_phiFine/a_phiCoar = nullptr if there is no finer/coarser level. The a_scaleCoar and a_scaleFine arguments
    do not matter if there is no coarser/finer level.
    @param[inout] a_phiCoar   Data on coarse level
    @param[inout] a_phi       Data on this level
    @param[inout] a_phiFine   Data on finer level
    @param[in]    a_deltaM    Redistribution data on this level
    @param[in]    a_scaleCoar Scaling factor for redistribution into the coarse level
    @param[in]    a_scale     Scaling factor for redistribution into this level
    @param[in]    a_scaleFine Scaling factor for redistribution into the fine level
    @param[in]    a_variables Variables to redistribute
  */
  virtual void
  redistributeAMR(LevelData<EBCellFAB>*             a_phiCoar,
                  LevelData<EBCellFAB>*             a_phi,
                  LevelData<EBCellFAB>*             a_phiFine,
                  const LevelData<BaseIVFAB<Real>>& a_deltaM,
                  const Real                        a_scaleCoar,
                  const Real                        a_scale,
                  const Real                        a_scaleFine,
                  const Interval&                   a_variables) const noexcept;

  /*!
    @brief Redistribute the input quantity to the coarse level.
    @param[inout] a_phiCoar   Data on coarse level
    @param[in]    a_deltaM    Redistribution data on this level
    @param[in]    a_scaleCoar Scaling factor for redistribution into the coarse level
    @param[in]    a_variables Variables to redistribute
  */
  virtual void
  redistributeCoar(LevelData<EBCellFAB>&             a_phiCoar,
                   const LevelData<BaseIVFAB<Real>>& a_deltaM,
                   const Real&                       a_scaleCoar,
                   const Interval&                   a_variables) const noexcept;

  /*!
    @brief Redistribute the input quantity to this level.
    @param[inout] a_phi       Data on this level
    @param[in]    a_deltaM    Redistribution data on this level
    @param[in]    a_scale     Scaling factor for redistribution on this level
    @param[in]    a_variables Variables to redistribute
  */
  virtual void
  redistributeLevel(LevelData<EBCellFAB>&             a_phi,
                    const LevelData<BaseIVFAB<Real>>& a_deltaM,
                    const Real&                       a_scale,
                    const Interval&                   a_variables) const noexcept;

  /*!
    @brief Redistribute the input quantity to the coarse level.
    @param[inout] a_phiCoar   Data on fine level
    @param[in]    a_deltaM    Redistribution data on this level
    @param[in]    a_scaleFine Scaling factor for redistribution into the fine level
    @param[in]    a_variables Variables to redistribute
  */
  virtual void
  redistributeFine(LevelData<EBCellFAB>&             a_phiFine,
                   const LevelData<BaseIVFAB<Real>>& a_deltaM,
                   const Real&                       a_scaleFine,
                   const Interval&                   a_variables) const noexcept;

protected:
  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief If true, redistribute to the outside. Otherwise, do interior redistribution.
  */
  bool m_redistributeOutside;

  /*!
    @brief Has coarse level or not
  */
  bool m_hasCoar;

  /*!
    @brief Has fine level or not
  */
  bool m_hasFine;

  /*!
    @brief Refinement factor between this level and the coarse level
  */
  int m_refToCoar;

  /*!
    @brief Refinement factor between this level and the fine level
  */
  int m_refToFine;

  /*!
    @brief Redistribution radius. Always one.
  */
  int m_redistRadius;

  /*!
    @brief Grids on coarser level
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Grids on this level that are a refinement of the coarser level
  */
  EBLevelGrid m_eblgCoarsened;

  /*!
    @brief Grids on this level
  */
  EBLevelGrid m_eblg;

  /*!
    @brief Grids on this level that are a coarsening of the fine level
  */
  EBLevelGrid m_eblgRefined;

  /*!
    @brief Fine grids
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Copier for copying from the coarsening of this level to the coarse level. 
  */
  Copier m_coarCopier;

  /*!
    @brief Copier for adding data from this level to this level (yes, really).
  */
  Copier m_levelCopier;

  /*!
    @brief Copier for copying from the refinement of this level to the fine level.
  */
  Copier m_fineCopier;

  /*!
    @brief Stencils for redistribution into valid grid cells on the coarse level. 
    @note Defined on this level (m_eblg), but the stencils contain VoFs on the coarse level
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_redistStencilsCoar;

  /*!
    @brief Stencils for redistribution into valid grid cells on this level. 
    @note Defined on this level (m_eblg), and the stencils contain VoFs on this level also.
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_redistStencilsLevel;

  /*!
    @brief Stencils for redistribution into valid grid cells on the fine level. 
    @note Defined on this level (m_eblg), but the stencils contain VoFs on the fine level
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_redistStencilsFine;

  /*!
    @brief Iterator for going through all cells on this level that we redistribute from
  */
  mutable LayoutData<VoFIterator> m_vofit;

  /*!
    @brief Define redistribution stencils
  */
  virtual void
  defineStencils() noexcept;

  /*!
    @brief Define buffer storages
  */
  virtual void
  defineBuffers() noexcept;

  /*!
    @brief Define the valid region of this level.
    @param[inout] a_validCells Valid grid cells on this level. 
  */
  virtual void
  defineValidCells(LevelData<BaseFab<bool>>& a_validCells) const noexcept;

  /*!
    @brief Define the cells on the coarse side of the CF interface that we would redistribute to. 
    @param[inout] a_validCells Valid grid cells on this level. 
  */
  virtual void
  defineInterfaceCells(LevelData<BaseFab<bool>>& a_interfaceCells) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
