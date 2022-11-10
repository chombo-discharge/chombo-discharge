/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBMGRestrict.H
  @brief  Declaration of a restriction operator for EB geometric multigrid
  @author Robert Marskar
*/

#ifndef CD_EBMGRestrict_H
#define CD_EBMGRestrict_H

// Chombo includes
#include <EBLevelGrid.H>
#include <ProblemDomain.H>
#include <LevelData.H>
#include <EBCellFAB.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for restricting multigrid residual onto a coarser grid.
*/
class EBMGRestrict
{
public:
  /*!
    @brief Default constructor. User must subsequently call the define function
  */
  EBMGRestrict() noexcept;

  /*!
    @brief Full constructor. Calls the define fucntion
  */
  EBMGRestrict(const EBLevelGrid&   a_eblgFine,
               const EBLevelGrid&   a_eblgCoar,
               const int&           a_refRat) noexcept;

  /*!
    @brief Destructor
  */
  virtual ~EBMGRestrict() noexcept;

  /*!
    @brief Define function.
  */
  virtual void
  define(const EBLevelGrid&   a_eblgFine,
         const EBLevelGrid&   a_eblgCoar,
         const int&           a_refRat) noexcept;

  /*!
    @brief Restrict residual onto the coarse grid
  */
  virtual void
  restrict(LevelData<EBCellFAB>&       a_coarData,
           const LevelData<EBCellFAB>& a_fineData,
           const Interval              a_variables) const noexcept;

protected:
  /*!
    @brief Defined or not
  */
  bool m_isDefined;

  /*!
    @brief Refinement ratio
  */
  int m_refRat;

  /*!
    @brief Fine grid
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Coarse grid
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Coarsened fine grids
  */
  EBLevelGrid m_eblgCoFi;

  /*!
    @brief Buffer data on the coarsened fine grids
  */
  mutable LevelData<EBCellFAB> m_dataCoFi;

  /*!
    @brief Iterator for iterating over cut-cells on the coarse level
  */
  mutable LayoutData<VoFIterator> m_vofitCoar;

  /*!
    @brief Restriction stencils near the EB.
    @details The stencils are defined on the coarse-grid, but reach into the fine grid.
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_restrictStencils;
};

#include <CD_NamespaceFooter.H>

#endif
