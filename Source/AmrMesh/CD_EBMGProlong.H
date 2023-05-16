/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBMGProlong.H
  @brief  Declaration of a prolongation operator for EB geometric multigrid
  @author Robert Marskar
*/

#ifndef CD_EBMGProlong_H
#define CD_EBMGProlong_H

// Chombo includes
#include <EBLevelGrid.H>
#include <ProblemDomain.H>
#include <LevelData.H>
#include <EBCellFAB.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for prolongation of multigrid residual onto a finer grid. 
  @note This class only does constant prolongation. Using Chombos EBMGInterp with piecewise linear (bi/tri-linear) did not
  improve multigrid convergence. 
*/
class EBMGProlong
{
public:
  /*!
    @brief Default constructor. User must subsequently call the define function
  */
  EBMGProlong() noexcept;

  /*!
    @brief Full constructor. Calls the define fucntion
    @param[in] a_eblgFine Fine level that we prolong to.
    @param[in] a_eblgCoar Coarse level that we prolong from.
    @param[in] a_refRat   Refinement factor between levels.
  */
  EBMGProlong(const EBLevelGrid& a_eblgFine, const EBLevelGrid& a_eblgCoar, const int& a_refRat) noexcept;

  /*!
    @brief Destructor
  */
  virtual ~EBMGProlong() noexcept;

  /*!
    @brief Define function.
    @param[in] a_eblgFine Fine level that we prolong to.
    @param[in] a_eblgCoar Coarse level that we prolong from.
    @param[in] a_refRat   Refinement factor between levels.
  */
  virtual void
  define(const EBLevelGrid& a_eblgFine, const EBLevelGrid& a_eblgCoar, const int& a_refRat) noexcept;

  /*!
    @brief Prolong residual onto the coarse grid. 
    @param[inout] a_fineData  Fine-grid residual.
    @param[in]    a_coarData  Coarse-grid residual
    @param[in]    a_variables Variables to prolong.
  */
  virtual void
  prolongResidual(LevelData<EBCellFAB>&       a_fineData,
                  const LevelData<EBCellFAB>& a_coarData,
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
    @brief Iterator for iterating over cut-cells on the fine level
  */
  mutable LayoutData<VoFIterator> m_vofitFine;

  /*!
    @brief Prolongion stencils near the EB.
    @details The stencils are defined on the fine grid but their vofs reach into the coarse grid. 
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_prolongStencils;
};

#include <CD_NamespaceFooter.H>

#endif
