/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBMultigridInterpolator.H
  @brief  Parent base class for multigrid interpolation. 
  @author Robert Marskar
*/

#ifndef CD_EBMultigridInterpolator_H
#define CD_EBMultigridInterpolator_H

// Chombo includes
#include <EBLevelGrid.H>
#include <EBCellFAB.H>
#include <DerivStencil.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Multigrid interpolator parent class. Virtual beause I want different implementations to do different things. 
*/
class EBMultigridInterpolator
{
public:
  /*!
    @brief Basic constructor. 
  */
  EBMultigridInterpolator();

  /*!
    @brief Disallowed copy constructor.
  */
  EBMultigridInterpolator(const EBMultigridInterpolator& a_other) = delete;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBMultigridInterpolator();

  /*!
    @brief Return number of ghost cells across CF interface.
  */
  virtual int
  getGhostCF() const = 0;

  /*!
    @brief Get an explicit stencil for a ghost cell near the EB.
    @details This should return the fine + coarse stencils required for interpolating to the fine ghost cell.
    @param[in] a_fineGhost Ghost cell on the fine side of the refinement boundary.
    @param[in] a_din Grid index.
    @return This first/second entry in the returned pair should be the fine/coarse stencil.
  */
  virtual std::pair<VoFStencil, VoFStencil>
  getInterpolationStencilEB(const VolIndex& a_fineGhost, const DataIndex& a_din) const noexcept = 0;

  /*!
    @brief Get the interpolation stencil near the refinement boundary.
    @details This version can assume that all cells are regular, and that returned stencils may later be replaced
    by getInterpolationStencilEB.
    @param[in] a_fineGhost Ghost cell index on the fine grid
    @param[in] a_din Grid index
    @param[in] a_dir Coordinate direction
    @param[in] a_side Low or high side of the grid patch
  */
  virtual std::pair<DerivStencil, DerivStencil>
  getInterpolationStencilRegular(const IntVect&       a_fineGhost,
                                 const DataIndex&     a_din,
                                 const int            a_dir,
                                 const Side::LoHiSide a_side) const noexcept = 0;

  /*!
    @brief Do inhomogeneous interpolation
    @param[inout] a_phiFine   Fine phi
    @param[in]    a_phiCoar   Coarse phi
    @param[in]    a_variables Variables to interpolate
  */
  virtual void
  coarseFineInterp(LevelData<EBCellFAB>&       a_phiFine,
                   const LevelData<EBCellFAB>& a_phiCoar,
                   const Interval              a_variables) const = 0;

  /*!
    @brief Do homogeneous interpolation
    @param[in] a_phiFine   Fine phi
    @param[in] a_variables Variables to interpolate
  */
  virtual void
  coarseFineInterpH(LevelData<EBCellFAB>& a_phiFine, const Interval a_variables) const = 0;

  /*!
    @brief Do homogeneous coarse-fine interpolation on a grid patch
    @param[in] a_phiFine   Fine-level data
    @param[in] a_variables Variables to be interpolated
    @param[in] a_din       Data index
    @note This is the kernel version. The constructor guarantees that the fine-grid stencil will not out of the ghost region defined by a_ghostVector. However, this is
    the direct kernel version which is public because multigrid operators might need it for performance reasons. This routine does NOT warn against input data that do not
    contain enough ghost cells. 
  */
  virtual void
  coarseFineInterpH(EBCellFAB& a_phiFine, const Interval a_variables, const DataIndex& a_din) const = 0;
};

#include <CD_NamespaceFooter.H>

#endif
