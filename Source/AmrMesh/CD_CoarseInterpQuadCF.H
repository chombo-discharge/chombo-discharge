/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CoarseInterpQuadCF.H
  @brief  Supporting class for doing the coarse-side interpolation for multigrid
  @author Robert Marskar
*/

#ifndef CD_CoarseInterpQuadCF_H
#define CD_CoarseInterpQuadCF_H

// Chombo includes
#include <IntVect.H>
#include <FArrayBox.H>
#include <LoHiSide.H>
#include <DerivStencil.H>
#include <DisjointBoxLayout.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for computing the required Taylor terms for doing the orthogonal extrapolation to the ghost cell position. This
  is used in conjuction with EBLeastSquaresMultigridInterpolator.H, and should be understood in that context only. 

  @details This should be defined over the fine grid. It will store stencils for computing the various required terms in the 
  Taylor series, and the user can call various functions for getting those terms. 

  @todo This class stores all stencils explicitly. We might want to consider a "default case" (e.g., centered differencing) for compressing memory.
*/
class CoarseInterpQuadCF
{
public:
  /*!
    @brief Default constructor. Must subsequently call the define function
  */
  CoarseInterpQuadCF() noexcept;

  /*!
    @brief Destructor.
  */
  virtual ~CoarseInterpQuadCF() noexcept;

  /*!
    @brief Define function. Puts object in usable state.
    @param[in] a_dblFine Fine grids
    @param[in] a_domainCoar Coarse domain
    @param[in] a_dit Data index on the fine grid
    @param[in] a_fineGhostCells Ghost cells to interpolate (on the fine level)
    @param[in] a_refRat Refinement factor between coarse and fine
    @param[in] a_ignoreDir Coordinate to ignore when doing the interpolation. I.e., we interpolate 
    in all coordinates except this one.
  */
  virtual void
  define(const DisjointBoxLayout& a_dblFine,
         const ProblemDomain&     a_domainCoar,
         const DataIndex&         a_dit,
         const Box&               a_fineGhostCells,
         const int                a_refRat,
         const int                a_ignoreDir) noexcept;

  /*!
    @brief Compute first-derivative Taylor term in direction a_dir. 
    @param[in] a_coarPhi Coarse-side phi. 
    @param[in] a_ivCoar Cell where we compute.
    @param[in] a_dir Direction of derivative. Must be different from m_ignoreDir
    @param[in] a_coarVar a_dir Direction of derivative. Must be different from m_ignoreDir
  */
  virtual Real
  computeFirstDeriv(const FArrayBox& a_coarPhi,
                    const IntVect&   a_ivCoar,
                    const int        a_dir,
                    const int        a_coarVar) const noexcept;

  /*!
    @brief Compute second-derivative Taylor term in direction a_dir. 
    @param[in] a_coarPhi Coarse-side phi. 
    @param[in] a_ivCoar Cell where we compute.
    @param[in] a_dir Direction of derivative. Must be different from m_ignoreDir
    @param[in] a_coarVar a_dir Direction of derivative. Must be different from m_ignoreDir
  */
  virtual Real
  computeSecondDeriv(const FArrayBox& a_coarPhi,
                     const IntVect&   a_ivCoar,
                     const int        a_dir,
                     const int        a_coarVar) const noexcept;

  /*!
    @brief Compute mixed-derivative Taylor term in directions orthogonal m_ignoreDir.
    @param[in] a_coarPhi Coarse-side phi. 
    @param[in] a_ivCoar Cell where we compute.
    @param[in] a_coarVar a_dir Direction of derivative. Must be different from m_ignoreDir
  */
  virtual Real
  computeMixedDeriv(const FArrayBox& a_coarPhi, const IntVect& a_ivCoar, const int a_coarVar) const noexcept;

protected:
  /*!
    @brief Simple enum for holding stencil types for the first derivative
  */
  enum class FirstDerivStencil : short
  {
    Centered2,
    Backward2,
    Forward2,
    Backward1,
    Forward1,
  };
  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Fine grids
  */
  DisjointBoxLayout m_dblFine;

  /*!
    @brief Coarse domain
  */
  ProblemDomain m_domainCoar;

  /*!
    @brief Data index in the fine grid
  */
  DataIndex m_dit;

  /*!
    @brief Interpolation direction to ignore. 
  */
  int m_ignoreDir;

  /*!
    @brief Refinement factor to coarse. 
  */
  int m_refRat;

  /*!
    @brief Coarse-grid cells corresponding to the fine-grid ghost cells
    @note This is the coarsened box of the fine-grid ghost cell layer, and is where we've defined all our stencils. 
  */
  Box m_stencilBox;

  /*!
    @brief Stencil types for first derivative
  */
  BaseFab<FirstDerivStencil> m_firstDerivStencils[SpaceDim];

  /*!
    @brief Second-derivative stencils (in all directions except m_dir)
  */
  BaseFab<DerivStencil> m_secondDerivStencils[SpaceDim];

  /*!
    @brief Mixed-derivative stencil. Only defined in 3D.
  */
  BaseFab<DerivStencil> m_mixedDerivStencils;

  /*!
    @brief Compute all first-derivative stencils in the coarse grid cells.
    @note Must be accurate to O(h^2).
  */
  virtual void
  defineFirstDerivStencils() noexcept;

  /*!
    @brief Compute all second-derivative stencils in the coarse grid cells. 
    @details This only needs to accurate to O(h).
  */
  virtual void
  defineSecondDerivStencils() noexcept;

  /*!
    @brief Compute all second-derivative stencils in the coarse grid cells. 
    @details This only needs to accurate to O(h).
  */
  virtual void
  defineMixedDerivStencils() noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
