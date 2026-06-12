/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_CoarseInterpQuadCF.H
  @brief  Supporting class for doing the coarse-side interpolation for multigrid
  @author Robert Marskar
*/

#ifndef CD_COARSEINTERPQUADCF_H
#define CD_COARSEINTERPQUADCF_H

// Chombo includes
#include <IntVect.H>
#include <FArrayBox.H>
#include <LoHiSide.H>
#include <DerivStencil.H>
#include <DisjointBoxLayout.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Class for computing the required Taylor terms for doing the orthogonal extrapolation to the ghost cell position. This
  is used in conjunction with EBLeastSquaresMultigridInterpolator.H, and should be understood in that context only. 

  @details This should be defined over the fine grid. It will store stencils for computing the various required terms in the 
  Taylor series, and the user can call various functions for getting those terms. For the record, all stencils in this object
  assumes that dx = 1, so the user have to normalize against his/her grid resolutions later. 
*/
class CoarseInterpQuadCF
{
public:
  /**
    @brief Default constructor. Must subsequently call the define function
  */
  CoarseInterpQuadCF() noexcept;

  /**
    @brief Destructor.
  */
  virtual ~CoarseInterpQuadCF() noexcept;

  /**
    @brief Copy constructor — deleted.
  */
  CoarseInterpQuadCF(const CoarseInterpQuadCF&) = delete;

  /**
    @brief Copy assignment — deleted.
  */
  CoarseInterpQuadCF&
  operator=(const CoarseInterpQuadCF&) = delete;

  /**
    @brief Move constructor.
  */
  CoarseInterpQuadCF(CoarseInterpQuadCF&&) = default;

  /**
    @brief Move assignment.
    @return Reference to this object.
  */
  CoarseInterpQuadCF&
  operator=(CoarseInterpQuadCF&&) noexcept = default;

  /**
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
         int                      a_refRat,
         int                      a_ignoreDir) noexcept;

  /**
    @brief Compute first-derivative Taylor term in direction a_dir. 
    @param[in] a_coarPhi Coarse-side phi. 
    @param[in] a_ivCoar Cell where we compute.
    @param[in] a_dir Direction of derivative. Must be different from m_ignoreDir
    @param[in] a_coarVar a_dir Direction of derivative. Must be different from m_ignoreDir
    @return Computed first deriv
  */
  virtual Real
  computeFirstDeriv(const FArrayBox& a_coarPhi, const IntVect& a_ivCoar, int a_dir, int a_coarVar) const noexcept;

  /**
    @brief Compute second-derivative Taylor term in direction a_dir.
    @param[in] a_coarPhi Coarse-side phi. 
    @param[in] a_ivCoar Cell where we compute.
    @param[in] a_dir Direction of derivative. Must be different from m_ignoreDir
    @param[in] a_coarVar a_dir Direction of derivative. Must be different from m_ignoreDir
    @return Computed second deriv
  */
  virtual Real
  computeSecondDeriv(const FArrayBox& a_coarPhi, const IntVect& a_ivCoar, int a_dir, int a_coarVar) const noexcept;

  /**
    @brief Compute mixed-derivative Taylor term in directions orthogonal m_ignoreDir.
    @param[in] a_coarPhi Coarse-side phi. 
    @param[in] a_ivCoar Cell where we compute.
    @param[in] a_coarVar a_dir Direction of derivative. Must be different from m_ignoreDir
    @return Computed mixed deriv
  */
  virtual Real
  computeMixedDeriv(const FArrayBox& a_coarPhi, const IntVect& a_ivCoar, int a_coarVar) const noexcept;

  /**
    @brief Get the first derivative stencil for the input coarse vof.
    @details This will construct the stencil for representing the derivative on the coarse cell center
    and return it. Useful when explicit stencil storage is required. 
    @param[in] a_ivCoar Coarse-grid cell index (i.e., the one underneath the fine-grid ghost)
    @param[in] a_dir Direction of derivative. 
    @return First deriv stencil
  */
  virtual DerivStencil
  getFirstDerivStencil(const IntVect& a_ivCoar, int a_dir) const noexcept;

  /**
    @brief Get the second derivative stencil for the input coarse vof.
    @details This will construct the stencil for representing the derivative on the coarse cell center
    and return it. Useful when explicit stencil storage is required.     
    @param[in] a_ivCoar Coarse-grid cell index (i.e., the one underneath the fine-grid ghost)
    @param[in] a_dir Direction of derivative.     
    @return Second deriv stencil
  */
  virtual DerivStencil
  getSecondDerivStencil(const IntVect& a_ivCoar, int a_dir) const noexcept;

  /**
    @brief Get the second derivative stencil for the input coarse vof.
    @details This will construct the stencil for representing the derivative on the coarse cell center
    and return it. Useful when explicit stencil storage is required.     
    @param[in] a_ivCoar Coarse-grid cell index (i.e., the one underneath the fine-grid ghost)
    @return Mixed deriv stencil
  */
  virtual DerivStencil
  getMixedDerivStencil(const IntVect& a_ivCoar) const noexcept;

protected:
  /**
    @brief Simple enum for holding stencil types for the first derivative
  */
  enum class FirstDerivStencil : unsigned short
  {
    Centered2,
    Backward2,
    Forward2,
    Backward1,
    Forward1,
  };

  /**
    @brief Simple enum for holding stencil types for the second derivative
  */
  enum class SecondDerivStencil : unsigned short
  {
    Centered2,
    Backward1,
    Forward1
  };

  /**
    @brief Simple enum for holding stencil types for the mixed derivative
  */
  enum class MixedDerivStencil : unsigned short
  {
    Standard,
    Explicit
  };

  /**
    @brief Is defined or not
  */
  bool m_isDefined;

  /**
    @brief Fine grids
  */
  DisjointBoxLayout m_dblFine;

  /**
    @brief Coarse domain
  */
  ProblemDomain m_domainCoar;

  /**
    @brief Data index in the fine grid
  */
  DataIndex m_dit;

  /**
    @brief Interpolation direction to ignore. 
  */
  int m_ignoreDir;

  /**
    @brief One tangential direction
  */
  int m_tanDir1;

  /**
    @brief One tangential direction
  */
  int m_tanDir2;

  /**
    @brief Refinement factor to coarse. 
  */
  int m_refRat;

  /**
    @brief Coarse-grid cells corresponding to the fine-grid ghost cells
    @note This is the coarsened box of the fine-grid ghost cell layer, and is where we've defined all our stencils. 
  */
  Box m_stencilBox;

  /**
    @brief Stencil types for first derivative
  */
  BaseFab<FirstDerivStencil> m_firstDerivStencils[SpaceDim];

  /**
    @brief Second-derivative stencils (in all directions except m_dir)
  */
  BaseFab<SecondDerivStencil> m_secondDerivStencils[SpaceDim];

  /**
    @brief Identifier for whether or not we should use a standard or explicit mixed derivative stencil
  */
  BaseFab<MixedDerivStencil> m_mixedDerivStencils;

  /**
    @brief Mixed-derivative stencil. Only defined in 3D.
  */
  BaseFab<DerivStencil> m_explicitMixedDerivStencils;

  /**
    @brief Compute all first-derivative stencils in the coarse grid cells.
    @note Must be accurate to O(h^2).
  */
  virtual void
  defineStencils() noexcept;

  /**
    @brief Compute all second-derivative stencils in the coarse grid cells. 
    @details This only needs to accurate to O(h).
  */
  virtual void
  defineMixedDerivStencils() noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
