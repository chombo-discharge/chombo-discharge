/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_EBGradient.H
 * @brief  Declaration of a class which computes gradients in an EBAMR context.
 * @author Robert Marskar
 */

#ifndef CD_EBGRADIENT_H
#define CD_EBGRADIENT_H

// Std includes
#include <array>

// Chombo includes
#include <EBLevelGrid.H>
#include <AggStencil.H>
#include <FaceIterator.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Class for computing gradients in an EBAMR context.
 * @details Although computing the gradient might seem like an easy thing to do, the presence of a coarse-fine boundary
 * together with an EB complicates both the stencil availability and the consistent handling of the gradient on the
 * coarse level. This class introduces a two-level view of such stencil computations. This class uses regular second
 * order finite differencing in regular regions, including those covered by a finer level. Near the EB, however, the
 * stencil might reach into cut-cells that are covered by the finer level. The invalid data might be filled with bogus
 * data, or at the very best coarsened from the fine level. Since there is no consistent way of doing that coarsening
 * (both conserative and averaging make sense), this class switches to a least squares based approach for such cases. In
 * those particular corner cases we only reach into valid regions in the domain. Because of that, the directional
 * derivatives in those cells are in the form
 *
 *      d(phi)/dx = sum(vofCoar) w(vofCoar)*phi(vofCoar) + sum(vofFine)[w(vofFine)*phi(vofFine)]
 *
 * i.e., the gradient is a two-level stencil.
 *
 * In cut cells and next to the domain boundary, the derivative in a direction is a centered difference if the cells on
 * both sides are available, and otherwise a second order one-sided difference that reaches two cells away. Only if the
 * second cell is not available does the stencil drop to a first order one-sided difference. Near the refinement
 * boundary these stencils only use valid cells, i.e. cells that are not covered by the finer level. If that costs
 * order, the cell gets a two-level least squares stencil instead, fitted to the valid coarse cells within the reach of
 * the finite difference stencils and to the fine cells under the covered ones.
 * @note This class requires at least one ghost cell. The second order one-sided differences require two ghost cells,
 * both in the data and in the EBISLayout; with fewer, the one-sided differences are first order.
 */
class EBGradient
{
public:
  /**
   * @brief Specification of cell-center or cell-centroid type of interpolation
   */
  using CellLocation = Location::Cell;

  /**
   * @brief Disallowed constructor.
   */
  EBGradient() = delete;

  /**
   * @brief Disallowed copy constructor.
   */
  EBGradient(const EBGradient& a_other) = delete;

  /**
   * @brief Full constructor. Calsl the define function.
   * @param[in] a_eblg         Grids
   * @param[in] a_eblgFine     Fine grids
   * @param[in] a_eblgFiCo     Refinement of coarse grids
   * @param[in] a_hasFine      Has fine level or not
   * @param[in] a_dx           Resolution
   * @param[in] a_refRat       Refinement ratio
   * @param[in] a_order        Interpolation order for modified gradients.
   * @param[in] a_weighting    Weighting factor least squares.
   * @param[in] a_ghostVector Ghost vector
   */
  EBGradient(const EBLevelGrid& a_eblg,
             const EBLevelGrid& a_eblgFine,
             const EBLevelGrid& a_eblgFiCo,
             const bool         a_hasFine,
             const Real         a_dx,
             const int          a_refRat,
             const int          a_order,
             const int          a_weighting,
             const IntVect      a_ghostVector) noexcept;

  /**
   * @brief Destructor (does nothing)
   */
  virtual ~EBGradient() noexcept;

  /**
   * @brief Define function.
   * @param[in] a_eblg         Grids
   * @param[in] a_eblgFine     Fine grids
   * @param[in] a_eblgFiCo     Refinement of coarse grids
   * @param[in] a_hasFine      Has fine level or not
   * @param[in] a_dx           Resolution
   * @param[in] a_refRat       Refinement ratio
   * @param[in] a_order        Interpolation order for modified gradients.
   * @param[in] a_weighting    Weighting factor least squares.
   * @param[in] a_ghostVector Ghost vector
   */
  virtual void
  define(const EBLevelGrid& a_eblg,
         const EBLevelGrid& a_eblgFine,
         const EBLevelGrid& a_eblgFiCo,
         const bool         a_hasFine,
         const Real         a_dx,
         const int          a_refRat,
         const int          a_order,
         const int          a_weighting,
         const IntVect      a_ghostVector) noexcept;

  /**
   * @brief Compute gradient using data on the input level only.
   * @param[out] a_gradient Gradient of input scalar
   * @param[in]  a_phi      Input scalar
   */
  virtual void
  computeLevelGradient(LevelData<EBCellFAB>& a_gradient, const LevelData<EBCellFAB>& a_phi) const noexcept;

  /**
   * @brief Compute the gradient component which is normal to the face, using centered differencing.
   * @details If the face is oriented along dir, this computes the dir-component of the gradient using
   * centered differences on the face. This routine is useful after e.g. averaging a cell-centered gradient
   * to faces (in order to shrink the stencil).
   * @param[in,out] a_gradient Gradient of input scalar
   * @param[in]    a_phi      Input scalar
   */
  virtual void
  computeNormalDerivative(LevelData<EBFluxFAB>& a_gradient, const LevelData<EBCellFAB>& a_phi) const noexcept;

  /**
   * @brief Compute gradient using two-level stencils (matching at EBCF).
   * @param[out] a_gradient Gradient of input scalar
   * @param[in]  a_phi      Input scalar
   * @param[in]  a_phiFine  Input scalar on finer level.
   */
  virtual void
  computeAMRGradient(LevelData<EBCellFAB>&       a_gradient,
                     const LevelData<EBCellFAB>& a_phi,
                     const LevelData<EBCellFAB>& a_phiFine) const noexcept;

protected:
  /**
   * @brief Scalar component to reach into
   */
  static constexpr int m_comp = 0;

  /**
   * @brief Number of components that are supported
   */
  static constexpr int m_nComp = 1;

  /**
   * @brief How to interpret data locations.
   */
  CellLocation m_dataLocation;

  /**
   * @brief Grids
   */
  EBLevelGrid m_eblg;

  /**
   * @brief Fine grids
   */
  EBLevelGrid m_eblgFine;

  /**
   * @brief Refined coarse grids
   */
  EBLevelGrid m_eblgFiCo;

  /**
   * @brief Is defined or not
   */
  bool m_isDefined;

  /**
   * @brief True if we have a finer level (m_eblgFine is defined)
   */
  bool m_hasFine;

  /**
   * @brief True if we have EBCF crossings that require special stencils.
   */
  bool m_hasEBCF;

  /**
   * @brief Resolution on this level.
   */
  Real m_dx;

  /**
   * @brief Resolution on finer level.
   */
  Real m_dxFine;

  /**
   * @brief Desired order for least squares stencils.
   */
  int m_order;

  /**
   * @brief Weight (for least squares)
   */
  int m_weighting;

  /**
   * @brief Refinement factor between fine and coarse level.
   */
  int m_refRat;

  /**
   * @brief Number of ghost cells in incoming data holders
   */
  IntVect m_ghostVector;

  /**
   * @brief Reach of the finite difference stencils in cut cells and next to the domain boundary.
   * @details Two if the data and the EBISLayout have at least two ghost cells, in which case one-sided differences are
   * second order. Otherwise one, in which case one-sided differences are first order.
   */
  int m_stencilWidth;

  /**
   * @brief Pre-built copier for making copying from m_eblg to m_eblgFiCo go faster
   */
  Copier m_copier;

  /**
   * @brief Stencils for boundaries and cut-cells
   * @note Lives on the coarse layout.
   */
  LayoutData<BaseIVFAB<VoFStencil>> m_levelStencils;

  /**
   * @brief Part of EBCF gradient stencil which reaches into the coarse level.
   */
  LayoutData<BaseIVFAB<VoFStencil>> m_ebcfStencilsCoar;

  /**
   * @brief Part of EBCF gradient stencil which reaches into the coarse level.
   */
  LayoutData<BaseIVFAB<VoFStencil>> m_ebcfStencilsFine;

  /**
   * @brief AggStencils for boundaries and cut-cells. Provides faster application of m_levelStencils.
   * @note Lives on the coarse layout and only has one component.
   */
  LayoutData<RefCountedPtr<AggStencil<EBCellFAB, EBCellFAB>>> m_aggLevelStencils[SpaceDim];

  /**
   * @brief Iterator for boundary cells (either domain or EB)
   * @note Lives on the coarse layout.
   */
  mutable LayoutData<VoFIterator> m_levelIterator;

  /**
   * @brief Iterator for EBCF cells.
   * @note Lives on the coarse layout.
   */
  mutable LayoutData<VoFIterator> m_ebcfIterator;

  /**
   * @brief Face iterators restricted to multiply-cut faces per direction (valid box, SurroundingNoBoundary).
   * @details Used in computeFaceGradient as a second pass after the regular box loop to avoid
   * double-processing singly-cut faces.
   */
  mutable LayoutData<std::array<FaceIterator, SpaceDim>> m_multiCutFaceIter;

  /**
   * @brief Define level stencils.
   */
  virtual void
  defineLevelStencils() noexcept;

  /**
   * @brief Make masks
   * @param[out] a_coarMaskCF      Mask holding a value of 1 in coarse cells within m_stencilWidth cells of the
   * refinement boundary, outside the finer level.
   * @param[out] a_coarMaskInvalid Mask holding a value of 1 in coarse cells that are covered by a finer grid level. It
   * has m_stencilWidth ghost cells.
   */
  virtual void
  defineMasks(LevelData<FArrayBox>& a_coarMaskCF, LevelData<FArrayBox>& a_coarMaskInvalid) noexcept;

  /**
   * @brief Define iterators and cells that need modified stencils
   * @param[in]  a_coarMaskCF            Mask for CF region on the coarse layout.
   * @param[in]  a_coarMaskInvalid       Mask for invalid regions on the coarse layout.
   */
  virtual void
  defineIteratorsEBCF(const LevelData<FArrayBox>& a_coarMaskCF, const LevelData<FArrayBox>& a_coarMaskInvalid) noexcept;

  /**
   * @brief Define stencils for interpolation across the EBCF refinement boundary.
   * @param[in] a_coarMaskInvalid Coar mask invalid
   */
  virtual void
  defineStencilsEBCF(const LevelData<FArrayBox>& a_coarMaskInvalid) noexcept;

  /**
   * @brief Define buffer data
   */
  virtual void
  defineBuffers() noexcept;

  /**
   * @brief Get the finite difference stencil for the gradient in a cell, assuming that the data is cell-centered.
   * @details In each direction this is a centered difference if the cells on both sides are available. Otherwise it is
   * a one-sided difference, second order if m_stencilWidth is two and the second cell is available, and first order if
   * not. Cells in a_invalidCells are treated as unavailable.
   * @param[out] a_stencil      Finite difference stencil. Component dir holds the derivative in direction dir.
   * @param[in]  a_vof          Grid cell
   * @param[in]  a_ebisBox      EBISBox
   * @param[in]  a_invalidCells Cells that the stencil may not use. If nullptr, all cells can be used.
   * @return True if leaving out a_invalidCells did not lower the order of the derivative in any direction.
   */
  virtual bool
  getFiniteDifferenceStencil(VoFStencil&     a_stencil,
                             const VolIndex& a_vof,
                             const EBISBox&  a_ebisBox,
                             IntVectSet*     a_invalidCells) const noexcept;

  /**
   * @brief Build a gradient stencil using least squares.
   * @param[out] a_stencilCoar       Coarse-grid stencil
   * @param[out] a_stencilFine       Fine-grid stencil
   * @param[in]  a_vofCoar           VoF on the coarse level.
   * @param[in]  a_dataLocation      Data location
   * @param[in]  a_ebislCoar         EBISLayout on the coarse level
   * @param[in]  a_ebislFine         EBISLayout on the fine level
   * @param[in]  a_dit               Data index
   * @param[in]  a_validCellsCoar    Map of cells that are valid on the coarse domain.
   * @param[in]  a_validCellsFine    Map of cells that are valid on the fine domain.
   * @param[in]  a_dxCoar            Grid resolution on the coarse level
   * @param[in]  a_dxFine            Grid resolution on the fine level
   * @param[in]  a_order             Interpolation order
   * @param[in]  a_weight            Equation weighting for least squares reconstruction.
   * @return Returns true if a stencil was found, in which case it is defined in a_stencilFine and a_stencilCoar.
   */
  virtual bool
  getLeastSquaresStencil(VoFStencil&            a_stencilCoar,
                         VoFStencil&            a_stencilFine,
                         const VolIndex&        a_vofCoar,
                         const CellLocation&    a_dataLocation,
                         const EBISLayout&      a_ebislCoar,
                         const EBISLayout&      a_ebislFine,
                         const DataIndex&       a_dit,
                         const DenseIntVectSet& a_validCellsCoar,
                         const DenseIntVectSet& a_validCellsFine,
                         const Real&            a_dxCoar,
                         const Real&            a_dxFine,
                         const int&             a_order,
                         const int&             a_weight) const noexcept;

  /**
   * @brief Check if the level finite difference stencils in a cell only reach into valid cells.
   * @details This will check the stencils for all vofs in a_ivCoar. Stencils in cut cells may not reach into any
   * covered cell. Stencils in regular cells may reach into covered regular cells, whose data is the average of the
   * finer cells, but not into covered cut cells.
   * @param[in] a_ivCoar        Coarse grid cell.
   * @param[in] a_ebisBox       EBISBox
   * @param[in] a_invalidRegion Map of cells that are covered by a finer cell.
   * @return True if finite difference stencil valid, false otherwise
   */
  virtual bool
  isFiniteDifferenceStencilValid(const IntVect&   a_ivCoar,
                                 const EBISBox&   a_ebisBox,
                                 const FArrayBox& a_invalidRegion) const noexcept;

  /**
   * @brief Turn raw stencils into optimized AggStencils
   */
  virtual void
  makeAggStencils() noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
