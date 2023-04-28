/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBGradient.H
  @brief  Declaration of a class which computes gradients in an EBAMR context. 
  @author Robert Marskar
*/

#ifndef CD_EBGradient_H
#define CD_EBGradient_H

// Chombo includes
#include <EBLevelGrid.H>
#include <AggStencil.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for computing gradients in an EBAMR context.
  @details Although computing the gradient might seem like an easy thing to do, the presence of a coarse-fine boundary together with an EB complicates
  both the stencil availability and the consistent handling of the gradient on the coarse level. This class introduces a two-level view of such stencil
  computations. This class uses regular second order finite differencing in regular regions, including those covered by a finer level. Near the EB, however,
  the stencil might reach into cut-cells that are covered by the finer level. The invalid data might be filled with bogus data, or at the very best coarsened
  from the fine level. Since there is no consistent way of doing that coarsening (both conserative and averaging make sense), this class switches to a least
  squares based approach for such cases. In those particular corner cases we only reach into valid regions in the domain. Because of that, the directional
  derivatives in those cells are in the form

       d(phi)/dx = sum(vofCoar) w(vofCoar)*phi(vofCoar) + sum(vofFine)[w(vofFine)*phi(vofFine)]

  i.e., the gradient is a two-level stencil.

  Near the domain boundaries this routine switches to first order stencils. 
  @note If there is no finer level, this class will compute regular second-order stencils away from the EB and switch to compact width-one stencils if
  there are not enough cells available near the EB. 
  @note This class requires at least one ghost cell. 
*/
class EBGradient
{
public:
  /*!
    @brief Specification of cell-center or cell-centroid type of interpolation
  */
  using CellLocation = Location::Cell;

  /*!
    @brief Disallowed constructor. 
  */
  EBGradient() = delete;

  /*!
    @brief Disallowed copy constructor.
  */
  EBGradient(const EBGradient& a_other) = delete;

  /*!
    @brief Full constructor. Calsl the define function.
    @param[in] a_eblg         Grids
    @param[in] a_eblgFine     Fine grids
    @param[in] a_eblgFiCo     Refinement of coarse grids
    @param[in] a_hasFine      Has fine level or not
    @param[in] a_dx           Resolution
    @param[in] a_refRat       Refinement ratio
    @param[in] a_order        Interpolation order for modified gradients. 
    @param[in] a_weighting    Weighting factor least squares. 
    @param[in] a_ghostVect    Number of ghost cells in the input data
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

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBGradient() noexcept;

  /*!
    @brief Define function.
    @param[in] a_eblg         Grids
    @param[in] a_eblgFine     Fine grids
    @param[in] a_eblgFiCo     Refinement of coarse grids
    @param[in] a_hasFine      Has fine level or not
    @param[in] a_dx           Resolution
    @param[in] a_refRat       Refinement ratio
    @param[in] a_order        Interpolation order for modified gradients. 
    @param[in] a_weighting    Weighting factor least squares. 
    @param[in] a_ghostVect    Number of ghost cells in the input data
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

  /*!
    @brief Compute gradient using data on the input level only. 
    @param[out] a_gradient Gradient of input scalar
    @param[in]  a_phi      Input scalar
  */
  virtual void
  computeLevelGradient(LevelData<EBCellFAB>& a_gradient, const LevelData<EBCellFAB>& a_phi) const noexcept;

  /*!
    @brief Compute the gradient component which is normal to the face, using centered differencing. 
    @details If the face is oriented along dir, this computes the dir-component of the gradient using
    centered differences on the face. This routine is useful after e.g. averaging a cell-centered gradient
    to faces (in order to shrink the stencil). 
    @param[inout] a_gradient Gradient of input scalar
    @param[in]    a_phi      Input scalar
  */
  virtual void
  computeNormalDerivative(LevelData<EBFluxFAB>& a_gradient, const LevelData<EBCellFAB>& a_phi) const noexcept;

  /*!
    @brief Compute gradient using two-level stencils (matching at EBCF). 
    @param[out] a_gradient Gradient of input scalar
    @param[in]  a_phi      Input scalar
    @param[in]  a_phiFine  Input scalar on finer level. 
  */
  virtual void
  computeAMRGradient(LevelData<EBCellFAB>&       a_gradient,
                     const LevelData<EBCellFAB>& a_phi,
                     const LevelData<EBCellFAB>& a_phiFine) const noexcept;

protected:
  /*!
    @brief Scalar component to reach into
  */
  static constexpr int m_comp = 0;

  /*!
    @brief Number of components that are supported
  */
  static constexpr int m_nComp = 1;

  /*!
    @brief How to interpret data locations. 
  */
  CellLocation m_dataLocation;

  /*!
    @brief Grids
  */
  EBLevelGrid m_eblg;

  /*!
    @brief Fine grids
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Refined coarse grids
  */
  EBLevelGrid m_eblgFiCo;

  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief True if we have a finer level (m_eblgFine is defined)
  */
  bool m_hasFine;

  /*!
    @brief True if we have EBCF crossings that require special stencils. 
  */
  bool m_hasEBCF;

  /*!
    @brief Resolution on this level. 
  */
  Real m_dx;

  /*!
    @brief Resolution on finer level. 
  */
  Real m_dxFine;

  /*!
    @brief Desired order for least squares stencils. 
  */
  int m_order;

  /*!
    @brief Weight (for least squares)
  */
  int m_weighting;

  /*!
    @brief Refinement factor between fine and coarse level. 
  */
  int m_refRat;

  /*!
    @brief Number of ghost cells in incoming data holders
  */
  IntVect m_ghostVector;

  /*!
    @brief Stencils for boundaries and cut-cells
    @note Lives on the coarse layout. 
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_levelStencils;

  /*!
    @brief Part of EBCF gradient stencil which reaches into the coarse level. 
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_ebcfStencilsCoar;

  /*!
    @brief Part of EBCF gradient stencil which reaches into the coarse level. 
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_ebcfStencilsFine;

  /*!
    @brief AggStencils for boundaries and cut-cells. Provides faster application of m_levelStencils.
    @note Lives on the coarse layout and only has one component. 
  */
  LayoutData<RefCountedPtr<AggStencil<EBCellFAB, EBCellFAB>>> m_aggLevelStencils[SpaceDim];

  /*!
    @brief Iterator for boundary cells (either domain or EB)
    @note Lives on the coarse layout. 
  */
  mutable LayoutData<VoFIterator> m_levelIterator;

  /*!
    @brief Iterator for EBCF cells. 
    @note Lives on the coarse layout. 
  */
  mutable LayoutData<VoFIterator> m_ebcfIterator;

  /*!
    @brief Define level stencils. 
  */
  virtual void
  defineLevelStencils() noexcept;

  /*!
    @brief Make masks
    @param[out] a_coarMaskCF      Mask holding a value of 1 in coarse cells that live on the refinement boundary. 
    @param[out] a_coarMaskInvalid Mask holding a value of 1 in coarse cells that are covered by a finer grid level. 
  */
  virtual void
  defineMasks(LevelData<FArrayBox>& a_coarMaskCF, LevelData<FArrayBox>& a_coarMaskInvalid) noexcept;

  /*!
    @brief Define iterators and cells that need modified stencils
    @param[in]  a_coarMaskCF            Mask for CF region on the coarse layout. 
    @param[in]  a_coarMaskInvalid       Mask for invalid regions on the coarse layout. 
  */
  virtual void
  defineIteratorsEBCF(const LevelData<FArrayBox>& a_coarMaskCF, const LevelData<FArrayBox>& a_coarMaskInvalid) noexcept;

  /*!
    @brief Define stencils for interpolation across the EBCF refinement boundary. 
    @param[in] a_bufferCoarMaskInvalid Map of coarse cells that are covered by finer cells. 
  */
  virtual void
  defineStencilsEBCF(const LevelData<FArrayBox>& a_coarMaskInvalid) noexcept;

  /*!
    @brief Get a finite difference stencil, assuming that the data is cell-centered
    @param[out] a_stencil       Finite difference stencil. 
    @param[in]  a_vof           Grid cell
    @param[in]  a_ebisBox       EBISBox
    @param[in]  a_invalidRegion Map of cells that can/can't be used. 
    @param[in]  a_dx            Grid resolution
  */
  virtual bool
  getFiniteDifferenceStencil(VoFStencil&            a_stencil,
                             const VolIndex&        a_vof,
                             const EBISBox&         a_ebisBox,
                             const DenseIntVectSet& a_invalidRegion,
                             const Real             a_dx) const noexcept;

  /*!
    @brief Build a gradient stencil using least squares. 
    @param[out] a_stencilCoar       Coarse-grid stencil
    @param[out] a_stencilFine       Fine-grid stencil
    @param[in]  a_vofCoar           VoF on the coarse level. 
    @param[in]  a_dataLocation      Data location
    @param[in]  a_ebisBoxCoar       Coarse grid EBISBox
    @param[in]  a_ebisBoxFine       Fine grid EBISBox
    @param[in]  a_validCellsCoar    Map of cells that are valid on the coarse domain. 
    @param[in]  a_validCellsFine    Map of cells that are valid on the coarse domain. 
    @param[in]  a_dxCoar            Grid resolution on the coarse level
    @param[in]  a_dxFine            Grid resolution on the fine level
    @param[in]  a_order             Interpolation order
    @param[in]  a_weight            Equation weighting for least squares reconstruction. 
    @return Returns true if a stencil was found, in which case it is defined in a_stencilFine and a_stencilCoar. 
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

  /*!
    @brief Check if a cell can produce valid finite difference stencils
    @details This will check for stencils for all vofs in a_ivCoar. 
    @param[in] a_ivCoar        Coarse grid cell. 
    @param[in] a_ebisBox       EBISBox
    @param[in] a_invalidRegion Map of cells that are covered by a finer cell. 
  */
  virtual bool
  isFiniteDifferenceStencilValid(const IntVect&   a_ivCoar,
                                 const EBISBox&   a_ebisBox,
                                 const FArrayBox& a_invalidRegion) const noexcept;

  /*!
    @brief Turn raw stencils into optimized AggStencils
  */
  virtual void
  makeAggStencils() noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
