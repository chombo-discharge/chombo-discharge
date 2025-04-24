/* chombo-discharge
 * Copyright © 2025 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CellCentroidInterpolation.H
  @brief  Declaration of class for interpolating cell-centered data to cell centroids. 
  @author Robert Marskar
*/

#ifndef CD_CellCentroidInterpolation_H
#define CD_CellCentroidInterpolation_H

// Chombo includes
#include <EBLevelGrid.H>
#include <LevelData.H>
#include <LayoutData.H>
#include <VoFIterator.H>
#include <Stencils.H>
#include <BaseIVFAB.H>
#include <EBCellFAB.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for interpolating cell-centered data to the cell centroid.

  @details Various interpolation methods are supported, where some are supported by explicitly
  expressable stencils and others through limiters. 

  This is a per-level class without AMR functionality. On coarse-fine interfaces near the EB the
  user must fill the ghost cells prior to the interpolation. 
*/
class CellCentroidInterpolation
{
public:
  /*!
    @brief Supported interpolation types
  */
  enum class Type
  {
    Constant,           // Uses cell-centered value
    Linear,             // Uses bi/tri-linear interpolation
    Taylor,             // Uses Taylor-extrapolation
    LeastSquares,       // Uses least-squares reconstruction
    PiecewiseLinear,    // Uses piecewise-linear interpolation
    MinMod,             // Uses slope limiter with minmod slope
    MonotonizedCentral, // Uses slope limiter with van Leer limiter
    Superbee            // Uses slope limiter with superbee limiter
  };

  /*!
    @brief Default constructor. Must call the define function afterwards
  */
  CellCentroidInterpolation() noexcept;

  /*!
    @brief Copy constructor. Not allowed.
  */
  CellCentroidInterpolation(const CellCentroidInterpolation& a_other) = delete;

  /*!
    @brief Defining constructor. This calls the define function.
    @param[in] a_eblg Grids
    @param[in] a_dx Grid resolution
    @param[in] a_interpolationType Interpolation type
  */
  CellCentroidInterpolation(const EBLevelGrid& a_eblg, const Real& a_dx, const Type& a_interpolationType) noexcept;

  /*!
    @brief Destructor. Does nothing.
  */
  virtual ~CellCentroidInterpolation() noexcept;

  /*!
    @brief Define function. Puts object in usable state. 
    @param[in] a_eblg Grids
    @param[in] a_dx Grid resolution
    @param[in] a_interpolationType Interpolation type
  */
  virtual void
  define(const EBLevelGrid& a_eblg, const Real& a_dx, const Type& a_interpolationType) noexcept;

  /*!
    @brief Function for interpolating data. 
    @param[out] a_centroidData Centroid-centered data.
    @param[in] a_cellData Cell-centered data.
  */
  virtual void
  interpolate(LevelData<BaseIVFAB<Real>>& a_centroidData, const LevelData<EBCellFAB>& a_cellData) const noexcept;

  /*!
    @brief Function for interpolating data.
    @details This version will copy the data from the regular cells and interpolate data in the cut-cells.
    @param[out] a_centroidData Centroid-centered data.
    @param[in] a_cellData Cell-centered data.
  */
  virtual void
  interpolate(LevelData<EBCellFAB>& a_centroidData, const LevelData<EBCellFAB>& a_cellData) const noexcept;

  /*!
    @brief Function for interpolating data.
    @param[inout] a_centroidData Cell-centered on input and centroid centered on output.
  */
  virtual void
  interpolate(LevelData<EBCellFAB>& a_data) const noexcept;

  /*!
    @brief Interpolation function. 
    @param[out] a_centroidData Centroid-centered data.
    @param[in] a_cellData Cell-centered data.
    @note T can be either a BaseIVFAB<Real> or an EBCellFAB.
  */
  template <typename T>
  void
  interpolate(T& a_centroidData, const EBCellFAB& a_cellData, const DataIndex& a_din) const noexcept;

protected:
  /*!
    @brief Defined or not
  */
  bool m_isDefined;

  /*!
    @brief Grids
  */
  EBLevelGrid m_eblg;

  /*!
    @brief Interpolation type
  */
  Type m_interpolationType;

  /*!
    @brief Grid resolution
  */
  Real m_dx;

  /*!
    @brief VoF iterator for iterating through the cut-cells. 
  */
  mutable LayoutData<VoFIterator> m_vofIterator;

  /*!
    @brief Stencils. This is only populated when the interpolation is expressable as a stencil.
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_interpStencils;

  /*!
    @brief Utility function for fetching a bilinear/trilinear stencil
    @param[out] a_stencil Stencil
    @param[in] a_vof VolIndex for cut-cell
    @param[in] a_ebisBox EB grid box
    @param[in] a_domain Grid domain
  */
  virtual bool
  getLinearStencil(VoFStencil&          a_stencil,
                   const VolIndex&      a_vof,
                   const EBISBox&       a_ebisBox,
                   const ProblemDomain& a_domain) const noexcept;

  /*!
    @brief Utility function for fetching a Taylor-extrapolation stencil.
    @param[out] a_stencil Stencil
    @param[in] a_vof VolIndex for cut-cell
    @param[in] a_ebisBox EB grid box
    @param[in] a_domain Grid domain    
  */
  virtual bool
  getTaylorExtrapolationStencil(VoFStencil&          a_stencil,
                                const VolIndex&      a_vof,
                                const EBISBox&       a_ebisBox,
                                const ProblemDomain& a_domain) const noexcept;

  /*!
    @brief Utility function for fetching a least-squares based interpolation stencil.
    @param[out] a_stencil Stencil
    @param[in] a_vof VolIndex for cut-cell
    @param[in] a_ebisBox EB grid box
    @param[in] a_domain Grid domain    
  */
  virtual bool
  getLeastSquaresStencil(VoFStencil&          a_stencil,
                         const VolIndex&      a_vof,
                         const EBISBox&       a_ebisBox,
                         const ProblemDomain& a_domain) const noexcept;

  /*!
    @brief Utility function for fetching a piecewise linear stencil
    @param[out] a_stencil Stencil
    @param[in] a_vof VolIndex for cut-cell
    @param[in] a_ebisBox EB grid box
    @param[in] a_domain Grid domain    
  */
  virtual bool
  getPiecewiseLinearStencil(VoFStencil&          a_stencil,
                            const VolIndex&      a_vof,
                            const EBISBox&       a_ebisBox,
                            const ProblemDomain& a_domain) const noexcept;

  /*!
    @brief Minmod slope function. 
    @param[in] a_dwl Left slope
    @param[in] a_dwr Right slope
  */
  inline Real
  MinMod(const Real& a_dwl, const Real& a_dwr) const noexcept;

  /*!
    @brief Monotonized central difference slope limiter.
    @param[in] a_dwl Left slope
    @param[in] a_dwr Right slope
  */
  inline Real
  MonotonizedCentral(const Real& a_dwl, const Real& a_dwr) const noexcept;

  /*!
    @brief Superbee slope limiter. 
    @param[in] a_dwl Left slope
    @param[in] a_dwr Right slope
  */
  inline Real
  Superbee(const Real& a_dwl, const Real& a_dwr) const noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_CellCentroidInterpolationImplem.H>

#endif
