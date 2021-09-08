/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbCentroidInterpolationStencil.H
  @brief  Declaration of IrregStencil-derived class that can from cell centers to embedded boundary centroids. 
  @author Robert Marskar
*/

#ifndef CD_EbCentroidInterpolationStencil_H
#define CD_EbCentroidInterpolationStencil_H

// Our includes
#include <CD_IrregStencil.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which can take cell-centered data and interpolate it to embedded boundary centroids. 
*/
class EbCentroidInterpolationStencil : public IrregStencil {
public:

  /*!
    @brief Disallowed weak constructor. 
  */
  EbCentroidInterpolationStencil() = delete;

  /*!
    @brief Full constructor. 
    @param[in] a_dbl    Grids
    @param[in] a_ebisl  EBIS layout
    @param[in] a_domain Problem domain
    @param[in] a_dx     Resolutions
    @param[in] a_order  Interpolation order
    @param[in] a_radius Radius for least squares
    @param[in] a_type   Stencil type
  */
  EbCentroidInterpolationStencil(const DisjointBoxLayout&        a_dbl,
				 const EBISLayout&               a_ebisl,
				 const ProblemDomain&            a_domain,
				 const Real&                     a_dx,
				 const int                       a_order,
				 const int                       a_radius,
				 const IrregStencil::StencilType a_type);

  /*!
    @brief Destructor (does nothing)
  */
  ~EbCentroidInterpolationStencil();

protected:
  
  /*!
    @brief Build the interpolation stencils
    @param[out] a_sten    Returned stencil
    @param[in]  a_vof     VolIndex where stencil is defined
    @param[in]  a_grids   Grids
    @param[in]  a_domain  Problem domain
    @param[in]  a_ebisbox EBISBox
    @param[in]  a_box     Box currently being evaluated
    @param[in]  a_dx      Grid resolution
    @param[in]  a_cfivs   Coarse-fine interface
    @details This calls one of the other function for computing the stencil. 
  */
  virtual void buildStencil(VoFStencil&              a_sten,
			    const VolIndex&          a_vof,
			    const DisjointBoxLayout& a_grids,
			    const ProblemDomain&     a_domain,
			    const EBISBox&           a_ebisbox,
			    const Box&               a_box,
			    const Real&              a_dx,
			    const IntVectSet&        a_cfivs) override;
  
  /*!
    @brief Build the interpolation stencils using a Taylor series expansion. 
    @param[out] a_sten    Returned stencil
    @param[in]  a_vof     VolIndex where stencil is defined
    @param[in]  a_grids   Grids
    @param[in]  a_domain  Problem domain
    @param[in]  a_ebisbox EBISBox
    @param[in]  a_box     Box currently being evaluated
    @param[in]  a_dx      Grid resolution
    @param[in]  a_cfivs   Coarse-fine interface
    @details This tries to compute up to second order derivatives for use in a Taylor series. 
  */
  bool getTaylorExtrapolationStencil(VoFStencil&              a_sten,
				     const VolIndex&          a_vof,
				     const DisjointBoxLayout& a_dbl,
				     const ProblemDomain&     a_domain,
				     const EBISBox&           a_ebisbox,
				     const Box&               a_box,
				     const Real&              a_dx,
				     const IntVectSet&        a_cfivs);

  /*!
    @brief Build the interpolation stencils using least squares interpolation.
    @param[out] a_sten    Returned stencil
    @param[in]  a_vof     VolIndex where stencil is defined
    @param[in]  a_grids   Grids
    @param[in]  a_domain  Problem domain
    @param[in]  a_ebisbox EBISBox
    @param[in]  a_box     Box currently being evaluated
    @param[in]  a_dx      Grid resolution
    @param[in]  a_cfivs   Coarse-fine interface
    @details This defines an overdetermined system of equations, using a large neighborhood of cells, which is solved with weighted least squarse. 
  */
  bool getLeastSquaresInterpolationStencil(VoFStencil&              a_sten,
					   const VolIndex&          a_vof,
					   const DisjointBoxLayout& a_dbl,
					   const ProblemDomain&     a_domain,
					   const EBISBox&           a_ebisbox,
					   const Box&               a_box,
					   const Real&              a_dx,
					   const IntVectSet&        a_cfivs);

  /*!
    @brief Build the interpolation stencils using piecewise linear reconstruction. 
    @param[out] a_sten    Returned stencil
    @param[in]  a_vof     VolIndex where stencil is defined
    @param[in]  a_grids   Grids
    @param[in]  a_domain  Problem domain
    @param[in]  a_ebisbox EBISBox
    @param[in]  a_box     Box currently being evaluated
    @param[in]  a_dx      Grid resolution
    @param[in]  a_cfivs   Coarse-fine interface
  */
  bool getPiecewiseLinearStencil(VoFStencil&              a_sten,
				 const VolIndex&          a_vof,
				 const DisjointBoxLayout& a_dbl,
				 const ProblemDomain&     a_domain,
				 const EBISBox&           a_ebisbox,
				 const Box&               a_box,
				 const Real&              a_dx,
				 const IntVectSet&        a_cfivs);
};

#include <CD_NamespaceFooter.H>

#endif
