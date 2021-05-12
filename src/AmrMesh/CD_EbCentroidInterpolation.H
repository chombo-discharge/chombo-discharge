/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   EbCentroidInterpolation.H
  @brief  Implementation of irreg_stencil that interpolates from cell centers to cell centroids
  @author Robert Marskar
*/

#ifndef CD_EbCentroidInterp_H
#define CD_EbCentroidInterp_H

// Our includes
#include <irreg_stencil.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which can take cell-centered data and interpolate it to embedded boundary centroids. 
*/
class EbCentroidInterpolation : public irreg_stencil {
public:

  /*!
    @brief Constructor. 
  */
  EbCentroidInterpolation() = delete;

  /*!
    @brief Full constructor. Subsequently calls define
    @param[in] a_dbl Grids
    @param[in] a_ebisl EBIS layout
    @param[in] a_domain Problem domain
    @param[in] a_dx Resolutions
    @param[in] a_order Interpolation order
    @param[in] a_radius Radius for least squares
    @param[in] a_type Stencil type
  */
  EbCentroidInterpolation(const DisjointBoxLayout& a_dbl,
			  const EBISLayout&        a_ebisl,
			  const ProblemDomain&     a_domain,
			  const Real&              a_dx,
			  const int                a_order,
			  const int                a_radius,
			  const stencil_type       a_type);

  /*!
    @brief Destructor
  */
  ~EbCentroidInterpolation();

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
  virtual void build_stencil(VoFStencil&              a_sten,
			     const VolIndex&          a_vof,
			     const DisjointBoxLayout& a_grids,
			     const ProblemDomain&     a_domain,
			     const EBISBox&           a_ebisbox,
			     const Box&               a_box,
			     const Real&              a_dx,
			     const IntVectSet&        a_cfivs) override;

private:
  
  /*!
    @brief Build the interpolation stencils using Taylor series. 
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
