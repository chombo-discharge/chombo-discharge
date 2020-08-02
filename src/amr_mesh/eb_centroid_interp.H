/*!
  @file eb_centroid_interp.H
  @brief Implementation of irreg_stencil that interpolates from cell centers to cell centroids
  @author Robert Marskar
  @date Aug. 2017
*/


#ifndef _EB_CENTROID_INTERP_
#define _EB_CENTROID_INTERP_

#include "type_definitions.H"
#include "irreg_stencil.H"

/*!
  @brief Class for holding stencils on irregular cells over a single AMR level
*/
class eb_centroid_interp : public irreg_stencil {
public:

  /*!
    @brief Constructor
  */
  eb_centroid_interp();

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
  eb_centroid_interp(const DisjointBoxLayout&       a_dbl,
		     const EBISLayout&              a_ebisl,
		     const ProblemDomain&           a_domain,
		     const Real&                    a_dx,
		     const int                      a_order,
		     const int                      a_radius,
		     const stencil_type             a_type);

  /*!
    @brief Destructor
  */
  ~eb_centroid_interp();

  /*!
    @brief Tolerance
  */
  static Real s_tolerance;

  /*!
    @brief Use quadrant based stencils or not
  */
  static bool s_quadrant_based;

protected:
  
  /*!
    @brief Build the interpolation stencils
  */
  virtual void build_stencil(VoFStencil&              a_sten,
			     const VolIndex&          a_vof,
			     const DisjointBoxLayout& a_dbl,
			     const ProblemDomain&     a_domain,
			     const EBISBox&           a_ebisbox,
			     const Box&               a_box,
			     const Real&              a_dx,
			     const IntVectSet&        a_cfivs);

  /*!
    @brief Get a Taylor-based stencil
  */
  virtual bool get_taylor_stencil(VoFStencil&              a_sten,
				  const VolIndex&          a_vof,
				  const DisjointBoxLayout& a_dbl,
				  const ProblemDomain&     a_domain,
				  const EBISBox&           a_ebisbox,
				  const Box&               a_box,
				  const Real&              a_dx,
				  const IntVectSet&        a_cfivs);

  /*!
    @brief Get a least-squares based first order stencil
  */
  virtual bool get_lsq_grad_stencil(VoFStencil&              a_sten,
				    const VolIndex&          a_vof,
				    const DisjointBoxLayout& a_dbl,
				    const ProblemDomain&     a_domain,
				    const EBISBox&           a_ebisbox,
				    const Box&               a_box,
				    const Real&              a_dx,
				    const IntVectSet&        a_cfivs);

  /*!
    @brief Get a piecewise linear stencil based on the cell-centered derivative
  */
  virtual bool get_pwl_stencil(VoFStencil&              a_sten,
			       const VolIndex&          a_vof,
			       const DisjointBoxLayout& a_dbl,
			       const ProblemDomain&     a_domain,
			       const EBISBox&           a_ebisbox,
			       const Box&               a_box,
			       const Real&              a_dx,
			       const IntVectSet&        a_cfivs);

  /*!
    @brief Get the average stencil
  */
  virtual bool get_ebavg_stencil(VoFStencil&              a_sten,
				 const VolIndex&          a_vof,
				 const DisjointBoxLayout& a_dbl,
				 const ProblemDomain&     a_domain,
				 const EBISBox&           a_ebisbox,
				 const Box&               a_box,
				 const Real&              a_dx,
				 const IntVectSet&        a_cfivs);
  
};
#endif
