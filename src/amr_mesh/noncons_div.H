/*!
  @file   noncons_div.H
  @brief  Implementation of irreg_stencil that can perform the nonconservative divergence averaging
  @author Robert Marskar
  @date   April 2020
*/


#ifndef _NONCONS_DIV_
#define _NONCONS_DIV_

#include "irreg_stencil.H"

/*!
  @brief Class for holding stencils on irregular cells over a single AMR level
*/
class noncons_div : public irreg_stencil {
public:

  /*!
    @brief Constructor
  */
  noncons_div();

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
  noncons_div(const DisjointBoxLayout&       a_dbl,
	      const EBISLayout&              a_ebisl,
	      const ProblemDomain&           a_domain,
	      const Real&                    a_dx,
	      const int                      a_order,
	      const int                      a_radius,
	      const stencil_type::which_type a_type);

  /*!
    @brief Destructor
  */
  ~noncons_div();
  
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
  
};
#endif
