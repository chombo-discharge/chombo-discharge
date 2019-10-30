/*!
  @file   ebavg_stencil.H
  @brief  Declaration of stencils that sets an EB value to be the weighted average of the connected EB cells
  @author Robert Marskar
  @date   Nov. 2019
*/


#ifndef _EBAVG_STENCIL_
#define _EBAVG_STENCIL_

#include "type_definitions.H"
#include "irreg_stencil.H"

/*!
  @brief Class for holding stencils on irregular cells over a single AMR level
*/
class ebavg_stencil : public irreg_stencil {
public:

  /*!
    @brief Constructor
  */
  ebavg_stencil();

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
  ebavg_stencil(const DisjointBoxLayout&       a_dbl,
		const EBISLayout&              a_ebisl,
		const ProblemDomain&           a_domain,
		const Real&                    a_dx,
		const int                      a_order,
		const int                      a_radius,
		const stencil_type::which_type a_type);

  /*!
    @brief Destructor
  */
  ~ebavg_stencil();

  /*!
    @brief Tolerance
  */
  static Real s_tolerance;

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
