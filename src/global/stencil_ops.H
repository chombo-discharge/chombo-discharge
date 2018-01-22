/*!
  @file stencil_ops.H
  @brief Agglomeration of stencil stuff
  @author Robert Marskar
  @date Jan. 2018
*/

#ifndef _STENCIL_OPS_
#define _STENCIL_OPS_

#include <Stencil.H>
#include <ProblemDomain.H>
#include <EBISBox.H>

/*!
  @brief Static class with routines for getting various kinds of stencils
*/
class stencil_ops {
public:

  static Real tolerance;

  /*!
    @brief Get a bi(tri)linear interpolation stencil in 2D(3D)
  */
  static bool get_linear_interp_stencil(VoFStencil&          a_stencil, 
					const RealVect&      a_centroid,
					const VolIndex&      a_vof,
					const ProblemDomain& a_domain,
					const EBISBox&       a_ebisbox);

  /*!
    @brief Compute a bilinear interpolation stencil
  */
  static bool compute_bilinear_stencil(VoFStencil&          a_stencil,
				       const RealVect&      a_centroid,
				       const VolIndex&      a_vof,
				       const ProblemDomain& a_domain,
				       const EBISBox&       a_ebisbox);

  /*!
    @brief Compute a bilinear interpolation stencil
  */
  static bool compute_interp_stencil_1D(VoFStencil&          a_stencil,
					const RealVect&      a_centroid,
					const VolIndex&      a_vof,
					const ProblemDomain& a_domain,
					const EBISBox&       a_ebisbox,
					const int            a_interp_dir);
  /*!
    @brief Compute a bilinear interpolation stencil
  */
  static bool compute_interp_stencil_2D(VoFStencil&          a_stencil,
					const RealVect&      a_centroid,
					const VolIndex&      a_vof,
					const ProblemDomain& a_domain,
					const EBISBox&       a_ebisbox,
					const int            a_no_interp_dir = 2);


  /*!
    @brief Compute a trilinear interpolation stencil
  */
  static bool compute_interp_stencil_3D(VoFStencil&          a_stencil,
					const RealVect&      a_centroid,
					const VolIndex&      a_vof,
					const ProblemDomain& a_domain,
					const EBISBox&       a_ebisbox);
};

#endif
