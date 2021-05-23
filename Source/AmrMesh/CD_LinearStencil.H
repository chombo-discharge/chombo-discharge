/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LinearStencil.H
  @brief  Routines for computing a bi/tri-linear stencil in cut-cells, with modifications for a bunch of corner cases involving multi-cells. 
  @author Robert Marskar
*/

#ifndef CD_LinearStencil_H
#define CD_LinearStencil_H

// Chombo includes
#include <Stencil.H>
#include <ProblemDomain.H>
#include <EBISBox.H>

// Our includes
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Static class with routines for getting various kinds of stencils
*/
class LinearStencil {
public:

  static Real tolerance;

  /*!
    @brief Get a bi(tri)linear interpolation stencil in 2D(3D)
  */
  static bool getLinearInterpStencil(VoFStencil&          a_stencil, 
				     const RealVect&      a_centroid,
				     const VolIndex&      a_vof,
				     const ProblemDomain& a_domain,
				     const EBISBox&       a_ebisbox);

protected:
  
  /*!
    @brief Compute a bilinear interpolation stencil
  */
  static bool computeInterpStencil1D(VoFStencil&          a_stencil,
				     const RealVect&      a_centroid,
				     const VolIndex&      a_vof,
				     const ProblemDomain& a_domain,
				     const EBISBox&       a_ebisbox,
				     const int            a_interp_dir);
  /*!
    @brief Compute a bilinear interpolation stencil
  */
  static bool computeInterpStencil2D(VoFStencil&          a_stencil,
				     const RealVect&      a_centroid,
				     const VolIndex&      a_vof,
				     const ProblemDomain& a_domain,
				     const EBISBox&       a_ebisbox,
				     const int            a_no_interp_dir = 2);


  /*!
    @brief Compute a trilinear interpolation stencil
  */
  static bool computeInterpStencil3D(VoFStencil&          a_stencil,
				     const RealVect&      a_centroid,
				     const VolIndex&      a_vof,
				     const ProblemDomain& a_domain,
				     const EBISBox&       a_ebisbox);
};

#include <CD_NamespaceFooter.H>

#endif
