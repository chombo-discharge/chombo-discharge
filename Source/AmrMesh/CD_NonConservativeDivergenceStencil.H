/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   NonConservativeDivergenceStencil.H
  @brief  Implementation of IrregStencil that can perform the nonconservative divergence averaging
  @author Robert Marskar
*/

#ifndef CD_NonConservativeDivergenceStencil_H
#define CD_NonConservativeDivergenceStencil_H

// Our includes
#include <CD_IrregStencil.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for holding stencils on irregular cells over a single AMR level
*/
class NonConservativeDivergenceStencil : public IrregStencil {
public:

  /*!
    @brief Disallowed weak constructor. 
  */
  NonConservativeDivergenceStencil() = delete;

  /*!
    @brief Full constructor. Subsequently calls define
    @param[in] a_dbl    Grids
    @param[in] a_ebisl  EBIS layout
    @param[in] a_domain Problem domain
    @param[in] a_dx     Resolutions
    @param[in] a_order  Interpolation order
    @param[in] a_radius Radius for least squares
    @param[in] a_type   Stencil type
  */
  NonConservativeDivergenceStencil(const DisjointBoxLayout&        a_dbl,
				   const EBISLayout&               a_ebisl,
				   const ProblemDomain&            a_domain,
				   const Real&                     a_dx,
				   const int                       a_order,
				   const int                       a_radius,
				   const IrregStencil::StencilType a_type);

  /*!
    @brief Destructor
  */
  ~NonConservativeDivergenceStencil();
  
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
			    const DisjointBoxLayout& a_dbl,
			    const ProblemDomain&     a_domain,
			    const EBISBox&           a_ebisbox,
			    const Box&               a_box,
			    const Real&              a_dx,
			    const IntVectSet&        a_cfivs);
  
};
  
#include <CD_NamespaceFooter.H>

#endif
