/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NonConservativeDivergenceStencil.H
  @brief  Implementation of IrregStencil that can perform the nonconservative divergence averaging
  @author Robert Marskar
*/

#ifndef CD_NonConservativeDivergenceStencil_H
#define CD_NonConservativeDivergenceStencil_H

// Our includes
#include <CD_IrregStencil.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Stencil class for computing non-conservative divergences using neighborhood algorithms. 
  @details This class will compute a stencil s(phi) = sum(kappa * phi)/sum(kappa) where
  the stencil points consist of cells that can be reached through a monotone path. 
*/
class NonConservativeDivergenceStencil : public IrregStencil
{
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
    @param[in] a_order  Stencil order (dummy argument)
    @param[in] a_radius Stencil radius
    @param[in] a_type   Stencil type (dummy argument)
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
  virtual ~NonConservativeDivergenceStencil();

protected:
  /*!
    @brief Build the interpolation stencils
    @param[out] a_sten    Returned stencil
    @param[in]  a_vof     VolIndex where stencil is defined
    @param[in]  a_dbl     Grids
    @param[in]  a_domain  Problem domain
    @param[in]  a_ebisbox EBISBox
    @param[in]  a_box     Box currently being evaluated
    @param[in]  a_dx      Grid resolution
    @param[in]  a_cfivs   Coarse-fine interface
    @details This calls one of the other function for computing the stencil. 
  */
  virtual void
  buildStencil(VoFStencil&              a_sten,
               const VolIndex&          a_vof,
               const DisjointBoxLayout& a_dbl,
               const ProblemDomain&     a_domain,
               const EBISBox&           a_ebisbox,
               const Box&               a_box,
               const Real&              a_dx,
               const IntVectSet&        a_cfivs) override;
};

#include <CD_NamespaceFooter.H>

#endif
