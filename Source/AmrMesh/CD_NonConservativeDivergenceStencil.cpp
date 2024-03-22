/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NonConservativeDivergenceStencil.cpp
  @brief  Implementation of CD_NonConservativeDivergenceStencil.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_NonConservativeDivergenceStencil.H>
#include <CD_VofUtils.H>
#include <CD_NamespaceHeader.H>

NonConservativeDivergenceStencil::NonConservativeDivergenceStencil(const DisjointBoxLayout&        a_dbl,
                                                                   const EBISLayout&               a_ebisl,
                                                                   const ProblemDomain&            a_domain,
                                                                   const Real&                     a_dx,
                                                                   const int                       a_order,
                                                                   const int                       a_radius,
                                                                   const IrregStencil::StencilType a_type)
  : IrregStencil()
{

  CH_TIME("NonConservativeDivergenceStencil::NonConservativeDivergenceStencil");

  // Order and radius are dummy arguments.
  this->define(a_dbl, a_ebisl, a_domain, a_dx, a_order, a_radius, IrregStencil::StencilType::Linear);
}

NonConservativeDivergenceStencil::~NonConservativeDivergenceStencil()
{
  CH_TIME("NonConservativeDivergenceStencil::~NonConservativeDivergenceStencil");
}

void
NonConservativeDivergenceStencil::buildStencil(VoFStencil&              a_sten,
                                               const VolIndex&          a_vof,
                                               const DisjointBoxLayout& a_dbl,
                                               const ProblemDomain&     a_domain,
                                               const EBISBox&           a_ebisbox,
                                               const Box&               a_box,
                                               const Real&              a_dx,
                                               const IntVectSet&        a_cfivs)
{
  CH_TIME("NonConservativeDivergenceStencil::buildStencil");

  a_sten.clear();

  Real sumKappa = 0.;

  const Vector<VolIndex> vofs = VofUtils::getVofsInRadius(a_vof,
                                                          a_ebisbox,
                                                          m_radius,
                                                          VofUtils::Connectivity::MonotonePath,
                                                          true);

  for (int i = 0; i < vofs.size(); i++) {
    const VolIndex& ivof  = vofs[i];
    const Real&     kappa = a_ebisbox.volFrac(ivof);

    sumKappa += kappa;

    a_sten.add(ivof, kappa);
  }

  a_sten *= 1. / sumKappa;
}

#include <CD_NamespaceFooter.H>
