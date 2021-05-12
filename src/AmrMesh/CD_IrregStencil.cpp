/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   IrregStencil.cpp
  @brief  Implementation of IrregStencil.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_IrregStencil.H>
#include <CD_NamespaceHeader.H>

IrregStencil::IrregStencil(){
  CH_TIME("IrregStencil::IrregStencil");
  m_defined = false;
}

IrregStencil::~IrregStencil(){

}

IrregStencil::IrregStencil(const DisjointBoxLayout&          a_dbl,
			     const EBISLayout&               a_ebisl,
			     const ProblemDomain&            a_domain,
			     const Real&                     a_dx,
			     const int                       a_order,
			     const int                       a_radius,
			     const IrregStencil::StencilType a_type){

  this->define(a_dbl, a_ebisl, a_domain, a_dx, a_order, a_radius, a_type);
}


void IrregStencil::define(const DisjointBoxLayout&         a_dbl,
			   const EBISLayout&               a_ebisl,
			   const ProblemDomain&            a_domain,
			   const Real&                     a_dx,
			   const int                       a_order,
			   const int                       a_radius,
			   const IrregStencil::StencilType a_type){

  m_dbl          = a_dbl;
  m_ebisl        = a_ebisl;
  m_dx           = a_dx;
  m_radius       = a_radius;
  m_order        = a_order;
  m_stencilType = a_type;
  
  const int ncomp = 1;

  m_stencils.define(m_dbl);
  for (DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit){
    const Box&     box     = m_dbl.get(dit());
    const EBISBox& ebisbox = m_ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs  = ebisbox.getIrregIVS(box);

    m_stencils[dit()] = RefCountedPtr<BaseIVFAB<VoFStencil> > (new BaseIVFAB<VoFStencil>(ivs, ebgraph, ncomp));

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      VoFStencil& stencil = (*m_stencils[dit()])(vof, 0);
      this->buildStencil(stencil, vof, m_dbl, m_domain, ebisbox, box, m_dx, IntVectSet());// Can respect the CFIVS if we want to, but we dont. cfivs[dit()]);

#if 0 // Safety test
      Real sum = 0.0;
      for (int i = 0; i < stencil.size(); i++){
	sum += stencil.weight(i);
      }

      if(Abs(sum - 1.0) > 1.E-6){
	MayDay::Warning("IrregStencil::define - weights do not sum to 1. Something has gone wrong with one of your stencils");
	stencil *= 1./sum;
      }
#endif
    }
  }  

}

#include <CD_NamespaceFooter.H>
