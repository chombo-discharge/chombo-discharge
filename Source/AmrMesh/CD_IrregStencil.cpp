/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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

constexpr int IrregStencil::m_defaultStenComp;

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

  LayoutData<IntVectSet> cfivs;
  EBArith::defineCFIVS(cfivs, a_dbl, a_domain);

  m_stencils.define(m_dbl);
  m_vofIter. define(m_dbl);
  for (DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit){
    const Box&     box     = m_dbl.get(dit());
    const EBISBox& ebisbox = m_ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs  = ebisbox.getIrregIVS(box);

    m_stencils[dit()] = RefCountedPtr<BaseIVFAB<VoFStencil> > (new BaseIVFAB<VoFStencil>(ivs, ebgraph, ncomp));


    VoFIterator& vofit = m_vofIter[dit()];
    vofit.define(ivs, ebgraph);
    
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      VoFStencil& stencil = (*m_stencils[dit()])(vof, 0);
      this->buildStencil(stencil, vof, m_dbl, m_domain, ebisbox, box, m_dx, cfivs[dit()]);

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

void IrregStencil::apply(EBCellFAB& a_dst, const EBCellFAB& a_src, const DataIndex& a_dit) const {
  CH_TIME("IrregStencil::apply");

  const BaseIVFAB<VoFStencil>& stencils = *m_stencils[a_dit];

  VoFIterator& vofit = m_vofIter[a_dit];
  
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex&   vof     = vofit();
    const VoFStencil& stencil = stencils(vof, m_defaultStenComp);

    for (int comp = 0; comp < a_dst.nComp(); comp++){

      a_dst(vof, comp) = 0.0;
      
      for (int i = 0; i < stencil.size(); i++){
	const VolIndex& ivof    = stencil.vof(i);
	const Real&     iweight = stencil.weight(i);

	a_dst(vof, comp) += iweight * a_src(ivof, comp);
      }
    }
  }
}

void IrregStencil::apply(BaseIVFAB<Real>& a_dst, const EBCellFAB& a_src, const DataIndex& a_dit) const {
  CH_TIME("IrregStencil::apply");

  const BaseIVFAB<VoFStencil>& stencils = *m_stencils[a_dit];

  VoFIterator& vofit = m_vofIter[a_dit];
  
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex&   vof     = vofit();
    const VoFStencil& stencil = stencils(vof, m_defaultStenComp);

    for (int comp = 0; comp < a_dst.nComp(); comp++){

      a_dst(vof, comp) = 0.0;
      
      for (int i = 0; i < stencil.size(); i++){
	const VolIndex& ivof    = stencil.vof(i);
	const Real&     iweight = stencil.weight(i);

	a_dst(vof, comp) += iweight * a_src(ivof, comp);
      }
    }
  }  
}

#include <CD_NamespaceFooter.H>
