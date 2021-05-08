/*!
  @file irreg_stencil.cpp
  @brief Implementation of irreg_stencil.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "irreg_stencil.H"

#include <EBArith.H>

#include "CD_NamespaceHeader.H"

irreg_stencil::irreg_stencil(){
  CH_TIME("irreg_stencil::irreg_stencil");
  m_defined = false;
}

irreg_stencil::~irreg_stencil(){

}

irreg_stencil::irreg_stencil(const DisjointBoxLayout&       a_dbl,
			     const EBISLayout&              a_ebisl,
			     const ProblemDomain&           a_domain,
			     const Real&                    a_dx,
			     const int                      a_order,
			     const int                      a_radius,
			     const stencil_type             a_type){

  this->define(a_dbl, a_ebisl, a_domain, a_dx, a_order, a_radius, a_type);
}


void irreg_stencil::define(const DisjointBoxLayout&       a_dbl,
			   const EBISLayout&              a_ebisl,
			   const ProblemDomain&           a_domain,
			   const Real&                    a_dx,
			   const int                      a_order,
			   const int                      a_radius,
			   const stencil_type             a_type){

  m_dbl          = a_dbl;
  m_ebisl        = a_ebisl;
  m_dx           = a_dx;
  m_radius       = a_radius;
  m_order        = a_order;
  m_stencil_type = a_type;
  
  const int ncomp = 1;

  // LayoutData<IntVectSet> cfivs;
  // EBArith::defineCFIVS(cfivs, m_dbl, m_domain);

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
      this->build_stencil(stencil, vof, m_dbl, m_domain, ebisbox, box, m_dx, IntVectSet());// Can respect the CFIVS if we want to, but we dont. cfivs[dit()]);

#if 0 // Safety test
      Real sum = 0.0;
      for (int i = 0; i < stencil.size(); i++){
	sum += stencil.weight(i);
      }

      if(Abs(sum - 1.0) > 1.E-6){
	MayDay::Warning("irreg_stencil::define - weights do not sum to 1. Something has gone wrong with one of your stencils");
	stencil *= 1./sum;
      }
#endif
    }
  }  

  m_defined = true;
}
#include "CD_NamespaceFooter.H"
