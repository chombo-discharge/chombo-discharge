/*!
  @file   noncons_div.cpp
  @brief  Implementation of noncons_div.H
  @author Robert Marskar
  @date   May 2020
*/

#include "noncons_div.H"
#include "CD_VofUtils.H"

#include "EBArith.H"

#include "CD_NamespaceHeader.H"

noncons_div::noncons_div() : IrregStencil(){
  CH_TIME("noncons_div::noncons_div");
}

noncons_div::noncons_div(const DisjointBoxLayout&       a_dbl,
			 const EBISLayout&              a_ebisl,
			 const ProblemDomain&           a_domain,
			 const Real&                    a_dx,
			 const int                      a_order,
			 const int                      a_radius,
			 const IrregStencil::StencilType             a_type) : IrregStencil() { 

  CH_TIME("noncons_div::noncons_div");

  // Order and radius are dummy arguments. 
  this->define(a_dbl, a_ebisl, a_domain, a_dx, a_order, a_radius, IrregStencil::StencilType::Linear);
}

noncons_div::~noncons_div(){
  CH_TIME("noncons_div::~noncons_div");
}

void noncons_div::buildStencil(VoFStencil&              a_sten,
				const VolIndex&          a_vof,
				const DisjointBoxLayout& a_dbl,
				const ProblemDomain&     a_domain,
				const EBISBox&           a_ebisbox,
				const Box&               a_box,
				const Real&              a_dx,
				const IntVectSet&        a_cfivs){
  CH_TIME("noncons_div::buildStencil");
  
  a_sten.clear();

  Real norm = 0.;

  const Vector<VolIndex> vofs = VofUtils::getAllConnectedVofsInRadius(a_vof, a_ebisbox, m_radius, IntVectSet());
  for (int i = 0; i < vofs.size(); i++){
    if(vofs[i] != a_vof){
      const VolIndex& ivof = vofs[i];
      const Real iweight   = a_ebisbox.volFrac(ivof);

      norm += iweight;
      a_sten.add(ivof, iweight);
    }
  }

  a_sten *= 1./norm;
}
#include "CD_NamespaceFooter.H"
