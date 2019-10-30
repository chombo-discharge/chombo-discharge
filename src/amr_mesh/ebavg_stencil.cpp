/*!
  @file   ebavg_stencil.cpp
  @brief  Implementation of ebavg_stencil.H
  @author Robert Marskar
  @date   Nov. 2019
*/

#include "type_definitions.H"
#include "ebavg_stencil.H"
#include "stencil_ops.H"

#include "EBArith.H"

Real ebavg_stencil::s_tolerance      = 1.E-8;

ebavg_stencil::ebavg_stencil() : irreg_stencil(){
  CH_TIME("ebavg_stencil::ebavg_stencil");
}

ebavg_stencil::ebavg_stencil(const DisjointBoxLayout&       a_dbl,
			     const EBISLayout&              a_ebisl,
			     const ProblemDomain&           a_domain,
			     const Real&                    a_dx,
			     const int                      a_order,
			     const int                      a_radius,
			     const stencil_type::which_type a_type) : irreg_stencil(){

  CH_TIME("ebavg_stencil::ebavg_stencil");

  this->define(a_dbl, a_ebisl, a_domain, a_dx, a_order, a_radius, a_type);
}

ebavg_stencil::~ebavg_stencil(){
  CH_TIME("ebavg_stencil::~ebavg_stencil");
}

void ebavg_stencil::build_stencil(VoFStencil&              a_sten,
				  const VolIndex&          a_vof,
				  const DisjointBoxLayout& a_dbl,
				  const ProblemDomain&     a_domain,
				  const EBISBox&           a_ebisbox,
				  const Box&               a_box,
				  const Real&              a_dx,
				  const IntVectSet&        a_cfivs){
  CH_TIME("ebavg_stencil::build_stencil");

  bool found_stencil = false;

  // VolIndices
  Vector<VolIndex> otherVoFs;
  EBArith::getAllVoFsWithinRadius(otherVoFs, a_vof, a_ebisbox, m_radius);

  // Make sure the VoFs lies on the problemdomain and on the DBL. Don't reach across CFIVS
  Vector<VolIndex> vofList;
  for (int i = 0; i < otherVoFs.size(); i++){
    const IntVect iv = otherVoFs[i].gridIndex();

    
  }
}
