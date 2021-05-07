/*!
  @file   centroid_interp.cpp
  @brief  Implementation of centroid_interp.H
  @author Robert Marskar
  @date   Nov. 2017
*/

#include "centroid_interp.H"
#include "LinearStencil.H"
#include "CD_LeastSquares.H"

#include "EBArith.H"

#define DEBUG_CENTROID_INTERP 1

using namespace ChomboDischarge;

Real centroid_interp::s_tolerance = 1.E-8;

centroid_interp::centroid_interp() : irreg_stencil(){
  CH_TIME("centroid_interp::centroid_interp");
}

centroid_interp::centroid_interp(const DisjointBoxLayout&       a_dbl,
				 const EBISLayout&              a_ebisl,
				 const ProblemDomain&           a_domain,
				 const Real&                    a_dx,
				 const int                      a_order,
				 const int                      a_radius,
				 const stencil_type             a_type) : irreg_stencil(){

  CH_TIME("centroid_interp::centroid_interp");

  this->define(a_dbl, a_ebisl, a_domain, a_dx, a_order, a_radius, a_type);
}

centroid_interp::~centroid_interp(){
  CH_TIME("centroid_interp::~centroid_interp");
}

void centroid_interp::build_stencil(VoFStencil&              a_sten,
				    const VolIndex&          a_vof,
				    const DisjointBoxLayout& a_dbl,
				    const ProblemDomain&     a_domain,
				    const EBISBox&           a_ebisbox,
				    const Box&               a_box,
				    const Real&              a_dx,
				    const IntVectSet&        a_cfivs){
  CH_TIME("centroid_interp::build_stencil");
  
  bool found_stencil = false;

  if(m_stencil_type == stencil_type::linear){
    found_stencil = LinearStencil::getLinearInterpStencil(a_sten, a_ebisbox.centroid(a_vof), a_vof, a_domain, a_ebisbox);
  }
  else if(m_stencil_type == stencil_type::taylor){
    found_stencil = this->get_taylor_stencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, a_cfivs);
  }
  else if(m_stencil_type == stencil_type::lsq){
    found_stencil = this->get_lsq_grad_stencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, a_cfivs);
  }
  else if(m_stencil_type == stencil_type::pwl){
    found_stencil = this->get_pwl_stencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, a_cfivs);
  }
  else{
    MayDay::Abort("centroid_interp::build_stencil - Unsupported stencil type");
  }

  // Drop to zeroth order if we couldn't find a stencil. 
  if(!found_stencil){ 
    a_sten.clear();
    a_sten.add(a_vof, 1.0);
  }

#if DEBUG_CENTROID_INTERP
  for (int i = 0; i < a_sten.size(); i++){
    const Real w = a_sten.weight(i);
    if(w < 0.0) MayDay::Warning("centroid_interp::build_stencils - I got negative weights!!!");
  }
#endif
}

bool centroid_interp::get_taylor_stencil(VoFStencil&              a_sten,
					 const VolIndex&          a_vof,
					 const DisjointBoxLayout& a_dbl,
					 const ProblemDomain&     a_domain,
					 const EBISBox&           a_ebisbox,
					 const Box&               a_box,
					 const Real&              a_dx,
					 const IntVectSet&        a_cfivs){
  CH_TIME("centroid_interp::get_taylor_stencil");

  const Real tol           = 1.E-6;
  const int comp           = 0;
  const RealVect& centroid = a_ebisbox.centroid(a_vof);
  IntVectSet* cfivs        = const_cast<IntVectSet*>(&a_cfivs);

  int order;               
  if(m_order == 1){
    order = EBArith::getFirstOrderExtrapolationStencil(a_sten, centroid*a_dx, a_dx*RealVect::Unit, a_vof, a_ebisbox, -1, cfivs, comp);
  }
  else if(m_order == 2){
    order = EBArith::getExtrapolationStencil(a_sten, centroid*a_dx, a_dx*RealVect::Unit, a_vof, a_ebisbox, -1, cfivs, comp);
  }
  else {
    MayDay::Abort("centroid_interp::get_taylor_stencil - bad order requested. Only first and second order is supported");
  }

  return (order > 0);
}

bool centroid_interp::get_lsq_grad_stencil(VoFStencil&              a_sten,
					   const VolIndex&          a_vof,
					   const DisjointBoxLayout& a_dbl,
					   const ProblemDomain&     a_domain,
					   const EBISBox&           a_ebisbox,
					   const Box&               a_box,
					   const Real&              a_dx,
					   const IntVectSet&        a_cfivs){

  const int weightingPower = 0;
  
  a_sten = LeastSquares::getInterpolationStencilUsingAllVofsInRadius(LeastSquares::CellPosition::Centroid,
  								     LeastSquares::CellPosition::Center,
  								     a_vof,
  								     a_ebisbox,
  								     a_dx,
  								     weightingPower,
  								     m_radius,
  								     m_order);

  return (a_sten.size() > 0);
}

bool centroid_interp::get_pwl_stencil(VoFStencil&              a_sten,
				      const VolIndex&          a_vof,
				      const DisjointBoxLayout& a_dbl,
				      const ProblemDomain&     a_domain,
				      const EBISBox&           a_ebisbox,
				      const Box&               a_box,
				      const Real&              a_dx,
				      const IntVectSet&        a_cfivs){
  a_sten.clear();
  a_sten.add(a_vof, 1.0);
  
  const RealVect centroid = a_ebisbox.centroid(a_vof);
  for (int dir = 0; dir < SpaceDim; dir++){

    bool hasLo, hasLower, hasHi, hasHigher;
    VolIndex loVoF, lowerVoF, hiVoF, higherVoF;
    EBArith::getVoFsDir(hasLo, loVoF, hasLower, lowerVoF,   a_ebisbox, a_vof, dir, Side::Lo, (IntVectSet*) &a_cfivs);
    EBArith::getVoFsDir(hasHi, hiVoF, hasHigher, higherVoF, a_ebisbox, a_vof, dir, Side::Hi, (IntVectSet*) &a_cfivs);

    if(hasLo || hasHi){
      if(centroid[dir] > s_tolerance && hasHi){ // Deriv is to the right
	a_sten.add(hiVoF,  centroid[dir]);
	a_sten.add(a_vof, -centroid[dir]);
      }
      else if(centroid[dir] < -s_tolerance && hasLo){ // Deriv is to the left
	a_sten.add(loVoF, -centroid[dir]);
	a_sten.add(a_vof,  centroid[dir]);
      }
      else if(Abs(centroid[dir]) < s_tolerance){
	// No deriv in this direction
      }
    }
  }


  return true;
}
