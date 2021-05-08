/*!
  @file eb_centroid_interp.cpp
  @brief Implementation of eb_centroid_interp.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "eb_centroid_interp.H"
#include "LinearStencil.H"
#include "CD_LeastSquares.H"

#include "EBArith.H"

#include "CD_NamespaceHeader.H"

Real eb_centroid_interp::s_tolerance      = 1.E-8;
bool eb_centroid_interp::s_quadrant_based = true;   // Use quadrant based stencils for least squares

eb_centroid_interp::eb_centroid_interp() : irreg_stencil(){
  CH_TIME("eb_centroid_interp::eb_centroid_interp");
}

eb_centroid_interp::eb_centroid_interp(const DisjointBoxLayout&       a_dbl,
				       const EBISLayout&              a_ebisl,
				       const ProblemDomain&           a_domain,
				       const Real&                    a_dx,
				       const int                      a_order,
				       const int                      a_radius,
				       const stencil_type             a_type) : irreg_stencil(){

  CH_TIME("eb_centroid_interp::eb_centroid_interp");

  this->define(a_dbl, a_ebisl, a_domain, a_dx, a_order, a_radius, a_type);
}

eb_centroid_interp::~eb_centroid_interp(){
  CH_TIME("eb_centroid_interp::~eb_centroid_interp");
}

void eb_centroid_interp::build_stencil(VoFStencil&              a_sten,
				       const VolIndex&          a_vof,
				       const DisjointBoxLayout& a_dbl,
				       const ProblemDomain&     a_domain,
				       const EBISBox&           a_ebisbox,
				       const Box&               a_box,
				       const Real&              a_dx,
				       const IntVectSet&        a_cfivs){
  CH_TIME("eb_centroid_interp::build_stencil");

  
  bool found_stencil = false;

  // Find the preferred stencil type
  if(m_stencil_type == stencil_type::linear){
    const RealVect centroid = a_ebisbox.bndryCentroid(a_vof);
    found_stencil = LinearStencil::getLinearInterpStencil(a_sten, centroid, a_vof, a_domain, a_ebisbox);
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
    MayDay::Abort("eb_centroid_interp::build_stencil - Unsupported stencil type");
  }

  // If we couldn't find a stencil, try other types in this order
  if(!found_stencil){
    const RealVect centroid = a_ebisbox.bndryCentroid(a_vof);
    found_stencil = LinearStencil::getLinearInterpStencil(a_sten, centroid, a_vof, a_domain, a_ebisbox);
  }
  if(!found_stencil){
    found_stencil = this->get_taylor_stencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, a_cfivs);
  }
  if(!found_stencil){
    found_stencil = this->get_lsq_grad_stencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, a_cfivs);
  }

#if 0 // Only meant for testing but left in place: This code averages several neighbors over the EB. Only ment for testing so far
  found_stencil = this->get_ebavg_stencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, a_cfivs);
#endif

  if(!found_stencil){ // Drop to zeroth order.
    a_sten.clear();
    a_sten.add(a_vof, 1.0);
  }
}

bool eb_centroid_interp::get_taylor_stencil(VoFStencil&              a_sten,
					    const VolIndex&          a_vof,
					    const DisjointBoxLayout& a_dbl,
					    const ProblemDomain&     a_domain,
					    const EBISBox&           a_ebisbox,
					    const Box&               a_box,
					    const Real&              a_dx,
					    const IntVectSet&        a_cfivs){
  CH_TIME("eb_centroid_interp::get_taylor_stencil");
  
  const int comp           = 0;
  const RealVect& centroid = a_ebisbox.bndryCentroid(a_vof);
  IntVectSet* cfivs        = const_cast<IntVectSet*>(&a_cfivs);
  
  if(m_order == 1){
    EBArith::getFirstOrderExtrapolationStencil(a_sten, centroid*a_dx, a_dx*RealVect::Unit, a_vof, a_ebisbox, -1, cfivs, comp);
  }
  else if(m_order == 2){
    EBArith::getExtrapolationStencil(a_sten, centroid*a_dx, a_dx*RealVect::Unit, a_vof, a_ebisbox, -1, cfivs, comp);
  }
  else {
    MayDay::Abort("eb_centroid_interp::get_taylor_stencil - bad order requested. Only first and second order is supported");
  }

  return a_sten.size() > 0;

}

bool eb_centroid_interp::get_lsq_grad_stencil(VoFStencil&              a_sten,
					      const VolIndex&          a_vof,
					      const DisjointBoxLayout& a_dbl,
					      const ProblemDomain&     a_domain,
					      const EBISBox&           a_ebisbox,
					      const Box&               a_box,
					      const Real&              a_dx,
					      const IntVectSet&        a_cfivs){
  CH_TIME("eb_centroid_interp::get_lsq_grad_stencil");

  const int weightingPower = 0;
  
  a_sten = LeastSquares::getInterpolationStencilUsingAllVofsInRadius(LeastSquares::CellPosition::Boundary,
								     LeastSquares::CellPosition::Center,
								     a_vof,
								     a_ebisbox,
								     a_dx,
								     weightingPower,
								     m_radius,
								     m_order);

  return (a_sten.size() > 0);
}

bool eb_centroid_interp::get_pwl_stencil(VoFStencil&              a_sten,
					 const VolIndex&          a_vof,
					 const DisjointBoxLayout& a_dbl,
					 const ProblemDomain&     a_domain,
					 const EBISBox&           a_ebisbox,
					 const Box&               a_box,
					 const Real&              a_dx,
					 const IntVectSet&        a_cfivs){

  a_sten.clear();
  a_sten.add(a_vof, 1.0);

  const Real thresh = 1.E-5;
  
  const RealVect centroid = a_ebisbox.bndryCentroid(a_vof);
  for (int dir = 0; dir < SpaceDim; dir++){

    bool hasLo, hasLower, hasHi, hasHigher;
    VolIndex loVoF, lowerVoF, hiVoF, higherVoF;
    EBArith::getVoFsDir(hasLo, loVoF, hasLower, lowerVoF,   a_ebisbox, a_vof, dir, Side::Lo, (IntVectSet*) &a_cfivs);
    EBArith::getVoFsDir(hasHi, hiVoF, hasHigher, higherVoF, a_ebisbox, a_vof, dir, Side::Hi, (IntVectSet*) &a_cfivs);

    if(hasLo || hasHi){
      if(centroid[dir] > thresh && hasHi){ // Deriv is to the right
	a_sten.add(hiVoF,  centroid[dir]);
	a_sten.add(a_vof, -centroid[dir]);
      }
      else if(centroid[dir] < -thresh && hasLo){ // Deriv is to the left
	a_sten.add(loVoF, -centroid[dir]);
	a_sten.add(a_vof,  centroid[dir]);
      }
      else if(Abs(centroid[dir]) < thresh){
	// No deriv in this direction
      }
    }
  }

  return true;
}
#include "CD_NamespaceFooter.H"
