/*!
  @file centroid_interp.cpp
  @brief Implementation of centroid_interp.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "centroid_interp.H"
#include "stencil_ops.H"
#include "LeastSquares.H"

#include "EBArith.H"

Real centroid_interp::s_tolerance      = 1.E-8;
bool centroid_interp::s_quadrant_based = true;   // Use quadrant based stencils for least squares

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

  // Find the preferred stencil type
  if(m_stencil_type == stencil_type::linear){
    const RealVect centroid = a_ebisbox.centroid(a_vof);
    found_stencil = stencil_ops::get_linear_interp_stencil(a_sten, centroid, a_vof, a_domain, a_ebisbox);
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

  // If we couldn't find a stencil, just take the cell center to be the "interpolation". 
#if 0
  if(!found_stencil){
    const RealVect centroid = a_ebisbox.centroid(a_vof);
    found_stencil = stencil_ops::get_linear_interp_stencil(a_sten, centroid, a_vof, a_domain, a_ebisbox);
  }
  if(!found_stencil){
    found_stencil = this->get_taylor_stencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, a_cfivs);
  }
  if(!found_stencil){
    found_stencil = this->get_lsq_grad_stencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, a_cfivs);
  }
#endif

  if(!found_stencil){ // Drop to zeroth order. 
    a_sten.clear();
    a_sten.add(a_vof, 1.0);
  }

#if 1 // Debug
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

  return order > 0;
}

bool centroid_interp::get_lsq_grad_stencil(VoFStencil&              a_sten,
					   const VolIndex&          a_vof,
					   const DisjointBoxLayout& a_dbl,
					   const ProblemDomain&     a_domain,
					   const EBISBox&           a_ebisbox,
					   const Box&               a_box,
					   const Real&              a_dx,
					   const IntVectSet&        a_cfivs){
  CH_TIME("centroid_interp::get_lsq_grad_stencil");
#if 0 // Old way
#if CH_SPACEDIM == 2
  const int minStenSize = 3;
#elif CH_SPACEDIM== 3
  const int minStenSize = 7;
#endif

  a_sten.clear();
    
  Real weight              = 0;
  const int comp           = 0;
  const RealVect& centroid = a_ebisbox.centroid(a_vof);
  const RealVect normal    = centroid/centroid.vectorLength();

  // Find the quadrant that the normal points into
  IntVect quadrant;
  for (int dir = 0; dir < SpaceDim; dir++){
    quadrant[dir] = (normal[dir] < 0) ? -1 : 1;
  }

  // Quadrant or monotone-path based stencils
  if(s_quadrant_based){
    EBArith::getLeastSquaresGradSten(a_sten,               // Find a stencil for the gradient at the cell center
				     weight,
				     normal,
				     RealVect::Zero,
				     quadrant,
				     a_vof,
				     a_ebisbox,
				     a_dx*RealVect::Unit,
				     a_domain,
				     comp);
  }
  else{
    EBArith::getLeastSquaresGradStenAllVoFsRad(a_sten,               // Find a stencil for the gradient at the cell center
					       weight,
					       normal,
					       RealVect::Zero,
					       a_vof,
					       a_ebisbox,
					       a_dx*RealVect::Unit,
					       a_domain,
					       comp,
					       m_radius);
  }

  if(a_sten.size() >= minStenSize){         // Do a Taylor expansion phi_centroid = phi_cell + grad_cell*(x_centroid - x_cell)
    a_sten.add(a_vof, weight);
    a_sten *= m_dx*centroid.vectorLength();
    a_sten.add(a_vof, 1.0);
      
    return true;
  }
  else{
    return false;
  }
#else // New, and much better way.

  const int weightP = 0;
  const int radius  = 1;
  const int order   = 1;
  
  a_sten = LeastSquares::getInterpolationStencilUsingAllVoFsInRadius(LeastSquares::CellPosition::Centroid,
  								     LeastSquares::CellPosition::Center,
  								     a_vof,
  								     a_ebisbox,
  								     a_dx,
  								     weightP,
  								     radius,
  								     order);

#endif
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

  const Real thresh = 1.E-5;
  
  const RealVect centroid = a_ebisbox.centroid(a_vof);
  for (int dir = 0; dir < SpaceDim; dir++){

    bool hasLo, hasLower, hasHi, hasHigher;
    VolIndex loVoF, lowerVoF, hiVoF, higherVoF;
    EBArith::getVoFsDir(hasLo, loVoF, hasLower, lowerVoF,   a_ebisbox, a_vof, dir, Side::Lo, NULL);
    EBArith::getVoFsDir(hasHi, hiVoF, hasHigher, higherVoF, a_ebisbox, a_vof, dir, Side::Hi, NULL);

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
