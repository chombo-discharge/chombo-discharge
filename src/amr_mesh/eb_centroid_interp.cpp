/*!
  @file eb_centroid_interp.cpp
  @brief Implementation of eb_centroid_interp.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "eb_centroid_interp.H"
#include "stencil_ops.H"

#include "EBArith.H"

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
    MayDay::Abort("eb_centroid_interp::build_stencil - Unsupported stencil type");
  }

  // If we couldn't find a stencil, try other types in this order
  if(!found_stencil){
    const RealVect centroid = a_ebisbox.bndryCentroid(a_vof);
    found_stencil = stencil_ops::get_linear_interp_stencil(a_sten, centroid, a_vof, a_domain, a_ebisbox);
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

#if CH_SPACEDIM == 2
  const int minStenSize = 3;
#elif CH_SPACEDIM== 3
  const int minStenSize = 7;
#endif

  a_sten.clear();
    
  Real weight              = 0;
  const int comp           = 0;
  const RealVect& centroid = a_ebisbox.bndryCentroid(a_vof);
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



bool eb_centroid_interp::get_ebavg_stencil(VoFStencil&              a_sten,
					   const VolIndex&          a_vof,
					   const DisjointBoxLayout& a_dbl,
					   const ProblemDomain&     a_domain,
					   const EBISBox&           a_ebisbox,
					   const Box&               a_box,
					   const Real&              a_dx,
					   const IntVectSet&        a_cfivs){

  const int avg_radius = 1;

  a_sten.clear();

  // Get all the other VoFs
  Vector<VolIndex> candidateVoFs;
  EBArith::getAllVoFsInMonotonePath(candidateVoFs, a_vof, a_ebisbox, avg_radius);
  //  std::cout << candidateVoFs.size() << std::endl;
  const IntVect IV = a_vof.gridIndex();

  // Make sure the other VoFs lie on the domain
  Vector<VolIndex> goodVoFs;
  for (int i = 0; i < candidateVoFs.size(); i++){
    const VolIndex vof = candidateVoFs[i];
    const IntVect iv   = vof.gridIndex();

    if(vof != a_vof){
      if(a_ebisbox.isIrregular(iv) && !a_cfivs.contains(iv)){
	goodVoFs.push_back(vof);
      }
    }
  }
  //  goodVoFs.push_back(a_vof);

  const RealVect c0 = a_ebisbox.bndryCentroid(a_vof);

  ///  std::cout << goodVoFs.size() << std::endl;
  // Build the basic stencil
  a_sten.clear();
  Real sumweights = 0.0;

  const bool found_stencil = this->get_pwl_stencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, a_cfivs);
  const Real w0 = a_ebisbox.bndryArea(a_vof);
  a_sten *= w0;
  sumweights += w0;

  // Now do the other stencils
  for (int i = 0; i < goodVoFs.size(); i++){
    VoFStencil curSten;
    VolIndex curVoF = goodVoFs[i];
    const bool found_stencil = this->get_pwl_stencil(curSten, curVoF, a_dbl, a_domain, a_ebisbox, a_box, a_dx, a_cfivs);

    if(found_stencil){
      Real weight = a_ebisbox.bndryArea(curVoF);
      RealVect c = a_ebisbox.bndryCentroid(curVoF);

      RealVect d = (RealVect(IV) + c0) - (RealVect(curVoF.gridIndex())+c);
      //      weight *= 1./(1+d.vectorLength());
      
      curSten *= weight;

      sumweights += weight;
      a_sten += curSten;
    }
  }

  a_sten *= 1./sumweights;


  return true;
}
