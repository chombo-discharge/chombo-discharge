/*!
  @file   robinconductivitydomainbc.cpp
  @brief  Implementation of robinconductivitydomainbc.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "robinconductivitydomainbc.H"
#include "EBArithF_F.H"
#include "DirichletPoissonDomainBCF_F.H"
#include "PolyGeom.H"

robinconductivitydomainbc::robinconductivitydomainbc(){
  this->set_coefs(1., 1., 0);
}

robinconductivitydomainbc::~robinconductivitydomainbc(){
}

void robinconductivitydomainbc::set_coefs(RefCountedPtr<robin_coef> a_robinco){
  m_robinco = a_robinco;

  m_const_coeff = false;
  m_func_coeff  = true;
  m_data_coeff  = false;

}

void robinconductivitydomainbc::set_coefs(const Real a_aco, const Real a_bco, const Real a_rhs){
  m_aco = a_aco;
  m_bco = a_bco;
  m_rhs = a_rhs;

  m_const_coeff = true;
  m_func_coeff  = false;
  m_data_coeff  = false;
}

void robinconductivitydomainbc::set_coefs(const RefCountedPtr<LevelData<EBFluxFAB> >& a_aco,
					  const RefCountedPtr<LevelData<EBFluxFAB> >& a_bco,
					  const RefCountedPtr<LevelData<EBFluxFAB> >& a_rhs){
  MayDay::Abort("robinconductivitydomainbc::set_coefs - data-based not supported (yet)");
  m_acodata = a_aco;
  m_bcodata = a_bco;
  m_rhsdata = a_rhs;

  m_const_coeff = false;
  m_func_coeff  = false;
  m_data_coeff  = true;
}

void robinconductivitydomainbc::getFaceFlux(BaseFab<Real>&        a_faceFlux, 
					    const BaseFab<Real>&  a_phi, 
					    const RealVect&       a_probLo, 
					    const RealVect&       a_dx, 
					    const int&            a_idir, 
					    const Side::LoHiSide& a_side, 
					    const DataIndex&      a_dit, 
					    const Real&           a_time, 
					    const bool&           a_useHomogeneous){
  CH_TIME("EddingtonSP1::EddingtonSP1DomainBC::getFaceFlux");
  CH_assert(a_phi.nComp() == 1);

  const int ncomp = 1;
  const int comp  = 0;
  const int iside = sign(flip(a_side)); // iside = 1 for Lo, iside = -1 for Hi. This is the normal vector sign.
  const Box& box  = a_faceFlux.box();

  for (BoxIterator bit(box); bit.ok(); ++bit){

    const IntVect iv     = bit();
    const IntVect ivNear = iv + 1*iside*BASISV(a_idir);
    const IntVect ivFar  = iv + 2*iside*BASISV(a_idir);

    const Real y0 = a_phi(iv,     comp); // Boundary node
    const Real y1 = a_phi(ivNear, comp); // Nearest neighbor
    const Real y2 = a_phi(ivFar,  comp); // Second nearest neighbor

    const Real x0 = 0.;                     
    const Real x1 = 1*iside*a_dx[a_idir];
    const Real x2 = 2*iside*a_dx[a_idir];

    const RealVect pos = a_probLo + a_dx*(RealVect(iv) - iside*0.5*RealVect(BASISV(a_idir)));
    Real aco, bco, rhs;
    if(m_const_coeff){
      aco = m_aco;
      bco = m_bco;
      rhs = m_rhs;
    }
    else if(m_func_coeff){
      aco = m_robinco->aco(pos);
      bco = m_robinco->aco(pos);
      rhs = m_robinco->rhs(pos);
    }
    else if(m_data_coeff){
      MayDay::Abort("robinconductivitydomainbc::getFaceFlux - data-based not supported (yet)");
    }

    if(a_useHomogeneous){ // Homogeneous bc
      rhs = 0.;
    }

    // Linear extrapolation
    const Real xb = -0.5*iside*a_dx[a_idir];
    const Real yb = y0 + (y1-y0)/(x1-x0)*(xb-x0);// + ((y2-y1)/((x2-x0)*(x2-x1)) - (y1-y0)/((x2-x0)*(x1-x0)))*(xb-x0)*(xb-x1);

    a_faceFlux(iv, comp) = rhs/bco + iside*aco*yb/bco; // dphi/dn = g/b - a*phi/b
  }


  // Since we are re-using the BC classes from EBAMRPoissonOp we must shift the flux box, the
  // input flux is cell-centered and the input is the box adjecent to the domain boundary on 
  // the valid side. We shift the flux box over and then multiply with the appropriate b-coefficient
  // and then shift back
  //
  // TLDR; we must shift the faceflux box when we multiply with the b-coefficient

  a_faceFlux.shiftHalf(a_idir, -sign(a_side));
  const Box& faceBox = a_faceFlux.box();
  const BaseFab<Real>& regCoef = (*m_bcoef)[a_dit][a_idir].getSingleValuedFAB();
  int  isrc = 0;
  int  idst = 0;
  int  inum = 1;
  FORT_MULTIPLYTWOFAB(CHF_FRA(a_faceFlux),
                      CHF_CONST_FRA(regCoef),
                      CHF_BOX(faceBox),
                      CHF_INT(isrc),CHF_INT(idst),CHF_INT(inum));
  a_faceFlux.shiftHalf(a_idir, sign(a_side));
}

void robinconductivitydomainbc::getFaceFlux(Real&                 a_faceFlux,
					    const VolIndex&       a_vof, 
					    const int&            a_comp, 
					    const EBCellFAB&      a_phi, 
					    const RealVect&       a_probLo, 
					    const RealVect&       a_dx, 
					    const int&            a_idir, 
					    const Side::LoHiSide& a_side, 
					    const DataIndex&      a_dit, 
					    const Real&           a_time, 
					    const bool&           a_useHomogeneous){
  CH_TIME("robinconductivitydomainbc::getFaceFlux");

  a_faceFlux = 0.0;

  const EBISBox& ebisbox          = a_phi.getEBISBox();
  const ProblemDomain& domain_box = ebisbox.getDomain();
  const Vector<FaceIndex> faces   = ebisbox.getFaces(a_vof,a_idir,a_side);
  
  if (faces.size() > 0) {
    if (faces.size() == 1) {
      IntVectSet cfivs;
      FaceStencil faceSten = EBArith::getInterpStencil(faces[0],
						       cfivs,
						       ebisbox,
						       domain_box);
      for (int isten = 0; isten < faceSten.size(); isten++) {
	const Real& weight      = faceSten.weight(isten);
	const FaceIndex& face   = faceSten.face(isten);
	const RealVect centroid = ebisbox.centroid(face);
	
	Real thisFaceFlux;

	this->getFaceGradPhi(thisFaceFlux,
			     face,
			     a_comp,
			     a_phi,
			     a_probLo,
			     a_dx,
			     a_idir,
			     a_side,
			     a_dit,
			     a_time,
			     false,
			     centroid,
			     a_useHomogeneous);

	a_faceFlux += thisFaceFlux*weight;
      }
      
      a_faceFlux *= ebisbox.areaFrac(faces[0]);
      
    }
    else {
      MayDay::Error("robinconductivitydomainbc::Multi-valued faces on domain edge");
    }
  }

  // Scale with area
  Real bcoave = 0;
  Real areaTot = 0;
  for (int iface = 0; iface < faces.size(); iface++) {
    Real areaFrac = ebisbox.areaFrac(faces[iface]);
    Real bcoFace  = (*m_bcoef)[a_dit][a_idir](faces[iface], 0);
    bcoave += areaFrac*bcoFace;
    areaTot += areaFrac;
  }
  if (areaTot > 1.0e-8) {
    bcoave /= areaTot;
  }
  a_faceFlux *= bcoave;
  
}

void robinconductivitydomainbc::getFaceGradPhi(Real&                 a_faceFlux, 
					       const FaceIndex&      a_face, 
					       const int&            a_comp, 
					       const EBCellFAB&      a_phi, 
					       const RealVect&       a_probLo, 
					       const RealVect&       a_dx, 
					       const int&            a_idir, 
					       const Side::LoHiSide& a_side, 
					       const DataIndex&      a_dit, 
					       const Real&           a_time, 
					       const bool&           a_useAreaFrac, 
					       const RealVect&       a_centroid, 
					       const bool&           a_useHomogeneous){
  const int ncomp        = 1;
  const int comp         = 0;
  const int iside        = sign(flip(a_side));
  const Box box          = a_phi.getRegion();
  const EBISBox& ebisbox = a_phi.getEBISBox();

  const VolIndex vof      = a_face.getVoF(flip(a_side));
  const IntVect iv        = vof.gridIndex();
  const IntVect iv_near   = iv + iside*BASISV(a_idir);
  const VolIndex vof_near = VolIndex(iv_near, comp); // Dangerous stuff if that is a multicell....
  CH_assert(!ebisbox.isMultiValued(iv_near));        // Warn about this. 

  // Linear extrapolation
  const Real x0 =  0.0;
  const Real x1 =  1.0*iside*a_dx[a_idir];
  const Real xb = -0.5*iside*a_dx[a_idir];
  const Real y0 =  a_phi(vof, comp);
  Real y1       =  y0;
  if(!ebisbox.isCovered(iv_near)){
    y1 = a_phi(vof_near, comp);
  }
  const Real yb = y0 + (y1 - y0)/(x1-x0)*(xb - x0);

  
  Real aco, bco, rhs;
  const RealVect pos = a_probLo + a_dx*(RealVect(iv) - iside*0.5*RealVect(BASISV(a_idir)));
  if(m_const_coeff){
    aco = m_aco;
    bco = m_bco;
    rhs = m_rhs;
  }
  else if(m_func_coeff){
    aco = m_robinco->aco(pos);
    bco = m_robinco->aco(pos);
    rhs = m_robinco->rhs(pos);
  }
  else if(m_data_coeff){
    MayDay::Abort("robinconductivitydomainbc::getFaceGradPhi - data-based not supported (yet)");
  }

  if(a_useHomogeneous){ // Homogeneous bc
    rhs = 0.;
  }

  a_faceFlux = rhs/bco - iside*aco*yb/bco;

  if(a_useAreaFrac){
    MayDay::Abort("robinconductivitydomainbc::getFaceGradPhi - useAreaFrac=true shouldn't happen");
  }
}

