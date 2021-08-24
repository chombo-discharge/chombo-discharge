/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzRobinDomainBC.cpp
  @brief  Implementation of CD_EBHelmholtzRobinDomainBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>
#include <CH_Timer.H>

// Our includes
#include <CD_EBHelmholtzOpF_F.H>
#include <CD_EBHelmholtzRobinDomainBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzRobinDomainBC::EBHelmholtzRobinDomainBC(){
  CH_TIME("EBHelmholtzRobinDomainBC::EBHelmholtzRobinDomainBC()");
  
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzRobinDomainBC::~EBHelmholtzRobinDomainBC(){
  CH_TIME("EBHelmholtzRobinDomainBC::~EBHelmholtzRobinDomainBC()");
}

void EBHelmholtzRobinDomainBC::setCoefficients(const Real a_A, const Real a_B, const Real a_C){
  CH_TIME("EBHelmholtzRobinDomainBC::setCoefficients(Real, Real, Real)");
  
  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
}


void EBHelmholtzRobinDomainBC::setCoefficients(const std::function<Real(const RealVect& a_pos) >& a_A,
					       const std::function<Real(const RealVect& a_pos) >& a_B,
					       const std::function<Real(const RealVect& a_pos) >& a_C){
  CH_TIME("EBHelmholtzRobinDomainBC::setCoefficients(3x std::function<Real(RealVect)>)");
  
  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;
}

void EBHelmholtzRobinDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
					   const BaseFab<Real>&  a_phi,
					   const int&            a_dir,
					   const Side::LoHiSide& a_side,
					   const DataIndex&      a_dit,
					   const bool            a_useHomogeneous) const {
  CH_TIME("EBHelmholtzRobinDomainBC::getFaceFlux(BaseFab<Real>, BaseFab<Real>, int, Side::LoHiSide, DataIndex, bool)");
  
  const Box cellbox = a_faceFlux.box();
  const int isign   = (a_side == Side::Lo) ? -1 : 1;
  
  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
  const EBISBox& ebgraph = m_eblg.getEBISL()[a_dit];

  for (BoxIterator bit(cellbox); bit.ok(); ++bit){
    const IntVect iv     = bit();
    const IntVect ivNear = iv - isign*BASISV(a_dir);

    const bool hasNear = !ebisbox.isCovered(ivNear);

    // Linear extrapolation if we can (should always be possible AFAIK, since we're dealing with regular cells). 
    Real phiExtrap;
    if(hasNear) { 
      const Real y0 = a_phi(iv,     m_comp);
      const Real y1 = a_phi(ivNear, m_comp);
      
      phiExtrap = 1.5*y0 - 0.5*y1;
    }
    else{
      phiExtrap = a_phi(iv, m_comp);
    }

    Real A;
    Real B;
    Real C;
    if(m_useConstant){
      A = m_constantA;
      B = m_constantB;
      C = a_useHomogeneous ? 0 : m_constantC;
    }
    else if(m_useFunction){
      const RealVect pos = this->getBoundaryPosition(iv, a_dir, a_side);

      A = this->m_functionA(pos);
      B = this->m_functionB(pos);
      C = a_useHomogeneous ? 0 : this->m_functionC(pos);
    }

    //    a_faceFlux(iv, m_comp) = C/B + isign*A*phiExtrap/B;
    a_faceFlux(iv, m_comp) = C/B + isign*A*phiExtrap/B;
  }

  // Multiplies by B-coefficient.   
  const BaseFab<Real>& Bco = (*m_Bcoef)[a_dit][a_dir].getSingleValuedFAB();
  FORT_HELMHOLTZMULTFLUXBYBCO(CHF_FRA1(a_faceFlux, m_comp),
			      CHF_CONST_FRA1(Bco, m_comp),
			      CHF_CONST_INT(a_dir),
			      CHF_CONST_INT(isign),
			      CHF_BOX(cellbox));      
}
  
Real EBHelmholtzRobinDomainBC::getFaceFlux(const VolIndex&       a_vof,
					   const EBCellFAB&      a_phi,
					   const int&            a_dir,
					   const Side::LoHiSide& a_side,
					   const DataIndex&      a_dit,
					   const bool            a_useHomogeneous) const {
  CH_TIME("EBHelmholtzRobinDomainBC::getFaceFlux(VolIndex, EBCellFAB, int, Side::LoHiSide, DataIndex, bool)");
  
  const int isign             = (a_side == Side::Lo) ? -1 : 1;
  const Real ihdx             = 2.0/m_dx;
  const IntVect iv            = a_vof.gridIndex();
  const EBISBox& ebisbox      = m_eblg.getEBISL()[a_dit];
  const ProblemDomain& domain = m_eblg.getDomain();

  Real centroidFlux = 0.0;
  
  const Vector<FaceIndex> faces = ebisbox.getFaces(a_vof, a_dir, a_side);

  if(faces.size() > 0){

    if(faces.size() == 1){ // Get an interpolation stencil, using centered differences
      const FaceStencil faceSten = EBArith::getInterpStencil(faces[0], IntVectSet(), ebisbox, domain);

      Real phiExtrap = 0.0;
      
      // This loop will extrapolate phi to the boundary faces in the interpolation stencil. If we don't have
      // enough cells we will just use the value a_phi(a_vof, m_comp), but linearly extrapolate otherwise. If
      // we encounter a multi-valued cell we take the average of the values in that cell. 
      for (int i = 0; i < faceSten.size(); i++){
	const Real& weight     = faceSten.weight(i);
	const FaceIndex& face  = faceSten.face(i);
	const VolIndex& curVof = face.getVoF(flip(a_side));

	const Vector<VolIndex> nearVofs = ebisbox.getVoFs(curVof, a_dir, flip(a_side), 1);

	Real curExtrap;
	if(nearVofs.size() > 0){

	  // In case we have a multi-cell, do the average. 
	  Real nearPhi = 0.0;
	  for (int i = 0; i < nearVofs.size(); i++){
	    nearPhi += a_phi(nearVofs[i], m_comp);
	  }
	  nearPhi = nearPhi/nearVofs.size();

	  // Linear extrapolation to boundary from curVof
	  curExtrap = 1.5*a_phi(curVof, m_comp) - 0.5*nearPhi;
	}
	else{
	  curExtrap = a_phi(curVof, m_comp);
	}

	// Apply interpolation weight. 
	phiExtrap += curExtrap * weight;
      }

      // Multiply by the appropriate Robin coefficients
      Real A;
      Real B;
      Real C;
      if(m_useConstant){
	A = m_constantA;
	B = m_constantB;
	C = a_useHomogeneous ? 0 : m_constantC;
      }
      else if(m_useFunction){
	const RealVect pos = this->getBoundaryPosition(iv, a_dir, a_side);

	A = this->m_functionA(pos);
	B = this->m_functionB(pos);
	C = a_useHomogeneous ? 0 : this->m_functionC(pos);
      }

      // A*phi + B*dphi/dn = C => dphi/dn = C/B - A*phi/B, with n pointing into the boundary. 
      centroidFlux = C/B + isign*A*phiExtrap/B;

      // Multiply by the Helmholtz b-coefficient and aperture.
      const FaceIndex& bndryFace = faces[0];
      const Real Bco  = (*m_Bcoef)[a_dit][a_dir](bndryFace, m_comp);
      const Real area = ebisbox.areaFrac(bndryFace);

      centroidFlux *= Bco*area;
    }
    else{
      MayDay::Error("EBHelmholtzDirichletDomainBC -- boundary face is multivalued and EBHelmholtzDirichletDomainBC does not support this (yet)");
    }
  }

  return centroidFlux;

}

#include <CD_NamespaceFooter.H>
