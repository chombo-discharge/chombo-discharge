/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletDomainBC.cpp
  @brief  Implementation of CD_EBHelmholtzDirichletDomainBC.cpp
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_EBHelmholtzOpF_F.H>
#include <CD_EBHelmholtzDirichletDomainBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzDirichletDomainBC::EBHelmholtzDirichletDomainBC(){
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzDirichletDomainBC::EBHelmholtzDirichletDomainBC(const Real a_value){
  this->setValue(a_value);
}

EBHelmholtzDirichletDomainBC::EBHelmholtzDirichletDomainBC(const std::function<Real(const RealVect& a_pos)>& a_value){
  this->setValue(a_value);
}

EBHelmholtzDirichletDomainBC::~EBHelmholtzDirichletDomainBC(){

}

void EBHelmholtzDirichletDomainBC::setValue(const Real a_value){
  m_useConstant = true;
  m_useFunction = false;
  
  m_constantValue = a_value;
}

void EBHelmholtzDirichletDomainBC::setValue(const std::function<Real(const RealVect& a_pos)>& a_value){
  m_useConstant = false;
  m_useFunction = true;
  
  m_functionValue = a_value;
}

void EBHelmholtzDirichletDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
					       const BaseFab<Real>&  a_phi,
					       const int&            a_dir,
					       const Side::LoHiSide& a_side,
					       const DataIndex&      a_dit,
					       const bool            a_useHomogeneous) const {
  const Box cellbox        = a_faceFlux.box();
  const BaseFab<Real>& Bco = (*m_Bcoef)[a_dit][a_dir].getSingleValuedFAB();
  const int isign          = (a_side == Side::Lo) ? -1 : 1;

  if(a_useHomogeneous){
    const Real value = 0.0;

    FORT_HELMHOLTZDIRICHLETFLUX(CHF_FRA1(a_faceFlux, m_comp),
				CHF_CONST_FRA1(a_phi, m_comp),
				CHF_CONST_REAL(value),
				CHF_CONST_REAL(m_dx),				  
				CHF_CONST_INT(a_dir),
				CHF_CONST_INT(isign),
				CHF_BOX(cellbox));
  }
  else{
    if(m_useConstant){
      const Real value = m_constantValue;

      FORT_HELMHOLTZDIRICHLETFLUX(CHF_FRA1(a_faceFlux, m_comp),
				  CHF_CONST_FRA1(a_phi, m_comp),
				  CHF_CONST_REAL(value),
				  CHF_CONST_REAL(m_dx),				    
				  CHF_CONST_INT(a_dir),
				  CHF_CONST_INT(isign),
				  CHF_BOX(cellbox));
    }
    else if(m_useFunction){
      const Real ihdx = 2.0/m_dx;
      for (BoxIterator bit(cellbox); bit.ok(); ++bit){
	const IntVect iv = bit();

	const RealVect pos = this->getBoundaryPosition(iv, a_dir, a_side);
	const Real value   = m_functionValue(pos);

	a_faceFlux(iv, m_comp) = isign * ihdx * (value - a_phi(iv, m_comp));
      }
    }
    else{
      MayDay::Abort("EBHelmholtzDirichletDomainBC::getFaceFlux -- logic bust");
    }
  }

  FORT_HELMHOLTZMULTFLUXBYBCO(CHF_FRA1(a_faceFlux, m_comp),
			      CHF_CONST_FRA1(Bco, m_comp),
			      CHF_CONST_INT(a_dir),
			      CHF_CONST_INT(isign),
			      CHF_BOX(cellbox));      
}

Real EBHelmholtzDirichletDomainBC::getFaceFlux(const VolIndex&       a_vof,
					       const EBCellFAB&      a_phi,
					       const int&            a_dir,
					       const Side::LoHiSide& a_side,
					       const DataIndex&      a_dit,
					       const bool            a_useHomogeneous) const {

  const int isign             = (a_side == Side::Lo) ? -1 : 1;
  const Real ihdx             = 2.0/m_dx;
  const IntVect iv            = a_vof.gridIndex();
  const EBISBox& ebisbox      = m_eblg.getEBISL()[a_dit];
  const ProblemDomain& domain = m_eblg.getDomain();

  Real centroidFlux = 0.0;
  
  const Vector<FaceIndex> faces = ebisbox.getFaces(a_vof, a_dir, a_side);

  if(faces.size() > 0){
    if(faces.size() == 1){ // Get an interpolation stencil, using centered differences
      IntVectSet cfivs;
      const FaceStencil faceSten = EBArith::getInterpStencil(faces[0], cfivs, ebisbox, domain);

      for (int i = 0; i < faceSten.size(); i++){
	const Real& weight    = faceSten.weight(i);
	const FaceIndex& face = faceSten.face(i);

	// Get dphi/dx on the boundary
	Real value;
	if(a_useHomogeneous){
	  value = 0.0;
	}
	else {
	  if(m_useConstant){
	    value = m_constantValue;
	  }
	  else if(m_useFunction){
	    value = m_functionValue(this->getBoundaryPosition(iv, a_dir, a_side));
	  }
	}

	const Real centeredFaceFlux = isign * ihdx * (value - a_phi(a_vof, m_comp));

	centroidFlux += centeredFaceFlux * weight;
      }

      // Multiply by b-coefficient and aperture.
      const FaceIndex& bndryFace = faces[0];
      const Real Bco  = (*m_Bcoef)[a_dit][a_dir](bndryFace, m_comp);
      const Real area = ebisbox.areaFrac(bndryFace);

      centroidFlux *= Bco*area;
    }
    else{
      MayDay::Error("EBHelmholtzDirichletDomainBC -- boundary face is multivalued and EBHelmholtzDirichletDomainBC does not supportd that (yet)");
    }
  }

  return centroidFlux;
}

#include <CD_NamespaceFooter.H>
