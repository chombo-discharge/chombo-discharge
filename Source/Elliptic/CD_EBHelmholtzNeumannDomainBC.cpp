/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzNeumannDomainBC.cpp
  @brief  Implementation of CD_EBHelmholtzNeumannDomainBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzOpF_F.H>
#include <CD_EBHelmholtzNeumannDomainBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzNeumannDomainBC::EBHelmholtzNeumannDomainBC(){
  m_multByBco   = true;
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzNeumannDomainBC::EBHelmholtzNeumannDomainBC(const Real a_DphiDn){
  this->setDphiDn(a_DphiDn);
}

EBHelmholtzNeumannDomainBC::EBHelmholtzNeumannDomainBC(const std::function<Real(const RealVect& a_pos)>& a_DphiDn){
  this->setDphiDn(a_DphiDn);
}

EBHelmholtzNeumannDomainBC::~EBHelmholtzNeumannDomainBC(){

}

void EBHelmholtzNeumannDomainBC::setDphiDn(const int a_DphiDn){
  m_multByBco   = false;
  m_useConstant = true;
  m_useFunction = false;
  
  m_constantDphiDn = a_DphiDn;
}


void EBHelmholtzNeumannDomainBC::setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn){
  m_multByBco   = false;
  m_useConstant = false;
  m_useFunction = true;
  
  m_functionDphiDn = a_DphiDn;
}

void EBHelmholtzNeumannDomainBC::setBxDphiDn(const int a_BxDphiDn){
  this->setDphiDn(a_BxDphiDn);
  m_multByBco = false;
}

void EBHelmholtzNeumannDomainBC::setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn){
  this->setDphiDn(a_BxDphiDn);
  m_multByBco = false;
}

void EBHelmholtzNeumannDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
					     const BaseFab<Real>&  a_phi,
					     const int&            a_dir,
					     const Side::LoHiSide& a_side,
					     const DataIndex&      a_dit,
					     const bool            a_useHomogeneous) const {

  const Box cellbox        = a_faceFlux.box();
  const BaseFab<Real>& Bco = (*m_Bcoef)[a_dit][a_dir].getSingleValuedFAB();
  const int isign          = (a_side == Side::Lo) ? -1 : 1;
  
  if(a_useHomogeneous){
    a_faceFlux.setVal(0.0);
  }
  else{
    // Note the conspicuous minus-sign. It's because we assume that a positive DphiDn puts a flux INTO the domain. 
    if(m_useConstant){
      a_faceFlux.setVal(-isign*m_constantDphiDn); 
    }
    else if(m_useFunction){
      for (BoxIterator bit(cellbox); bit.ok(); ++bit){
	const IntVect iv   = bit();
	const RealVect pos = this->getBoundaryPosition(iv, a_dir, a_side);

	const Real DphiDn  = m_functionDphiDn(pos);

	a_faceFlux(iv, m_comp) = -isign * DphiDn;
      }
    }

    // Multiply by B-coefficient. 
    if(m_multByBco){
      FORT_HELMHOLTZMULTFLUXBYBCO(CHF_FRA1(a_faceFlux, m_comp),
				  CHF_CONST_FRA1(Bco, m_comp),
				  CHF_CONST_INT(a_dir),
				  CHF_CONST_INT(isign),
				  CHF_BOX(cellbox));      
    }
  }

}


Real EBHelmholtzNeumannDomainBC::getFaceFlux(const VolIndex&       a_vof,
					     const EBCellFAB&      a_phi,
					     const int&            a_dir,
					     const Side::LoHiSide& a_side,
					     const DataIndex&      a_dit,
					     const bool            a_useHomogeneous) const {

  return 0.0;
}

#include <CD_NamespaceFooter.H>
