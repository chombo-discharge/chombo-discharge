/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletDomainBC.cpp
  @brief  Implementation of CD_EBHelmholtzDirichletDomainBC.cpp
  @author Robert Marskar
*/

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

void EBHelmholtzDirichletDomainBC::getFaceFlux(FArrayBox&            a_faceFlux,
					       const FArrayBox&      a_phi,
					       const int&            a_dir,
					       const Side::LoHiSide& a_side,
					       const DataIndex&      a_dit,
					       const bool            a_useHomogeneous) const {
  const Box cellbox        = a_faceFlux.box();
  const BaseFab<Real>& Bco = (*m_Bcoef)[a_dit][a_dir].getSingleValuedFAB();
  const int isign          = sign(a_side); 

  if(a_useHomogeneous){
    const Real value = 0.0;

    FORT_HELMHOLTZDIRICHLETFLUX(CHF_FRA1(a_faceFlux, m_comp),
				CHF_CONST_FRA1(a_phi, m_comp),
				CHF_CONST_FRA1(Bco, m_comp),
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
				  CHF_CONST_FRA1(Bco, m_comp),
				  CHF_CONST_REAL(value),
				  CHF_CONST_REAL(m_dx),				    
				  CHF_CONST_INT(a_dir),
				  CHF_CONST_INT(isign),
				  CHF_BOX(cellbox));
    }
    else if(m_useFunction){
      MayDay::Abort("EBHelmholtzDirichletDomainBC::getFaceFlux -- function not supported (yet)");
    }
    else{
      MayDay::Abort("EBHelmholtzDirichletDomainBC::getFaceFlux -- logic bust");
    }
  }
}

void EBHelmholtzDirichletDomainBC::getFaceFlux(Real&                 a_faceFlux,
					       const VolIndex&       a_vof,
					       const EBCellFAB&      a_phi,
					       const int&            a_dir,
					       const Side::LoHiSide& a_side,
					       const DataIndex&      a_dit,
					       const bool            a_useHomogeneous) const {

  MayDay::Error("EBHelmholtzDirichletDomainBC::getFaceFlux - not implemented");
}

#include <CD_NamespaceFooter.H>
