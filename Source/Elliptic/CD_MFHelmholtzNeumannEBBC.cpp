/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzNeumannEBBC.cpp
  @brief  Implementation of CD_MFHelmholtzNeumannEBBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzNeumannEBBC.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzNeumannEBBC::MFHelmholtzNeumannEBBC(const int a_phase, const RefCountedPtr<JumpBC>& a_jumpBC) : MFHelmholtzEBBC(a_phase, a_jumpBC) {
  m_useConstant = false;
  m_useFunction = false;
}

MFHelmholtzNeumannEBBC::~MFHelmholtzNeumannEBBC(){

}

void MFHelmholtzNeumannEBBC::setDphiDn(const int a_DphiDn){
  m_multByBco   = true;
  m_useConstant = true;
  m_useFunction = false;
  
  m_constantDphiDn = a_DphiDn;
}

void MFHelmholtzNeumannEBBC::setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn){
  m_multByBco   = true;
  m_useConstant = false;
  m_useFunction = true;
  
  m_functionDphiDn = a_DphiDn;
}

void MFHelmholtzNeumannEBBC::setBxDphiDn(const int a_BxDphiDn){
  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}

void MFHelmholtzNeumannEBBC::setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn){
  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}
  
void MFHelmholtzNeumannEBBC::defineSinglePhase() {
  if(!(m_useConstant || m_useFunction)) MayDay::Error("MFHelmholtzNeumannEBBC - not using constant or function!");
}

void MFHelmholtzNeumannEBBC::applyEBFluxSinglePhase(VoFIterator&       a_singlePhaseVofs,
						    EBCellFAB&         a_Lphi,
						    const EBCellFAB&   a_phi,
						    const DataIndex&   a_dit,
						    const Real&        a_beta,
						    const bool&        a_homogeneousPhysBC) const {

  // TLDR: For Neumann, we want to add the flux beta*bco*area*(dphi/dn)/dx where the
  //       dx comes from the fact that the term we are computing will be added to kappa*div(F)
  for (a_singlePhaseVofs.reset(); a_singlePhaseVofs.ok(); ++a_singlePhaseVofs){
    const VolIndex& vof = a_singlePhaseVofs();
    
    if(!a_homogeneousPhysBC){
      Real value;
      if(m_useConstant){
	value = m_constantDphiDn;
      }
      else if(m_useFunction){
	const RealVect pos = this->getBoundaryPosition(vof, a_dit);
	value = m_functionDphiDn(pos);
      }

      // B-coefficient, area fraction, and division by dx (from Div(F)) already a part of the boundary weights, but
      // beta is not.
      const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
      const Real areaFrac    = ebisbox.bndryArea(vof);
      const Real B           = m_multByBco ? (*m_Bcoef)[a_dit](vof, m_comp) : 1;
      const Real kappaDivF   = a_beta*B*value*areaFrac/m_dx;
  
      a_Lphi(vof, m_comp) += kappaDivF;
    }
  }
  
  return;
}

#include <CD_NamespaceFooter.H>
