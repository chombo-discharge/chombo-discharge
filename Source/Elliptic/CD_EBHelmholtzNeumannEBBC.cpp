/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzNeumannEBBC.cpp
  @brief  Implementation of CD_EBHelmholtzNeumannEBBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzNeumannEBBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzNeumannEBBC::EBHelmholtzNeumannEBBC(){
  m_multByBco   = true;
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzNeumannEBBC::~EBHelmholtzNeumannEBBC(){

}

void EBHelmholtzNeumannEBBC::setDphiDn(const int a_DphiDn){
  m_multByBco   = true;
  m_useConstant = true;
  m_useFunction = false;
  
  m_constantDphiDn = a_DphiDn;
}

void EBHelmholtzNeumannEBBC::setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn){
  m_multByBco   = true;
  m_useConstant = false;
  m_useFunction = true;
  
  m_functionDphiDn = a_DphiDn;
}

void EBHelmholtzNeumannEBBC::setBxDphiDn(const int a_BxDphiDn){
  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}

void EBHelmholtzNeumannEBBC::setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn){
  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}
  
void EBHelmholtzNeumannEBBC::define() {
  if(!(m_useConstant || m_useFunction)) MayDay::Error("EBHelmholtzNeumannEBBC::define - logic bist");
  
  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  // Reset the stencil everywhere. 
  m_kappaDivFStencils.define(dbl);
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box box          = dbl[dit()];
    const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs  = ebisbox.getIrregIVS(box);

    BaseIVFAB<VoFStencil>& stencils = m_kappaDivFStencils[dit()];

    stencils.define(ivs, ebgraph, m_nComp);

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      stencils(vofit(), m_comp).clear();
    }
  }
}

void EBHelmholtzNeumannEBBC::applyEBFlux(VoFIterator&       a_vofit,
					 EBCellFAB&         a_Lphi,
					 const EBCellFAB&   a_phi,
					 const DataIndex&   a_dit,
					 const Real&        a_beta,
					 const bool&        a_homogeneousPhysBC) const {

  // TLDR: For Neumann, we want to add the flux beta*bco*area*(dphi/dn)/dx where the
  //       dx comes from the fact that the term we are computing will be added to kappa*div(F)
  if(!a_homogeneousPhysBC){
    for (a_vofit.reset(); a_vofit.ok(); ++a_vofit){
      const VolIndex& vof = a_vofit();
    
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
