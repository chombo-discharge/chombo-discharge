/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzRobinEBBC.cpp
  @brief  Implementation of CD_EBHelmholtzRobinEBBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_EBHelmholtzRobinEBBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzRobinEBBC::EBHelmholtzRobinEBBC(){
  m_useConstant = false;
  m_useFunction = false;
}


EBHelmholtzRobinEBBC::~EBHelmholtzRobinEBBC(){

}
  
void EBHelmholtzRobinEBBC::setCoefficients(const Real a_A, const Real a_B, const Real a_C){
  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
}

void EBHelmholtzRobinEBBC::setCoefficients(const std::function<Real(const RealVect& a_pos) >& a_A,
					   const std::function<Real(const RealVect& a_pos) >& a_B,
					   const std::function<Real(const RealVect& a_pos) >& a_C){
  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;
}

VoFStencil EBHelmholtzRobinEBBC::getExtrapolationStencil(const VolIndex& a_vof, const DataIndex& a_dit) const {
  VoFStencil extrapStencil;
  
  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
  const RealVect dist    = ebisbox.bndryCentroid(a_vof)*m_dx;

  const int order = EBArith::getFirstOrderExtrapolationStencil(extrapStencil, dist, m_dx*RealVect::Unit, a_vof, ebisbox, -1, &(*m_eblg.getCFIVS())[a_dit], m_comp);

  if(order == 0) MayDay::Error("EBHelmholtzRobinEBBC::getExtrapolationStencil - could not find stencil!");

  extrapStencil.clear();
  extrapStencil.add(a_vof, 1.0);
  
  return extrapStencil;
}

void EBHelmholtzRobinEBBC::define() {
  if(!(m_useConstant || m_useFunction)) MayDay::Error("EBHelmholtzRobinEBBC::define - not using constant or function!");

  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  m_kappaDivFStencils.define(dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box box          = dbl[dit()];
    const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs  = ebisbox.getIrregIVS(box);

    BaseIVFAB<VoFStencil>& stencils = m_kappaDivFStencils[dit()];

    stencils.define(ivs, ebgraph, m_nComp);

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const Real areaFrac = ebisbox.bndryArea(vof);
      const Real helmBco  = (*m_Bcoef)[dit()](vof, m_comp);

      stencils(vof, m_comp) = this->getExtrapolationStencil(vof, dit());

      Real A;
      Real B;
      if(m_useConstant){
	A = m_constantA;
	B = m_constantB;
      }
      else if(m_useFunction){
	const RealVect pos = this->getBoundaryPosition(vof, dit());
	A = m_functionA(pos);
	B = m_functionB(pos);
      }

      // The normal derivative is dphi/dn = (A*phi - C)/B and the (stencil) flux is
      // kappaDivF = area*b*dphidn/Delta x. Scale accordingly.
      if(std::abs(B) > 0.0){
	stencils(vof, m_comp) *= A*areaFrac*helmBco/(B*m_dx);
      }
      else{
	stencils(vof, m_comp).clear();
      }
    }
  }
}

void EBHelmholtzRobinEBBC::applyEBFlux(EBCellFAB&         a_Lphi,
				       const EBCellFAB&   a_phi,
				       const VolIndex&    a_vof,
				       const DataIndex&   a_dit,
				       const Real&        a_beta) const {

  // Recall that the "flux" is kappaDivF = area*dphi/dn/DeltaX where dphi/dn = (A*phi - C)/B. We already have the phi
  // term in the stencil so only need to add -C/B.

  Real B;
  Real C;
  if(m_useConstant){
    B = m_constantB;
    C = m_constantC;
  }
  else if(m_useFunction){
    const RealVect pos = this->getBoundaryPosition(a_vof, a_dit);
    B = m_functionB(pos);
    C = m_functionC(pos);
  }

  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
  const Real areaFrac    = ebisbox.bndryArea(a_vof);
  const Real helmBco     = (*m_Bcoef)[a_dit](a_vof, m_comp);
  const Real kappaDivF   = -a_beta*helmBco*areaFrac*C/(m_dx*B);

  if(std::abs(B) > 0.0){
    a_Lphi(a_vof, m_comp) += kappaDivF;
  }
}

#include <CD_NamespaceFooter.H>
