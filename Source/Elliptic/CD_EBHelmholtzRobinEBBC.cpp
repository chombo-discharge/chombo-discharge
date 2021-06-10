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
  m_useData     = false;
}


EBHelmholtzRobinEBBC::~EBHelmholtzRobinEBBC(){

}
  
void EBHelmholtzRobinEBBC::setCoefficients(const Real a_A, const Real a_B, const Real a_C){
  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
  m_useData     = false;
}

void EBHelmholtzRobinEBBC::setCoefficients(const std::function<Real(const RealVect& a_pos) >& a_A,
					   const std::function<Real(const RealVect& a_pos) >& a_B,
					   const std::function<Real(const RealVect& a_pos) >& a_C){
  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;
  m_useData     = false;
}

void EBHelmholtzRobinEBBC::setCoefficients(const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_A,
					   const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_B,
					   const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_C){
  m_dataA = a_A;
  m_dataB = a_B;
  m_dataC = a_C;

  m_useConstant = false;
  m_useFunction = false;
  m_useData     = true;
}

VoFStencil EBHelmholtzRobinEBBC::getExtrapolationStencil(const VolIndex& a_vof, const DataIndex& a_dit) const {
  VoFStencil extrapStencil;
  
  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
  const RealVect dist    = ebisbox.bndryCentroid(a_vof)*m_dx;

  IntVectSet ivs = IntVectSet();

  EBArith::getExtrapolationStencil(extrapStencil, dist, m_dx*RealVect::Unit, a_vof, ebisbox, -1, &ivs, m_comp);

  return extrapStencil;
}

void EBHelmholtzRobinEBBC::define() {
  if(!m_useConstant || !m_useFunction || !m_useData) MayDay::Error("EBHelmholtzRobinEBBC::define - not using constant, function, or data!");

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
      else if(m_useData){
	A = (*m_dataA)[dit()](vof, m_comp);
	B = (*m_dataA)[dit()](vof, m_comp);
      }

      // The normal derivative is dphi/dn = (A*phi - C)/B and the (stencil) flux is
      // kappaDivF = area*b*dphidn/Delta x. Scale accordingly. 
      stencils(vof, m_comp) *= A*areaFrac*helmBco/(B*m_dx);
    }
  }
}

void EBHelmholtzRobinEBBC::applyEBFlux(EBCellFAB&         a_Lphi,
				       const EBCellFAB&   a_phi,
				       const VolIndex&    a_vof,
				       const DataIndex&   a_dit,
				       const Real&        a_beta) const {

}

#include <CD_NamespaceFooter.H>
