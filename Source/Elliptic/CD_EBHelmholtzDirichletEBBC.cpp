/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletEBBC.cpp
  @brief  Implementation of CD_EBHelmholtzDiricheltEBBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_LeastSquares.H>
#include <CD_EBHelmholtzDirichletEBBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzDirichletEBBC::EBHelmholtzDirichletEBBC(){
  m_order       = -1;
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzDirichletEBBC::~EBHelmholtzDirichletEBBC(){

}

void EBHelmholtzDirichletEBBC::setOrder(const int a_order){
  m_order = a_order;
}

void EBHelmholtzDirichletEBBC::setValue(const int a_value){
  m_useConstant = true;
  m_useFunction = false;
  
  m_constantValue = a_value;
}

void EBHelmholtzDirichletEBBC::setValue(const std::function<Real(const RealVect& a_pos)>& a_value){
  m_useConstant = false;
  m_useFunction = true;
  
  m_functionValue = a_value;
}
  
void EBHelmholtzDirichletEBBC::define() {
  if(!(m_order == 1  || m_order == 2 )) MayDay::Error("EBHelmholtzDirichletEBBC - order is not 1 or 2!");
  if(!(m_useConstant || m_useFunction)) MayDay::Error("EBHelmholtzDirichletEBBC - not using constant or function!");

  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  m_boundaryWeights.define(dbl);
  m_kappaDivFStencils.define(dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box box = dbl[dit()];
    const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs  = ebisbox.getIrregIVS(box);

    BaseIVFAB<Real>&       weights  = m_boundaryWeights  [dit()];
    BaseIVFAB<VoFStencil>& stencils = m_kappaDivFStencils[dit()];

    weights. define(ivs, ebgraph, m_nComp);
    stencils.define(ivs, ebgraph, m_nComp);

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const Real areaFrac = ebisbox.bndryArea(vof);
      
      // Stencil and weight must be scaled by dx and the area fraction. 
      weights (vof, m_comp) *= areaFrac/m_dx;
      stencils(vof, m_comp) *= areaFrac/m_dx;
    }
  }
}

void EBHelmholtzDirichletEBBC::applyEBFlux(EBCellFAB&         a_Lphi,
					   const EBCellFAB&   a_phi,
					   const VolIndex&    a_vof,
					   const Real&        a_beta) const {

  Real flux = 0.0;

  // B-coefficient already a part of the stencil and weight, but beta is not. 
  flux *= a_beta;

  // 
  a_Lphi(a_vof, m_comp) += flux;
  
  return;
}

#include <CD_NamespaceFooter.H>
