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
    const Box box          = dbl[dit()];
    const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs  = ebisbox.getIrregIVS(box);

    BaseIVFAB<Real>&       weights  = m_boundaryWeights  [dit()];
    BaseIVFAB<VoFStencil>& stencils = m_kappaDivFStencils[dit()];

    const BaseIVFAB<Real>& Bcoef    = (*m_Bcoef)[dit()];

    weights. define(ivs, ebgraph, m_nComp);
    stencils.define(ivs, ebgraph, m_nComp);

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const Real areaFrac = ebisbox.bndryArea(vof);
      const Real B        = Bcoef(vof, m_comp);

      bool foundStencil = false;
      std::pair<Real, VoFStencil> pairSten;

      if(m_order == 1){
	foundStencil = this->getFirstOrderStencil(pairSten, vof, dit());
      }
      else if(m_order == 2){
	foundStencil = this->getSecondOrderStencil(pairSten, vof, dit());
	if(!foundStencil){
	  foundStencil = this->getFirstOrderStencil(pairSten, vof, dit());
	}
      }

      if(foundStencil){
	weights (vof, m_comp) = pairSten.first;
	stencils(vof, m_comp) = pairSten.second;
	
	// Stencil and weight must also be scaled by the b-coefficient, dx and the area fraction. 
	weights (vof, m_comp) *= B*areaFrac/m_dx;
	stencils(vof, m_comp) *= B*areaFrac/m_dx;
      }
      else{
	MayDay::Warning("EBHelmholtzDirichletEBBC::define - did not find stencil and this vof will default to homogeneous Neumann!");
	
	weights(vof, m_comp) = 0.0;
	stencils(vof, m_comp).clear();
      }
    }
  }
}

bool EBHelmholtzDirichletEBBC::getFirstOrderStencil(std::pair<Real, VoFStencil>& a_stencil, const VolIndex& a_vof, const DataIndex& a_dit){
  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
  const RealVect normal  = ebisbox.normal(a_vof);
  const RealVect point   = ebisbox.bndryCentroid(a_vof);

  Real w;
  VoFStencil sten;

  EBArith::getLeastSquaresGradSten(sten, w, a_vof, ebisbox, m_dx*RealVect::Unit, m_eblg.getDomain(), m_comp);

  a_stencil = std::make_pair(w, sten);

  return true;
}

bool EBHelmholtzDirichletEBBC::getSecondOrderStencil(std::pair<Real, VoFStencil>& a_stencil, const VolIndex& a_vof, const DataIndex& a_dit){
  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];

  bool foundStencil = false;

  bool dropOrder;
  Vector<VoFStencil> pointStencils;
  Vector<Real>       distanceAlongLine;
  EBArith::johanStencil(dropOrder, pointStencils, distanceAlongLine, a_vof, ebisbox, m_dx*RealVect::Unit, (*m_eblg.getCFIVS())[a_dit]);

  if(!dropOrder){ // We know we can at least a second order approximation

    const Real x1    = distanceAlongLine[0];
    const Real x2    = distanceAlongLine[1];
    const Real denom = x2*x2*x1 - x1*x1*x2;

    const Real weight = -x1*x1/denom + x2*x2/denom;

    VoFStencil sten;
    VoFStencil sten1 = pointStencils[0];
    VoFStencil sten2 = pointStencils[1];
    
    sten1 *= -x2*x2/denom;
    sten2 *= +x1*x1/denom;
    
    sten  += sten1;
    sten  += sten2;

    a_stencil = std::make_pair(weight, sten);

    foundStencil = true;
  }

  return foundStencil;
}

void EBHelmholtzDirichletEBBC::applyEBFlux(EBCellFAB&         a_Lphi,
					   const EBCellFAB&   a_phi,
					   const VolIndex&    a_vof,
					   const DataIndex&   a_dit,
					   const Real&        a_beta) const {

  Real value;
  if(m_useConstant){
    value = m_constantValue;
  }
  else if(m_useFunction){
    const RealVect pos = this->getBoundaryPosition(a_vof, a_phi, m_dx);
    value = m_functionValue(pos);
  }

  // B-coefficient, area fraction, and division by dx (from Div(F)) already a part of the boundary weights, but
  // beta is not. 
  a_Lphi(a_vof, m_comp) += a_beta*value*m_boundaryWeights[a_dit](a_vof, m_comp);
  
  return;
}

#include <CD_NamespaceFooter.H>
