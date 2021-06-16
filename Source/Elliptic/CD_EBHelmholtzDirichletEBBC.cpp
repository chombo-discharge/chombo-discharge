/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletEBBC.cpp
  @brief  Implementation of CD_EBHelmholtzDirichletEBBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>
#include <NeighborIterator.H>

// Our includes
#include <CD_LeastSquares.H>
#include <CD_EBHelmholtzDirichletEBBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzDirichletEBBC::EBHelmholtzDirichletEBBC(){
  m_order       = -1;
  m_weight      = -1;
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzDirichletEBBC::~EBHelmholtzDirichletEBBC(){

}

void EBHelmholtzDirichletEBBC::setOrder(const int a_order){
  CH_assert(a_order > 0);
  m_order = a_order;
}

void EBHelmholtzDirichletEBBC::setWeight(const int a_weight){
  CH_assert(a_weight > 0);
  m_weight = a_weight;
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
  if(  m_order <= 0  || m_weight <= 0 ) MayDay::Error("EBHelmholtzDirichletEBBC - must have order > 0 and weight > 0");
  if(!(m_useConstant || m_useFunction)) MayDay::Error("EBHelmholtzDirichletEBBC - not using constant or function!");

  const DisjointBoxLayout& dbl = m_eblg.getDBL();
  const ProblemDomain& domain  = m_eblg.getDomain();
  
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

      int order;
      bool foundStencil = false;
      std::pair<Real, VoFStencil> pairSten;

      // Try quadrants first. 
      order = m_order;
      while(!foundStencil && order > 0){
      	foundStencil = this->getLeastSquaresStencil(pairSten, vof, VofUtils::Neighborhood::Quadrant, dit(), order);
      	order--;

	// Check if stencil reaches too far across CF
	if(foundStencil) {
	  foundStencil = this->isStencilValidCF(pairSten.second, dit());
	}
      }

      // If we couldn't find in a quadrant, try a larger neighborhood
      order = m_order;
      while(!foundStencil && order > 0){
      	foundStencil = this->getLeastSquaresStencil(pairSten, vof, VofUtils::Neighborhood::Radius, dit(), order);
      	order--;

	// Check if stencil reaches too far across CF
	if(foundStencil) {
	  foundStencil = this->isStencilValidCF(pairSten.second, dit());
	}
      }

      if(foundStencil){
	weights (vof, m_comp) = pairSten.first;
	stencils(vof, m_comp) = pairSten.second;
	
	// Stencil and weight must also be scaled by the B-coefficient, dx (because it's used in kappa*Div(F)) and the area fraction. 
	weights (vof, m_comp) *= B*areaFrac/m_dx;
	stencils(vof, m_comp) *= B*areaFrac/m_dx;
      }
      else{
	// Dead cell. No flux. 
	weights (vof, m_comp) = 0.0;
	stencils(vof, m_comp).clear();
      }
    }
  }
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
    const RealVect pos = this->getBoundaryPosition(a_vof, a_dit);
    value = m_functionValue(pos);
  }

  // B-coefficient, area fraction, and division by dx (from Div(F)) already a part of the boundary weights, but
  // beta is not. 
  a_Lphi(a_vof, m_comp) += a_beta*value*m_boundaryWeights[a_dit](a_vof, m_comp);
  
  return;
}

bool EBHelmholtzDirichletEBBC::getLeastSquaresStencil(std::pair<Real, VoFStencil>& a_stencil,
						      const VolIndex&              a_vof,
						      const VofUtils::Neighborhood a_neighborhood,
						      const DataIndex&             a_dit,
						      const int                    a_order) const {
  bool foundStencil = false;
  
  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
  const RealVect normal  = ebisbox.normal(a_vof);  
    
  const VoFStencil gradientStencil = LeastSquares::getBndryGradSten(a_vof, a_neighborhood, ebisbox, m_dx, a_order, m_weight, a_order);

  if(gradientStencil.size() > 0 && normal != RealVect::Zero){
    
    const VoFStencil DphiDnStencil =  LeastSquares::projectGradSten(gradientStencil, -normal);
    const Real boundaryWeight      = -LeastSquares::sumAllWeights(DphiDnStencil);

    a_stencil = std::make_pair(boundaryWeight, DphiDnStencil);

    foundStencil = true;
  }

  return foundStencil;
}

bool EBHelmholtzDirichletEBBC::getJohansenStencil(std::pair<Real, VoFStencil>& a_stencil, const VolIndex& a_vof, const DataIndex& a_dit, const int a_order) const {
  if(!(a_order == 1 || a_order == 2)) MayDay::Abort("EBHelmholtzDirichletEBBC::getJohansenStencil - only order 1 and 2 is supported!");
  
  bool foundStencil = false;

  const EBISBox& ebisbox  = m_eblg.getEBISL()[a_dit];
  const RealVect normal   = ebisbox.normal(a_vof);
  const RealVect centroid = ebisbox.bndryCentroid(a_vof);

  bool dropOrder;
  Vector<VoFStencil> pointStencils;
  Vector<Real>       distanceAlongLine;
  EBArith::dataRayCast(dropOrder, pointStencils, distanceAlongLine, normal, centroid, a_vof, ebisbox, m_dx*RealVect::Unit, IntVectSet(), m_comp, a_order);

  if(!dropOrder){
    Real       weight;
    VoFStencil stencil;
    
    if(a_order == 1){
      const Real x1 = distanceAlongLine[0];
      const Real di = 1./x1;

      weight   = 1./x1;
      stencil  = pointStencils[0];
      stencil *= -1./x1;
    }
    else if(a_order == 2){
      const Real x1    = distanceAlongLine[0];
      const Real x2    = distanceAlongLine[1];
      const Real denom = x2*x2*x1 - x1*x1*x2;

      VoFStencil phi1Sten = pointStencils[0];
      VoFStencil phi2Sten = pointStencils[1];

      phi1Sten *= -x2*x2/denom;
      phi2Sten *=  x1*x1/denom;

      weight   = -x1*x1/denom + x2*x2/denom;
      stencil += phi1Sten;
      stencil += phi2Sten;
    }

    a_stencil = std::make_pair(weight, stencil);

    foundStencil = true;
  }

  return foundStencil;
}

bool EBHelmholtzDirichletEBBC::getChomboLsqStencil(std::pair<Real, VoFStencil>& a_stencil, const VolIndex& a_vof, const DataIndex& a_dit) const {
  VoFStencil sten;
  Real bndryWeight;

  EBArith::getLeastSquaresGradSten(sten, bndryWeight, a_vof, m_eblg.getEBISL()[a_dit], m_dx*RealVect::Unit, m_eblg.getDomain(), m_comp);

  a_stencil = std::make_pair(bndryWeight, sten);

  return (sten.size() > 0);
}

#include <CD_NamespaceFooter.H>
