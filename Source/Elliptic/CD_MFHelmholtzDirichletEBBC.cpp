/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzDirichletEBBC.cpp
  @brief  Implementation of CD_MFHelmholtzDirichletEBBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzDirichletEBBC.H>
#include <CD_LeastSquares.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzDirichletEBBC::MFHelmholtzDirichletEBBC(const int a_phase, const RefCountedPtr<JumpBC>& a_jumpBC) : MFHelmholtzEBBC(a_phase, a_jumpBC) {
  m_order       = -1;
  m_weight      = -1;
  m_useConstant = false;
  m_useFunction = false;
}

MFHelmholtzDirichletEBBC::~MFHelmholtzDirichletEBBC(){

}

void MFHelmholtzDirichletEBBC::setOrder(const int a_order){
  CH_assert(a_order > 0);
  m_order = a_order;
}

void MFHelmholtzDirichletEBBC::setWeight(const int a_weight){
  CH_assert(a_weight > 0);
  m_weight = a_weight;
}

void MFHelmholtzDirichletEBBC::setValue(const int a_value){
  m_useConstant = true;
  m_useFunction = false;
  
  m_constantValue = a_value;
}

void MFHelmholtzDirichletEBBC::setValue(const std::function<Real(const RealVect& a_pos)>& a_value){
  m_useConstant = false;
  m_useFunction = true;
  
  m_functionValue = a_value;
}

void MFHelmholtzDirichletEBBC::define() {
  if(  m_order <= 0  || m_weight <= 0 ) MayDay::Error("MFHelmholtzDirichletEBBC - must have order > 0 and weight > 0");
  if(!(m_useConstant || m_useFunction)) MayDay::Error("MFHelmholtzDirichletEBBC - not using constant or function!");

  const DisjointBoxLayout& dbl = m_eblg.getDBL();
  const ProblemDomain& domain  = m_eblg.getDomain();
  
  m_boundaryWeights.define(dbl);
  m_kappaDivFStencils.define(dbl);
  m_gradStencils.define(dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box box          = dbl[dit()];
    const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs  = ebisbox.getIrregIVS(box);

    BaseIVFAB<Real>&       weights  = m_boundaryWeights  [dit()];
    BaseIVFAB<VoFStencil>& stencils = m_gradStencils[dit()];

    const BaseIVFAB<Real>& Bcoef    = (*m_Bcoef)[dit()];

    weights. define(ivs, ebgraph, m_nComp);
    stencils.define(ivs, ebgraph, m_nComp);
    m_kappaDivFStencils[dit()].define(ivs, ebgraph, m_nComp); // Left undefined because full flux is applied in applyEBFlux

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const Real areaFrac = ebisbox.bndryArea(vof);
      const Real B        = Bcoef(vof, m_comp);

      int order;
      bool foundStencil = false;
      std::pair<Real, VoFStencil> pairSten;

      // Try quadrants first. 
      order = m_isMGLevel ? 1 : m_order; 
      while(!foundStencil && order > 0){
      	foundStencil = this->getLeastSquaresBoundaryGradStencil(pairSten, vof, VofUtils::Neighborhood::Quadrant, dit(), order);
      	order--;

	// Check if stencil reaches too far across CF
	if(foundStencil) {
	  foundStencil = this->isStencilValidCF(pairSten.second, dit());
	}
      }

      // If we couldn't find in a quadrant, try a larger neighborhood
      order = m_isMGLevel ? 1 : m_order;      
      while(!foundStencil && order > 0){
      	foundStencil = this->getLeastSquaresBoundaryGradStencil(pairSten, vof, VofUtils::Neighborhood::Radius, dit(), order);
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

void MFHelmholtzDirichletEBBC::applyEBFlux(VoFIterator&       a_singlePhaseVofs,
					   VoFIterator&       a_multiPhaseVofs,
					   EBCellFAB&         a_Lphi,
					   const EBCellFAB&   a_phi,
					   const DataIndex&   a_dit,
					   const Real&        a_beta,
					   const bool&        a_homogeneousPhysBC) const {
  // Apply the stencil for computing the contribution to kappaDivF. Note divF is sum(faces) B*grad(Phi)/dx and that this
  // is the contribution from the EB face. B/dx is already included in the stencils and boundary weights, but beta is not.

  // Do single phase cells
  for(a_singlePhaseVofs.reset(); a_singlePhaseVofs.ok(); ++a_singlePhaseVofs){
    const VolIndex& vof = a_singlePhaseVofs();

    // Homogeneous contribution
    a_Lphi(vof, m_comp) += a_beta*this->applyStencil(m_gradStencils[a_dit](vof, m_comp), a_phi);

    // Inhomogeneous contribution. 
    if(!a_homogeneousPhysBC){
      Real value = 0.0;
    
      if(m_useConstant){
	value = m_constantValue;
      }
      else if(m_useFunction){
	const RealVect pos = this->getBoundaryPosition(vof, a_dit);
	value = m_functionValue(pos);
      }

      a_Lphi(vof, m_comp) += a_beta*value*m_boundaryWeights[a_dit](vof, m_comp);
    }
  }

  // Do multi-phasephase cells. Currently just a copy of the other.
  for(a_multiPhaseVofs.reset(); a_multiPhaseVofs.ok(); ++a_multiPhaseVofs){
    const VolIndex& vof = a_multiPhaseVofs();

    // Homogeneous contribution
    a_Lphi(vof, m_comp) += a_beta*this->applyStencil(m_gradStencils[a_dit](vof, m_comp), a_phi);

    // Inhomogeneous contribution. 
    if(!a_homogeneousPhysBC){
      Real value = 0.0;
    
      if(m_useConstant){
	value = m_constantValue;
      }
      else if(m_useFunction){
	const RealVect pos = this->getBoundaryPosition(vof, a_dit);
	value = m_functionValue(pos);
      }

      value = m_jumpBC->getBndryPhi(m_phase, a_dit)(vof, m_comp);

      a_Lphi(vof, m_comp) += a_beta*value*m_boundaryWeights[a_dit](vof, m_comp);
    }
  }

  
  
  return;
}

bool MFHelmholtzDirichletEBBC::getLeastSquaresBoundaryGradStencil(std::pair<Real, VoFStencil>& a_stencil,
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

#include <CD_NamespaceFooter.H>
