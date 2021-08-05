/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzEBBC.cpp
  @brief  Implementation of CD_MFHelmholtzEBBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzEBBC.H>
#include <CD_LeastSquares.H>
#include <CD_NamespaceHeader.H>


MFHelmholtzEBBC::MFHelmholtzEBBC(const int a_phase, const RefCountedPtr<JumpBC>& a_jumpBC){
  m_jumpBC = a_jumpBC;
  m_phase  = a_phase;
}

MFHelmholtzEBBC::~MFHelmholtzEBBC(){

}

void MFHelmholtzEBBC::setOrder(const int a_order){
  CH_assert(a_order > 0);
  m_order = a_order;
}

void MFHelmholtzEBBC::setWeight(const int a_weight){
  CH_assert(a_weight > 0);
  m_weight = a_weight;
}

void MFHelmholtzEBBC::define(){
  this->defineMultiPhase();
  this->defineSinglePhase();
}

void MFHelmholtzEBBC::defineMultiPhase(){
  const DisjointBoxLayout& dbl = m_eblg.getDBL();
  const ProblemDomain& domain  = m_eblg.getDomain();

  m_boundaryWeights.  define(dbl);
  m_kappaDivFStencils.define(dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box box                = dbl[dit()];
    const EBISBox& ebisbox       = m_eblg.getEBISL()[dit()];
    const EBGraph& ebgraph       = ebisbox.getEBGraph();
    const IntVectSet& ivs        = ebisbox.getIrregIVS(box);
    const BaseIVFAB<Real>& Bcoef = (*m_Bcoef)[dit()];

    // These are used to reconstruct gradients at the EB. They are left undefined for single-phase cells. 
    m_boundaryWeights[dit()]  .define(ivs, ebgraph, m_nComp); 
    m_kappaDivFStencils[dit()].define(ivs, ebgraph, m_nComp); 
    
    BaseIVFAB<Real>&       weights  = m_boundaryWeights  [dit()];
    BaseIVFAB<VoFStencil>& stencils = m_kappaDivFStencils[dit()];

    VoFIterator& multiPhaseVofs  = m_jumpBC->getMultiPhaseVofs(m_phase, dit());

    // Define over the full ivs but only fill on multi-phase cells. 
    weights. define(ivs, ebgraph, m_nComp);
    stencils.define(ivs, ebgraph, m_nComp);

    // Build stencils for each vof. The order for the multiphase VoFs should follow the order for jumpBC, I think. 
    for (multiPhaseVofs.reset(); multiPhaseVofs.ok(); ++multiPhaseVofs){
      const VolIndex& vof = multiPhaseVofs();
      const Real areaFrac = ebisbox.bndryArea(vof);
      const Real B        = Bcoef(vof, m_comp);

      int order;
      bool foundStencil = false;
      std::pair<Real, VoFStencil> pairSten;

      // Try quadrants first. 
      order = m_jumpBC->getOrder();
      while(!foundStencil && order > 0){
      	foundStencil = this->getLeastSquaresBoundaryGradStencil(pairSten, vof, VofUtils::Neighborhood::Quadrant, dit(), order);
      	order--;

	// Check if stencil reaches too far across CF
	if(foundStencil) {
	  foundStencil = this->isStencilValidCF(pairSten.second, dit());
	}
      }

      // If we couldn't find in a quadrant, try a larger neighborhood
      order = m_jumpBC->getOrder();      
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
  
void MFHelmholtzEBBC::applyEBFlux(VoFIterator&       a_vofit,
				  EBCellFAB&         a_Lphi,
				  const EBCellFAB&   a_phi,
				  const DataIndex&   a_dit,
				  const Real&        a_beta,
				  const bool&        a_homogeneousPhysBC) const {
  
  VoFIterator& singlePhaseVofs = m_jumpBC->getSinglePhaseVofs(m_phase, a_dit);
  VoFIterator& multiPhaseVofs  = m_jumpBC->getMultiPhaseVofs (m_phase, a_dit);

  this->applyEBFluxSinglePhase(singlePhaseVofs, a_Lphi, a_phi, a_dit, a_beta, a_homogeneousPhysBC);
  this->applyEBFluxMultiPhase (multiPhaseVofs,  a_Lphi, a_phi, a_dit, a_beta, a_homogeneousPhysBC);
}

void MFHelmholtzEBBC::applyEBFluxMultiPhase(VoFIterator&       a_multiPhaseVofs,
					    EBCellFAB&         a_Lphi,
					    const EBCellFAB&   a_phi,
					    const DataIndex&   a_dit,
					    const Real&        a_beta,
					    const bool&        a_homogeneousPhysBC) const {
  // Apply the stencil for computing the contribution to kappaDivF. Note divF is sum(faces) B*grad(Phi)/dx and that this
  // is the contribution from the EB face. B/dx is already included in the stencils and boundary weights, but beta is not.
  for(a_multiPhaseVofs.reset(); a_multiPhaseVofs.ok(); ++a_multiPhaseVofs){
    const VolIndex& vof = a_multiPhaseVofs();

    // Homogeneous contribution
    const Real phiB = m_jumpBC->getBndryPhi(m_phase, a_dit)(vof, m_comp);
    a_Lphi(vof, m_comp) += a_beta*phiB*m_boundaryWeights[a_dit](vof, m_comp);
  }
  
  return;
}

bool MFHelmholtzEBBC::getLeastSquaresBoundaryGradStencil(std::pair<Real, VoFStencil>& a_stencil,
							 const VolIndex&              a_vof,
							 const VofUtils::Neighborhood a_neighborhood,
							 const DataIndex&             a_dit,
							 const int                    a_order) const {
  bool foundStencil = false;
  
  const bool addStartVof = false;
  
  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
  const RealVect normal  = ebisbox.normal(a_vof);  
    
  const VoFStencil gradStencil = LeastSquares::getBndryGradSten(a_vof,
								a_neighborhood,
								m_dataLocation,
								ebisbox,
								m_dx,
								a_order,
								m_weight,
								a_order,
								addStartVof);

  if(gradStencil.size() > 0 && normal != RealVect::Zero){
    
    const VoFStencil DphiDnStencil =  LeastSquares::projectGradSten(gradStencil, -normal);
    const Real boundaryWeight      = -LeastSquares::sumAllWeights(DphiDnStencil);

    a_stencil = std::make_pair(boundaryWeight, DphiDnStencil);

    foundStencil = true;
  }

  return foundStencil;
}

#include <CD_NamespaceFooter.H>
