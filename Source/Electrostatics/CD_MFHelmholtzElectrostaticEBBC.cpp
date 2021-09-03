/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzElectrostaticEBBC.cpp
  @brief  Implementation of CD_MFHelmholtzElectrostaticEBBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzElectrostaticEBBC.H>
#include <CD_LeastSquares.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzElectrostaticEBBC::MFHelmholtzElectrostaticEBBC(const int a_phase, const ElectrostaticEbBc& a_electrostaticBCs, const RefCountedPtr<JumpBC>& a_jumpBC)
  : MFHelmholtzEBBC(a_phase, a_jumpBC) {
  CH_TIME("MFHelmholtzElectrostaticEBBC::MFHelmholtzElectrostaticEBBC");

  m_order  = -1;
  m_weight = -1;

  m_electrostaticBCs = a_electrostaticBCs;
}

MFHelmholtzElectrostaticEBBC::~MFHelmholtzElectrostaticEBBC(){
  CH_TIME("MFHelmholtzElectrostaticEBBC::~MFHelmholtzElectrostaticEBBC");
}

void MFHelmholtzElectrostaticEBBC::setOrder(const int a_order){
  CH_TIME("MFHelmholtzElectrostaticEBBC::setOrder(int)");
  
  CH_assert(a_order > 0);

  m_order = a_order;
}

void MFHelmholtzElectrostaticEBBC::setWeight(const int a_weight){
  CH_TIME("MFHelmholtzElectrostaticEBBC::setWeight(int)");
  
  CH_assert(a_weight >= 0);

  m_weight = a_weight;
}

void MFHelmholtzElectrostaticEBBC::defineSinglePhase() {
  CH_TIME("MFHelmholtzElectrostaticEBBC::defineSinglePhase()");

  CH_assert(m_order  >  0);
  CH_assert(m_weight >= 0);  

  const DisjointBoxLayout& dbl = m_eblg.getDBL();
  const ProblemDomain& domain  = m_eblg.getDomain();

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box box          = dbl[dit()];
    const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs  = ebisbox.getIrregIVS(box);

    BaseIVFAB<Real>&       weights  = m_boundaryWeights [dit()];
    BaseIVFAB<VoFStencil>& stencils = m_kappaDivFStencils[dit()];
    
    const BaseIVFAB<Real>& Bcoef    = (*m_Bcoef)[dit()];

    VoFIterator& singlePhaseVofs = m_jumpBC->getSinglePhaseVofs(m_phase, dit());

    for (singlePhaseVofs.reset(); singlePhaseVofs.ok(); ++singlePhaseVofs){
      const VolIndex& vof = singlePhaseVofs();
      const Real areaFrac = ebisbox.bndryArea(vof);
      const Real B        = Bcoef(vof, m_comp);

      int order;
      bool foundStencil = false;
      std::pair<Real, VoFStencil> pairSten;

      // Try quadrants first.
      order = m_order;
      while(!foundStencil && order > 0){
      	foundStencil = this->getLeastSquaresBoundaryGradStencil(pairSten, vof, VofUtils::Neighborhood::Quadrant, dit(), order, m_weight);
      	order--;

	// Check if stencil reaches too far across CF
	if(foundStencil) {
	  foundStencil = this->isStencilValidCF(pairSten.second, dit());
	}
      }

      // If we couldn't find in a quadrant, try a larger neighborhood
      order = m_order;
      while(!foundStencil && order > 0){
      	foundStencil = this->getLeastSquaresBoundaryGradStencil(pairSten, vof, VofUtils::Neighborhood::Radius, dit(), order, m_weight);
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

void MFHelmholtzElectrostaticEBBC::applyEBFluxSinglePhase(VoFIterator&       a_singlePhaseVofs,
							  EBCellFAB&         a_Lphi,
							  const EBCellFAB&   a_phi,
							  const DataIndex&   a_dit,
							  const Real&        a_beta,
							  const bool&        a_homogeneousPhysBC) const {
  CH_TIME("MFHelmholtzElectrostaticEBBC::applyEBFluxSinglePhase(VoFIterator, EBCellFAB, EBCellFAB, DataIndex, Real, bool)");
  
  // Apply the stencil for computing the contribution to kappaDivF. Note divF is sum(faces) B*grad(Phi)/dx and that this
  // is the contribution from the EB face. B/dx is already included in the stencils and boundary weights, but beta is not.

  // Do single phase cells
  if(!a_homogeneousPhysBC){  
    for(a_singlePhaseVofs.reset(); a_singlePhaseVofs.ok(); ++a_singlePhaseVofs){
      const VolIndex& vof  = a_singlePhaseVofs();
      
      const RealVect pos   = this->getBoundaryPosition(vof, a_dit);      
      const Real     value = this->getElectrodePotential(pos);      

      a_Lphi(vof, m_comp) += a_beta*value*m_boundaryWeights[a_dit](vof, m_comp);
    }
  }
  
  return;
}

#include <CD_NamespaceFooter.H>
