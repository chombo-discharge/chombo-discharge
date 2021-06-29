/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzRobinEBBC.cpp
  @brief  Implementation of CD_MFHelmholtzRobinEBBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_MFHelmholtzRobinEBBC.H>
#include <CD_VofUtils.H>
#include <CD_LeastSquares.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzRobinEBBC::MFHelmholtzRobinEBBC(const int a_phase, const RefCountedPtr<JumpBC>& a_jumpBC) : MFHelmholtzEBBC(a_phase, a_jumpBC) {
  m_order       = -1;
  m_weight      = -1;
  m_useConstant = false;
  m_useFunction = false;
}

MFHelmholtzRobinEBBC::~MFHelmholtzRobinEBBC(){

}

void MFHelmholtzRobinEBBC::setCoefficients(const Real a_A, const Real a_B, const Real a_C){
  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
}

void MFHelmholtzRobinEBBC::setCoefficients(const std::function<Real(const RealVect& a_pos) >& a_A,
					   const std::function<Real(const RealVect& a_pos) >& a_B,
					   const std::function<Real(const RealVect& a_pos) >& a_C){
  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;
}

VoFStencil MFHelmholtzRobinEBBC::getInterpolationStencil(const VolIndex& a_vof, const DataIndex& a_dit) const {
  VoFStencil stencil;

  // First, try direct interpolation using least squares. Then, try Chombo taylor extrapolation. 
  // I don't think it should be necessary to check if the stencils reach across the CF interface because
  // MFHelmholtzOp always interpolates at least one ghost cell, and these stencils should really have a width of 1. 
  if(stencil.size() == 0) stencil = this->getMonoPathStencil(a_vof, a_dit);
  if(stencil.size() == 0) stencil = this->getTaylorStencil  (a_vof, a_dit);
  
  return stencil;
}

VoFStencil MFHelmholtzRobinEBBC::getMonoPathStencil(const VolIndex& a_vof, const DataIndex& a_dit) const {
  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
  const int pow          = 0;  // Can't run with weighting when the starting vof is included
  const int order        = 1;
  const int radius       = 1;
  const bool useStartVof = true;

  const VoFStencil stencil = LeastSquares::getInterpolationStencil(Location::Cell::Boundary,
								   Location::Cell::Center,
								   LeastSquares::Connectivity::MonotonePath,
								   a_vof,
								   ebisbox,
								   m_dx,
								   pow,
								   radius,
								   order,
								   useStartVof);


  return stencil;
}

VoFStencil MFHelmholtzRobinEBBC::getTaylorStencil(const VolIndex& a_vof, const DataIndex& a_dit) const {
  VoFStencil stencil;
  
  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
  const RealVect dist    = ebisbox.bndryCentroid(a_vof)*m_dx;
  IntVectSet ivs = IntVectSet();

  EBArith::getFirstOrderExtrapolationStencil(stencil, dist, m_dx*RealVect::Unit, a_vof, ebisbox, -1, &ivs, m_comp);

  return stencil;
}
  
void MFHelmholtzRobinEBBC::defineSinglePhase() {
  if(  m_order <= 0  || m_weight <= 0 ) MayDay::Error("MFHelmholtzRobinEBBC - must have order > 0 and weight > 0");
  if(!(m_useConstant || m_useFunction)) MayDay::Error("MFHelmholtzRobinEBBC - not using constant or function!");

  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box box          = dbl[dit()];
    const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs  = ebisbox.getIrregIVS(box);

    BaseIVFAB<Real>&       weights  = m_boundaryWeights [dit()];
    BaseIVFAB<VoFStencil>& stencils = m_kappaDivFStencils[dit()];
    

    // Build interpolation stencils. 
    VoFIterator& singlePhaseVofs = m_jumpBC->getSinglePhaseVofs(m_phase, dit());
    for (singlePhaseVofs.reset(); singlePhaseVofs.ok(); ++singlePhaseVofs){
      const VolIndex& vof = singlePhaseVofs();
      const Real areaFrac = ebisbox.bndryArea(vof);
      const Real helmBco  = (*m_Bcoef)[dit()](vof, m_comp);

      weights(vof,  m_comp) = 0.0;
      stencils(vof, m_comp) = this->getInterpolationStencil(vof, dit());

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

void MFHelmholtzRobinEBBC::applyEBFluxSinglePhase(VoFIterator&       a_singlePhaseVofs,
						  EBCellFAB&         a_Lphi,
						  const EBCellFAB&   a_phi,
						  const DataIndex&   a_dit,
						  const Real&        a_beta,
						  const bool&        a_homogeneousPhysBC) const {

  // TLDR: For Robin, the flux is b*dphi/dn = beta*b*A*phi/B - beta*b*C/B and we have stored
  //       the term b*A*phi/B in the interpolation stencil. The other term we compute below. 
  for (a_singlePhaseVofs.reset(); a_singlePhaseVofs.ok(); ++a_singlePhaseVofs){
    const VolIndex& vof = a_singlePhaseVofs();

    // Homogeneous contribution
    //    a_Lphi(vof, m_comp) += a_beta*this->applyStencil(m_interpolationStencils[a_dit](vof, m_comp), a_phi);

    // Inhomogeneous contribution
    if(!a_homogeneousPhysBC){
      Real B;
      Real C;
      if(m_useConstant){
	B = m_constantB;
	C = m_constantC;
      }
      else if(m_useFunction){
	const RealVect pos = this->getBoundaryPosition(vof, a_dit);
	B = m_functionB(pos);
	C = m_functionC(pos);
      }

      const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
      const Real areaFrac    = ebisbox.bndryArea(vof);
      const Real helmBco     = (*m_Bcoef)[a_dit](vof, m_comp);
      const Real kappaDivF   = -a_beta*helmBco*areaFrac*C/(m_dx*B);

      if(std::abs(B) > 0.0){
	a_Lphi(vof, m_comp) += kappaDivF;
      }
    }
  }

  
  return;
}

#include <CD_NamespaceFooter.H>
