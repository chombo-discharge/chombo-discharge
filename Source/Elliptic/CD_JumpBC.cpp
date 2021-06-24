/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_JumpBC.cpp
  @brief  Implementatino of CD_JumpBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_JumpBC.H>
#include <CD_VofUtils.H>
#include <CD_LeastSquares.H>
#include <CD_NamespaceHeader.H>

constexpr int JumpBC::m_comp;
constexpr int JumpBC::m_nComp;

JumpBC::JumpBC(const MFLevelGrid& a_mflg,
	       const BcoefPtr&    a_Bcoef,
	       const Real         a_dx,
	       const int          a_order,
	       const int          a_weight,
	       const int          a_radius,
	       const int          a_ghostCF){
  CH_TIME("JumpBC::JumpBC");
  
  m_mflg      = a_mflg;
  m_Bcoef     = a_Bcoef;
  m_dx        = a_dx;
  m_weight    = a_weight;
  m_order     = a_order;
  m_radius    = a_radius;
  m_ghostCF   = a_ghostCF;
  m_numPhases = m_mflg.numPhases();
  m_isMGLevel = false;

  this->define();
}

JumpBC::~JumpBC(){

}

void JumpBC::setMG(const bool a_isMGLevel){
  m_isMGLevel = a_isMGLevel;
}

void JumpBC::define(){
  if(m_numPhases > 1) { // JumpBC will never be called unless it's a multiphase problem.
    
    const DisjointBoxLayout& dbl = m_mflg.getGrids();

    m_boundaryPhi.define(dbl);
    m_avgStencils.define(dbl);
    m_avgWeights.define(dbl);
    m_avgBco.define(dbl);

    for (DataIterator dit(dbl); dit.ok(); ++dit){

      m_boundaryPhi[dit()].define(m_mflg, dit());
      m_avgStencils[dit()].define(m_mflg, dit());
      m_avgWeights [dit()].define(m_mflg, dit());
      m_avgBco     [dit()].define(m_mflg, dit());
    }

    // Build stencils and weights for each phase. 
    for (int iphase = 0; iphase < m_numPhases; iphase++){
      const EBLevelGrid& eblg = m_mflg.getEBLevelGrid(iphase);
      const EBISLayout& ebisl = eblg.getEBISL();

      for (DataIterator dit(dbl); dit.ok(); ++dit){
	const Box box          = dbl[dit()];
	const EBISBox& ebisbox = ebisl[dit()];
	const EBGraph& ebgraph = ebisbox.getEBGraph();

	// Build the average stencils. Only matters if the cell is a multi-valued cell. 
	BaseIVFAB<VoFStencil>& avgStencils = m_avgStencils[dit()].getIVFAB(iphase);
	BaseIVFAB<Real>& avgWeights        = m_avgWeights [dit()].getIVFAB(iphase);
	BaseIVFAB<Real>& avgBco            = m_avgBco     [dit()].getIVFAB(iphase);
	const BaseIVFAB<Real>& Bcoef       = (*m_Bcoef)   [dit()].getIVFAB(iphase);

	// Matching grid cells. 
	const IntVectSet& ivs = avgWeights.getIVS();

	// Build stencils like we always do
	BaseIVFAB<VoFStencil> gradStencils(ivs, ebgraph, m_nComp);
	BaseIVFAB<Real>       bndryWeights(ivs, ebgraph, m_nComp);
	for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();

	  int order;
	  bool foundStencil = false;
	  std::pair<Real, VoFStencil> pairSten;
	  
	  // Try quadrants first. 
	  order = m_isMGLevel ? 1 : m_order;
	  while(!foundStencil && order > 0){
	    foundStencil = this->getLeastSquaresBoundaryGradStencil(pairSten, vof, ebisbox, VofUtils::Neighborhood::Quadrant, order);
	    order--;

	    // Check if stencil reaches too far across CF
	    if(foundStencil) {
	      foundStencil = this->isStencilValidCF(pairSten.second, dit());
	    }
	  }

	  // If we couldn't find in a quadrant, try a larger neighborhood
	  order = m_isMGLevel ? 1 : m_order;
	  while(!foundStencil && order > 0){
	    foundStencil = this->getLeastSquaresBoundaryGradStencil(pairSten, vof, ebisbox, VofUtils::Neighborhood::Radius, order);
	    order--;

	    // Check if stencil reaches too far across CF
	    if(foundStencil) {
	      foundStencil = this->isStencilValidCF(pairSten.second, dit());
	    }
	  }

	  if(foundStencil){
	    bndryWeights(vof, m_comp) = pairSten.first;
	    gradStencils(vof, m_comp) = pairSten.second;
	  }
	  else{
	    bndryWeights(vof, m_comp) = 0.0;
	    gradStencils(vof, m_comp).clear();
	  }
	}

	// Build average stencils
	for (IVSIterator ivsIt(ivs); ivsIt.ok(); ++ivsIt){
	  const IntVect iv = ivsIt();

	  const VolIndex curVof(iv, 0);
	  const Vector<VolIndex> allVofs = ebisbox.getVoFs(iv);

	  Real& curBco           = avgBco     (curVof, m_comp);
	  Real& curWeight        = avgWeights (curVof, m_comp);
	  VoFStencil& curStencil = avgStencils(curVof, m_comp);

	  curBco    = 0.0;
	  curWeight = 0.0;
	  curStencil.clear();

	  for (int ivof = 0; ivof < allVofs.size(); ivof++){
	    const VolIndex& vof = allVofs[ivof];

	    curBco     += Bcoef       (vof, m_comp);
	    curWeight  += bndryWeights(vof, m_comp);
	    curStencil += gradStencils(vof, m_comp);
	  }

	  const Real invNum = 1./allVofs.size();

	  curBco     *= invNum;
	  curWeight  *= invNum;
	  curStencil *= invNum;
	}

      }
    }

    this->resetBC();
  }
}

const BaseIVFAB<Real>& JumpBC::getBndryPhi(const int a_phase, const DataIndex& a_dit) const {
  return m_boundaryPhi[a_dit].getIVFAB(a_phase);
}

bool JumpBC::getLeastSquaresBoundaryGradStencil(std::pair<Real, VoFStencil>& a_stencil,
						const VolIndex&              a_vof,
						const EBISBox&               a_ebisbox,
						const VofUtils::Neighborhood a_neighborhood,
						const int                    a_order) const {
  bool foundStencil = false;

  const RealVect normal   = a_ebisbox.normal(a_vof);  
    
  const VoFStencil gradientStencil = LeastSquares::getBndryGradSten(a_vof, a_neighborhood, a_ebisbox, m_dx, a_order, m_weight, a_order);

  if(gradientStencil.size() > 0 && normal != RealVect::Zero){
    
    const VoFStencil DphiDnStencil =  LeastSquares::projectGradSten(gradientStencil, -normal);
    const Real boundaryWeight      = -LeastSquares::sumAllWeights(DphiDnStencil);

    a_stencil = std::make_pair(boundaryWeight, DphiDnStencil);

    foundStencil = true;
  }

  return foundStencil;
}

void JumpBC::resetBC() const {
  for (DataIterator dit = m_boundaryPhi.dataIterator(); dit.ok(); ++dit){
    for (int iphase = 0; iphase < m_numPhases; iphase++){
      BaseIVFAB<Real>& bndryPhi = m_boundaryPhi[dit()].getIVFAB(iphase);
      bndryPhi.setVal(0.0);
    }
  }
}

void JumpBC::matchBC(const LevelData<MFCellFAB>&        a_phi,
		     const LevelData<BaseIVFAB<Real> >& a_jump,
		     const bool                         a_homogeneousPhysBC) const {
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit){
    this->matchBC(a_phi[dit()], a_jump[dit()], a_homogeneousPhysBC, dit());
  }
}

#include <CD_NamespaceFooter.H>
