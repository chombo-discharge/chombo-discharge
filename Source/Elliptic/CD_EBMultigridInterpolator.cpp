/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBMultigridInterpolator.cpp
  @brief  Implementation of CD_EBMultigridInterpolator.H
  @author Robert Marskar
  @todo   Consider making a completely new layout near the EB, so we don't define the temp holders for the full domain. 
  @todo   When we do the interpolation, we should copy from coarse to m_coarData and from fine to m_fineData and then run everything locally on each patch. 
  @todo   Check the required number of ghost cells -- not sure if 2*ghostCF is correct.
  @todo   Think about what we need to do if we have multi-cells. We should probably only include coarse-grid cells that can be reached with a monotone path from the fine vof. 
  @todo   defineGhosts should only do this near the EB
*/

// Chombo includes
#include <EBCellFactory.H>
#include <EBLevelDataOps.H>
#include <NeighborIterator.H>
#include <InterpF_F.H>

// Our includes
#include <CD_EBMultigridInterpolator.H>
#include <CD_VofUtils.H>
#include <CD_LeastSquares.H>
#include <CD_NamespaceHeader.H>

constexpr int EBMultigridInterpolator::m_stenComp;
constexpr int EBMultigridInterpolator::m_nStenComp;

EBMultigridInterpolator::EBMultigridInterpolator(){
  m_isDefined    = false;
  m_order        = 2;
  m_weight       = 2;
  m_cellLocation = Location::Cell::Center;
}

EBMultigridInterpolator::EBMultigridInterpolator(const EBLevelGrid& a_eblgFine,
						 const EBLevelGrid& a_eblgCoar,
						 const int          a_refRat,
						 const int          a_nVar,
						 const int          a_ghostCF) : EBMultigridInterpolator() {
  CH_assert(a_ghostCF   > 0);
  CH_assert(a_nVar      > 0);
  CH_assert(a_refRat%2 == 0);

  // Build the CFIVS for EBQuadCFInterp. 
  const DisjointBoxLayout& dblFine = a_eblgFine.getDBL();
  const ProblemDomain& domainFine  = a_eblgFine.getDomain();
  
  LayoutData<IntVectSet> cfivs(dblFine);
  for (DataIterator dit(dblFine); dit.ok(); ++dit){
    Box grownBox = grow(dblFine[dit()], 1);
    grownBox &= domainFine;

    cfivs[dit()] = IntVectSet(grownBox);
    
    NeighborIterator nit(dblFine);
    for (nit.begin(dit()); nit.ok(); ++ nit){
      cfivs[dit()] -= dblFine[nit()];
    }
  }

  EBQuadCFInterp::define(a_eblgFine.getDBL(),
			 a_eblgCoar.getDBL(),
			 a_eblgFine.getEBISL(),
			 a_eblgCoar.getEBISL(),
			 a_eblgCoar.getDomain(),
			 a_refRat,
			 a_nVar,
			 cfivs,
			 a_eblgFine.getEBIS(),
			 true);
  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;
  this->defineCFIVS();

// NEW CONSTRUCTOR IMPLEMENTATION GOES BELOW HERE
  m_refRat   = a_refRat;
  m_nComp    = a_nVar;
  m_ghostCF  = a_ghostCF;
  
  this->defineGrids(a_eblgFine, a_eblgCoar);
  this->defineGhosts();
  this->defineData();
  this->defineStencils();
}

int EBMultigridInterpolator::getGhostCF() const{
  return m_ghostCF;
}

void EBMultigridInterpolator::defineCFIVS(){

  
  for (int dir = 0; dir < SpaceDim; dir++){
    m_loCFIVS[dir].define(m_eblgFine.getDBL());
    m_hiCFIVS[dir].define(m_eblgFine.getDBL());

    for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit){
      m_loCFIVS[dir][dit()].define(m_eblgFine.getDomain(), m_eblgFine.getDBL()[dit()], m_eblgFine.getDBL(), dir, Side::Lo);
      m_hiCFIVS[dir][dit()].define(m_eblgFine.getDomain(), m_eblgFine.getDBL()[dit()], m_eblgFine.getDBL(), dir, Side::Hi);
    }
  }
}

EBMultigridInterpolator::~EBMultigridInterpolator(){

}

void EBMultigridInterpolator::coarseFineInterp(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phiCoar, const Interval a_variables){
  EBQuadCFInterp::interpolate(a_phiFine, a_phiCoar, a_variables);
}

void EBMultigridInterpolator::coarseFineInterpH(LevelData<EBCellFAB>& a_phiFine, const Interval a_variables){
  for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit){
    this->coarseFineInterpH(a_phiFine[dit()], a_variables, dit());
  }
}

void EBMultigridInterpolator::coarseFineInterpH(EBCellFAB& a_phi, const Interval a_variables, const DataIndex& a_dit){
  const Real m_dx     = 1.0;
  const Real m_dxCoar = m_dx*m_refRat;

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){

      const CFIVS* cfivsPtr = nullptr;
      if(sit() == Side::Lo) {
	cfivsPtr = &m_loCFIVS[dir][a_dit];
      }
      else{
	cfivsPtr = &m_hiCFIVS[dir][a_dit];
      }

      const IntVectSet& ivs = cfivsPtr->getFineIVS();
      if (cfivsPtr->isPacked() ){
	const int ihiorlo = sign(sit());
	FORT_INTERPHOMO(CHF_FRA(a_phi.getSingleValuedFAB()),
			CHF_BOX(cfivsPtr->packedBox()),
			CHF_CONST_REAL(m_dx),
			CHF_CONST_REAL(m_dxCoar),
			CHF_CONST_INT(dir),
			CHF_CONST_INT(ihiorlo));
      }
      else {
	if(!ivs.isEmpty()){
	  
	  Real halfdxcoar = m_dxCoar/2.0;
	  Real halfdxfine = m_dx/2.0;
	  Real xg = halfdxcoar -   halfdxfine;
	  Real xc = halfdxcoar +   halfdxfine;
	  Real xf = halfdxcoar + 3*halfdxfine;
	  Real hf = m_dx;
	  Real denom = xf*xc*hf;

	  const EBISBox& ebisBox = m_eblgFine.getEBISL()[a_dit];
	  const EBGraph& ebgraph = ebisBox.getEBGraph();
	      
	  for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	    const VolIndex& VoFGhost = vofit();

	    IntVect ivGhost  = VoFGhost.gridIndex();

	    Vector<VolIndex> farVoFs;
	    Vector<VolIndex> closeVoFs = ebisBox.getVoFs(VoFGhost, dir, flip(sit()), 1);
	    
	    bool hasClose = (closeVoFs.size() > 0);
	    bool hasFar = false;
	    Real phic = 0.0;
	    Real phif = 0.0;
	    
	    if (hasClose){
	      const int& numClose = closeVoFs.size();
	      for (int iVof=0;iVof<numClose;iVof++){
		const VolIndex& vofClose = closeVoFs[iVof];
		phic += a_phi(vofClose,0);
	      }
	      phic /= Real(numClose);

	      farVoFs = ebisBox.getVoFs(VoFGhost, dir, flip(sit()), 2);
	      hasFar   = (farVoFs.size()   > 0);
	      if (hasFar){
		const int& numFar = farVoFs.size();
		for (int iVof=0;iVof<numFar;iVof++){
		  const VolIndex& vofFar = farVoFs[iVof];
		  phif += a_phi(vofFar,0);
		}
		phif /= Real(numFar);
	      }
	    }

	    Real phiGhost;
	    if (hasClose && hasFar){
	      // quadratic interpolation  phi = ax^2 + bx + c
	      Real A = (phif*xc - phic*xf)/denom;
	      Real B = (phic*hf*xf - phif*xc*xc + phic*xf*xc)/denom;

	      phiGhost = A*xg*xg + B*xg;
	    }
	    else if (hasClose){
	      //linear interpolation
	      Real slope =  phic/xc;
	      phiGhost   =  slope*xg;
	    }
	    else{
	      phiGhost = 0.0; //nothing to interpolate from
	    }
	    a_phi(VoFGhost, 0) = phiGhost;
	  }
	}
      }
    }
  }
}

void EBMultigridInterpolator::defineGhosts(){

  const DisjointBoxLayout& dbl = m_eblgFine.getDBL();
  const ProblemDomain& domain  = m_eblgFine.getDomain();
  const EBISLayout& ebisl      = m_eblgFine.getEBISL();

  m_ghostCells.define(dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box cellBox      = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];
    Box grownBox = grow(cellBox, m_ghostCF);
    grownBox &= domain;

    m_ghostCells[dit()] = IntVectSet(grownBox);

    NeighborIterator nit(dbl);
    for (nit.begin(dit()); nit.ok(); ++nit){
	m_ghostCells[dit()] -= dbl[nit()];
    }

    m_ghostCells[dit()] -= cellBox;
  }
}

void EBMultigridInterpolator::defineGrids(const EBLevelGrid& a_eblgFine, const EBLevelGrid& a_eblgCoar){
  // Define the fine level grids. Use a shallow copy if we can.
  const int ghostReq  = 2*m_ghostCF;
  const int ghostFine = a_eblgFine.getGhost();
  if(ghostFine >= ghostReq){ 
    m_eblgFine = a_eblgFine;
  }
  else{ // Need to define from scratch
    m_eblgFine.define(m_eblgFine.getDBL(), m_eblgFine.getDomain(), ghostReq, m_eblgFine.getEBIS());
  }

  // Coarse is just a copy.
  m_eblgCoar = a_eblgCoar;

  // Define the coarsened fine grids -- needs same number of ghost cells as m_eblgFine. 
  coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
}

void EBMultigridInterpolator::defineData(){
  const DisjointBoxLayout& fineGrids = m_eblgFine.getDBL();
  const DisjointBoxLayout& coFiGrids = m_eblgCoFi.getDBL();

  const EBISLayout& fineEBISL        = m_eblgFine.getEBISL();
  const EBISLayout& coFiEBISL        = m_eblgCoFi.getEBISL();

  const ProblemDomain& fineDomain    = m_eblgFine.getDomain();
  const ProblemDomain& coarDomain    = m_eblgCoFi.getDomain();

  // Coarsen the fine-grid boxes and make a coarse layout, grown by a pretty big number of ghost cells...
  LayoutData<Box> grownFineBoxes(fineGrids);
  LayoutData<Box> grownCoarBoxes(coFiGrids);
  
  for (DataIterator dit(coFiGrids); dit.ok(); ++dit){
    Box coarBox = grow(coFiGrids[dit()], 2*m_ghostCF);
    Box fineBox = grow(fineGrids[dit()], 2*m_ghostCF);

    grownCoarBoxes[dit()] = coarBox & coarDomain;
    grownFineBoxes[dit()] = fineBox & fineDomain;
  }

  m_fineBoxes.define(grownFineBoxes);
  m_coarBoxes.define(grownCoarBoxes);

  m_fineData.define(m_fineBoxes, m_nComp, EBCellFactory(fineEBISL));
  m_coarData.define(m_coarBoxes, m_nComp, EBCellFactory(coFiEBISL));
}

void EBMultigridInterpolator::defineStencils(){
  const int comp  = 0;
  const int nComp = 1;

  const Real dxFine = 1.0;
  const Real dxCoar = dxFine*m_refRat;

  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();
  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();

  m_fineStencils.define(dblFine);
  m_coarStencils.define(dblFine);

  for (DataIterator dit(dblFine); dit.ok(); ++dit){
    const Box origFineBox  = dblFine[dit()];
    const Box grownFineBox = m_fineBoxes[dit()];
    const Box grownCoarBox = m_coarBoxes[dit()];

    // Define the valid regions such that the interpolation does not include coarse grid cells that fall beneath the fine level,
    // and no fine cells outside the CF.
    IntVectSet validFineCells = IntVectSet();
    IntVectSet validCoarCells = IntVectSet(grownCoarBox);

    // Coar mask is false everywhere under the fine grid box, and opposite for the fine mask. 
    for (BoxIterator bit(origFineBox); bit.ok(); ++bit){
      const IntVect fineIV = bit();
      const IntVect coarIV = coarsen(fineIV, m_refRat);

      validFineCells |= fineIV;
      validCoarCells -= coarIV;
    }

    // Same for parts of the current (grown) patch that overlaps with neighboring boxes.
    NeighborIterator nit(dblFine);
    for (nit.begin(dit()); nit.ok(); ++nit){
      const Box overlap  = grownFineBox & dblFine[nit()];
      
      for (BoxIterator bit(overlap); bit.ok(); ++bit){
	const IntVect fineIV = bit();
	validFineCells |= fineIV;
      }

      validCoarCells -= dblCoar[nit()];
    }

    // Now go through each ghost cell and get an interpolation stencil to specified order. 
    const EBISBox& ebisboxFine = m_eblgFine.getEBISL()[dit()];
    const EBISBox& ebisboxCoar = m_eblgCoFi.getEBISL()[dit()];

    const EBGraph& fineGraph   = ebisboxFine.getEBGraph();
    const EBGraph& coarGraph   = ebisboxCoar.getEBGraph();

    m_fineStencils[dit()].define(m_ghostCells[dit()], fineGraph, nComp);
    m_coarStencils[dit()].define(m_ghostCells[dit()], fineGraph, nComp);

    for (VoFIterator vofit(m_ghostCells[dit()], fineGraph); vofit.ok(); ++vofit){
      const VolIndex& ghostVof = vofit();

      VoFStencil& fineSten = m_fineStencils[dit()](ghostVof, comp);
      VoFStencil& coarSten = m_coarStencils[dit()](ghostVof, comp);

      int order         = m_order;
      bool foundStencil = false;
      
      while(order > 0 && !foundStencil){
	foundStencil = this->getStencil(fineSten,
					coarSten,
					m_cellLocation,					
					ghostVof,
					ebisboxFine,
					ebisboxCoar,
					validFineCells,
					validCoarCells,
					dxFine,
					dxCoar,
					order,
					m_weight);
	
	order--;
      }

      // Drop to order 0 if we never found a stencil, and issue an error code. 
      if(!foundStencil){
	fineSten.clear();

	int numCoarsen = 1;
	VolIndex ghostVofCoar = ghostVof;
	while(numCoarsen < m_refRat){
	  ghostVofCoar = ebisboxFine.coarsen(ghostVofCoar);
	  numCoarsen *= 2;
	}

	coarSten.add(ghostVofCoar, 1.0);
	
	MayDay::Warning("EBMultigridInterpolator::defineStencils -- could not find stencil and dropping to order 0");
      }
    }
  }
}

bool EBMultigridInterpolator::getStencil(VoFStencil&         a_stencilFine,
					 VoFStencil&         a_stencilCoar,
					 const CellLocation& a_cellLocation,
					 const VolIndex&     a_ghostVof,
					 const EBISBox&      a_ebisboxFine,
					 const EBISBox&      a_ebisboxCoar,
					 const IntVectSet&   a_validFineCells,
					 const IntVectSet&   a_validCoarCells,
					 const Real&         a_dxFine,
					 const Real&         a_dxCoar,
					 const int&          a_order,
					 const int&          a_weight){
  bool foundStencil = false;
  
  const int fineRadius = a_order;
  const int coarRadius = std::max(1, fineRadius/m_refRat);

  Vector<VolIndex> fineVofs;
  Vector<VolIndex> coarVofs;

  // Get the coarse vof which we obtain by coarsening a_ghostVof
  int numCoarsen = 1;
  VolIndex ghostVofCoar = a_ghostVof;
  while(numCoarsen < m_refRat){
    ghostVofCoar = a_ebisboxFine.coarsen(ghostVofCoar);
    numCoarsen *= 2;
  }

  // Get all Vofs in specified radii. Don't use cells that are not in a_validFineCells.
  fineVofs = VofUtils::getVofsInRadius(a_ghostVof,   a_ebisboxFine, fineRadius, VofUtils::Connectivity::MonotonePath, false);
  coarVofs = VofUtils::getVofsInRadius(ghostVofCoar, a_ebisboxCoar, coarRadius, VofUtils::Connectivity::MonotonePath, true );
  
  VofUtils::includeCells(fineVofs, a_validFineCells);
  VofUtils::includeCells(coarVofs, a_validCoarCells);
  
  VofUtils::onlyUnique(fineVofs);
  VofUtils::onlyUnique(coarVofs);

  const int numEquations = coarVofs.size() + fineVofs.size();
  const int numUnknowns  = LeastSquares::getTaylorExpansionSize(a_order);

  if(numEquations >= numUnknowns) { // We have enough equations to get a stencil

    // Build displacement vectors
    Vector<RealVect> fineDisplacements;
    Vector<RealVect> coarDisplacements;

    for (const auto& fineVof : fineVofs.stdVector()){
      fineDisplacements.push_back(LeastSquares::displacement(a_cellLocation, a_cellLocation, a_ghostVof, fineVof, a_ebisboxFine, a_dxFine));
    }

    for (const auto& coarVof : coarVofs.stdVector()){
      coarDisplacements.push_back(LeastSquares::displacement(a_cellLocation, a_cellLocation, a_ghostVof, coarVof, a_ebisboxFine, a_ebisboxCoar, a_dxFine, a_dxCoar));
    }


    // LeastSquares computes all unknown terms in a Taylor expansion up to specified order. We want the 0th order term, i.e. the interpolated value,
    // which in multi-index notation is the term (0,0), i.e. IntVect::Zero. The format of the two-level least squares routine is such that the
    // fine stencil lies on the first index. This can be confusing, but the LeastSquares uses a very compact notation...
    IntVect interpStenIndex = IntVect::Zero;
    IntVectSet derivs       = IntVectSet(interpStenIndex);
    IntVectSet knownTerms   = IntVectSet();

    std::map<IntVect, std::pair<VoFStencil, VoFStencil> > stencils = LeastSquares::computeDualLevelStencils(derivs,
													    knownTerms,
													    fineVofs,
													    coarVofs,
													    fineDisplacements,
													    coarDisplacements,
													    a_weight,
													    a_order);

    a_stencilFine = stencils.at(interpStenIndex).first;
    a_stencilCoar = stencils.at(interpStenIndex).second;

    if(a_stencilFine.size() != 0 || a_stencilCoar.size() != 0){
      foundStencil = true;
    }
    else{
      foundStencil = false;
    }
  }

  return foundStencil;
}

void EBMultigridInterpolator::interp(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phiCoar, const Interval a_variables) const {
  CH_TIME("EBMultigridInterpolator::interp");

  // Copy data to scratch data holders.
  a_phiFine.copyTo(m_fineData);
  a_phiCoar.copyTo(m_coarData);

  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit){
          EBCellFAB& dstFine = a_phiFine[dit()];
    const EBCellFAB& srcFine = m_fineData[dit()];
    const EBCellFAB& srcCoar = m_coarData[dit()];

    const BaseIVFAB<VoFStencil>& fineStencils = m_fineStencils[dit()];
    const BaseIVFAB<VoFStencil>& coarStencils = m_coarStencils[dit()];

    for (VoFIterator vofit(fineStencils.getIVS(), fineStencils.getEBGraph()); vofit.ok(); ++vofit){
      const VolIndex& ghostVoF = vofit();

      const VoFStencil& fineSten = fineStencils(ghostVoF, m_stenComp);
      const VoFStencil& coarSten = coarStencils(ghostVoF, m_stenComp);

      for (int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++){
	dstFine(ghostVoF, icomp) = 0.0;

	for (int ifine = 0; ifine < fineSten.size(); ifine++){
	  dstFine(ghostVoF, icomp) += fineSten.weight(ifine) * srcFine(fineSten.vof(ifine), icomp);
	}
	
	for (int icoar = 0; icoar < coarSten.size(); icoar++){
	  dstFine(ghostVoF, icomp) += coarSten.weight(icoar) * srcCoar(coarSten.vof(icoar), icomp);
	}
      }
    }
  }
}

void EBMultigridInterpolator::interpH(LevelData<EBCellFAB>& a_phiFine, const Interval a_variables) const {
  CH_TIME("EBMultigridInterpolator::interp");

  // Copy data to scratch data holders.
  a_phiFine.copyTo(m_fineData);

  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit){
          EBCellFAB& dstFine = a_phiFine[dit()];
    const EBCellFAB& srcFine = m_fineData[dit()];

    const BaseIVFAB<VoFStencil>& fineStencils = m_fineStencils[dit()];

    for (VoFIterator vofit(fineStencils.getIVS(), fineStencils.getEBGraph()); vofit.ok(); ++vofit){
      const VolIndex& ghostVoF = vofit();

      const VoFStencil& fineSten = fineStencils(ghostVoF, m_stenComp);

      for (int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++){
	dstFine(ghostVoF, icomp) = 0.0;

	for (int ifine = 0; ifine < fineSten.size(); ifine++){
	  dstFine(ghostVoF, icomp) += fineSten.weight(ifine) * srcFine(fineSten.vof(ifine), icomp);
	}
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
