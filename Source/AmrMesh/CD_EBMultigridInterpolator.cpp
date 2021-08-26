/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBMultigridInterpolator.cpp
  @brief  Implementation of CD_EBMultigridInterpolator.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBCellFactory.H>
#include <EBLevelDataOps.H>
#include <NeighborIterator.H>
#include <EBAlias.H>
#include <ParmParse.H>

// Our includes
#include <CD_Timer.H>
#include <CD_EBMultigridInterpolator.H>
#include <CD_EBMultigridInterpolatorF_F.H>
#include <CD_VofUtils.H>
#include <CD_LeastSquares.H>
#include <CD_NamespaceHeader.H>

constexpr int EBMultigridInterpolator::m_stenComp;
constexpr int EBMultigridInterpolator::m_nStenComp;
constexpr int EBMultigridInterpolator::m_comp;

EBMultigridInterpolator::EBMultigridInterpolator(const EBLevelGrid& a_eblgFine,
						 const EBLevelGrid& a_eblgCoar,
						 const CellLocation a_cellLocation,
						 const IntVect&     a_minGhost,
						 const int          a_refRat,
						 const int          a_nVar,
						 const int          a_ghostCF,
						 const int          a_order,
						 const int          a_weighting){
  CH_assert(a_ghostCF   > 0);
  CH_assert(a_nVar      > 0);
  CH_assert(a_refRat%2 == 0);

  const DisjointBoxLayout& gridsFine = a_eblgFine.getDBL();
  const DisjointBoxLayout& gridsCoar = a_eblgCoar.getDBL();

  Timer timer("EBMultigridInterpolator::EBMultigridInterpolator");

  timer.startEvent("Define QuadCFInterp");
  QuadCFInterp::define(gridsFine, &gridsCoar, 1.0, a_refRat, a_nVar, a_eblgFine.getDomain());
  timer.stopEvent("Define QuadCFInterp");  
  
  m_refRat        = a_refRat;
  m_nComp         = a_nVar;
  m_ghostCF       = a_ghostCF;
  m_order         = a_order;
  m_minGhost      = a_minGhost;
  m_cellLocation  = a_cellLocation;
  m_weight        = a_weighting;

  timer.startEvent("Define grids");  
  this->defineGrids(a_eblgFine, a_eblgCoar);
  timer.stopEvent("Define grids");    

  timer.startEvent("Define ghost regions");    
  this->defineGhostRegions();
  timer.stopEvent("Define ghost regions");      

  timer.startEvent("Define buffers");      
  this->defineBuffers();
  timer.stopEvent("Define buffers");        

  timer.startEvent("Define stencils");        
  this->defineStencils();
  timer.stopEvent("Define stencils");          

  timer.startEvent("Set coarsening ratio");          
  if(m_eblgFine.getMaxCoarseningRatio() < m_refRat){
    m_eblgFine.setMaxCoarseningRatio(m_refRat, m_eblgFine.getEBIS());
  }
  timer.stopEvent("Set coarsening ratio");            

  ParmParse pp("EBMultigridInterpolator");
  bool profile = false;
  pp.query("profile", profile);
  if(profile){
    timer.eventReport();
  }

  m_isDefined = true;
}

int EBMultigridInterpolator::getGhostCF() const{
  return m_ghostCF;
}

EBMultigridInterpolator::~EBMultigridInterpolator(){

}

void EBMultigridInterpolator::coarseFineInterp(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phiCoar, const Interval a_variables){
  CH_assert(m_ghostCF <= a_phiFine.ghostVect().max());
  CH_assert(m_ghostCF <= a_phiCoar.ghostVect().max());

  CH_TIME("EBMultigridInterpolator::interp");

  LevelData<FArrayBox> fineAlias;
  LevelData<FArrayBox> coarAlias;

  aliasEB(fineAlias, (LevelData<EBCellFAB>&) a_phiFine);
  aliasEB(coarAlias, (LevelData<EBCellFAB>&) a_phiCoar);

  QuadCFInterp::coarseFineInterp(fineAlias, coarAlias);

  // Do corrections near the EBCF.
  for (int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++){

    // Copy data to scratch data holders.
    const Interval srcInterv = Interval(icomp,   icomp);
    const Interval dstInterv = Interval(m_comp,  m_comp);
    
    // Need coarse data to be accesible by fine data, so I always have to copy here. 
    a_phiCoar.copyTo(srcInterv, m_grownCoarData, dstInterv);

    for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit){
      EBCellFAB& dstFine       = a_phiFine[dit()];

      const EBCellFAB& srcFine = a_phiFine[dit()];
      const EBCellFAB& srcCoar = m_grownCoarData[dit()];

      const BaseIVFAB<VoFStencil>& fineStencils = m_fineStencils[dit()];
      const BaseIVFAB<VoFStencil>& coarStencils = m_coarStencils[dit()];

      for (VoFIterator vofit(fineStencils.getIVS(), fineStencils.getEBGraph()); vofit.ok(); ++vofit){
	const VolIndex& ghostVoF = vofit();

	dstFine(ghostVoF, icomp) = 0.0;

	// Fine and coarse stencils for this ghost vof
	const VoFStencil& fineSten = fineStencils(ghostVoF, m_stenComp);
	const VoFStencil& coarSten = coarStencils(ghostVoF, m_stenComp);

	// Apply fine stencil
	for (int ifine = 0; ifine < fineSten.size(); ifine++){
	  dstFine(ghostVoF, icomp) += fineSten.weight(ifine) * srcFine(fineSten.vof(ifine), icomp);
	}

	// Apply coarse stencil
	for (int icoar = 0; icoar < coarSten.size(); icoar++){
	  dstFine(ghostVoF, icomp) += coarSten.weight(icoar) * srcCoar(coarSten.vof(icoar), m_comp);
	}
      }
    }
  }
}

void EBMultigridInterpolator::coarseFineInterpH(LevelData<EBCellFAB>& a_phiFine, const Interval a_variables){
  CH_assert(m_ghostCF <= a_phiFine.ghostVect().max());
  for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit){
    this->coarseFineInterpH(a_phiFine[dit()], a_variables, dit());
  }
}

void EBMultigridInterpolator::coarseFineInterpH(EBCellFAB& a_phi, const Interval a_variables, const DataIndex& a_dit){
  const Real dxFine = 1.0;
  const Real dxCoar = 1.0 * m_refRat;
  
  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++){
    
    // Do the regular interp on all sides of the current patch. 
    for (int dir = 0; dir < SpaceDim; dir++){
      for (SideIterator sit; sit.ok(); ++sit){

	// Regular homogeneous interpolation on side sit() in direction dir
	const Box ghostBox = m_cfivs[a_dit].at(std::make_pair(dir, sit()));
	
	if(!ghostBox.isEmpty()){
	  const int hiLo = sign(sit()); // Low side => -1, high side => 1
	  BaseFab<Real>& phiReg = a_phi.getSingleValuedFAB();
	  FORT_MGINTERPHOMO(CHF_FRA1(phiReg, ivar),
	  		    CHF_CONST_REAL(dxFine),
	  		    CHF_CONST_REAL(dxCoar),
	  		    CHF_CONST_INT(dir),
	  		    CHF_CONST_INT(hiLo),
	  		    CHF_BOX(ghostBox));
	}
      }
    }

    // Apply fine stencil near the EB. This might look weird but recall that ghostVof is a ghostCell on the other side of
    // the refinment boundary, whereas the stencil only reaches into cells on the finel level. So, we are not
    // writing to data that we will later fetch.
    const BaseIVFAB<VoFStencil>& fineStencils = m_fineStencils[a_dit];
    
    for (VoFIterator vofit(fineStencils.getIVS(), fineStencils.getEBGraph()); vofit.ok(); ++vofit){
      const VolIndex& ghostVoF   = vofit();
      const VoFStencil& fineSten = fineStencils(ghostVoF, m_stenComp);
      
      a_phi(ghostVoF, ivar) = 0.0;

      for (int ifine = 0; ifine < fineSten.size(); ifine++){
	CH_assert(ghostVoF != fineSten.vof(ifine));
	a_phi(ghostVoF, ivar) += fineSten.weight(ifine) * a_phi(fineSten.vof(ifine), ivar);
      }
    }
  }
}

void EBMultigridInterpolator::defineGhostRegions(){
  const DisjointBoxLayout& dbl = m_eblgFine.getDBL();
  const ProblemDomain& domain  = m_eblgFine.getDomain();
  const EBISLayout& ebisl      = m_eblgFine.getEBISL();

  // Define the "regular" ghost interpolation regions. This is just one cell wide since the operator stencil
  // has a width of 1 in regular cells. 
  m_cfivs.define(dbl);
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box cellBox = dbl[dit()];

    std::map<std::pair<int, Side::LoHiSide>, Box>& cfivsBoxes = m_cfivs[dit()];
    
    for (int dir = 0; dir < SpaceDim; dir++){
      for (SideIterator sit; sit.ok(); ++sit){

	IntVectSet cfivs = IntVectSet(adjCellBox(cellBox, dir, sit(), 1));

	NeighborIterator nit(dbl); // Subtract the other boxes if they intersect this box. 
	for (nit.begin(dit()); nit.ok(); ++nit){
	  cfivs -= dbl[nit()];
	}

	cfivs &= domain;
	cfivs.recalcMinBox();
	
	cfivsBoxes.emplace(std::make_pair(dir, sit()), cfivs.minBox());
      }
    }
  }

  // Define ghost cells to be interpolated near the EB. 
  m_ghostCells.define(dbl);
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box cellBox      = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];

    if(ebisbox.isAllRegular() || ebisbox.isAllCovered()){
      m_ghostCells[dit()] = IntVectSet();
    }
    else{

      // 1. Define the width of the ghost layer region around current (irregular grid) patch
      Box grownBox = grow(cellBox, m_ghostCF);
      grownBox &= domain;

      m_ghostCells[dit()] = IntVectSet(grownBox);

      NeighborIterator nit(dbl);
      for (nit.begin(dit()); nit.ok(); ++nit){
	m_ghostCells[dit()] -= dbl[nit()];
      }
      m_ghostCells[dit()] -= cellBox;

      // 2. Only include ghost cells that are within range m_ghostCF of an irregular grid cell
      IntVectSet irreg = ebisbox.getIrregIVS(cellBox);
      irreg.grow(m_ghostCF);
      m_ghostCells[dit()] &= irreg;
    }
  }
}

void EBMultigridInterpolator::defineGrids(const EBLevelGrid& a_eblgFine, const EBLevelGrid& a_eblgCoar){
  // Define the fine level grids. Use a shallow copy if we can.
  const int ghostReq  = m_ghostCF + m_order;
  const int ghostFine = a_eblgFine.getGhost();

  if(ghostFine >= ghostReq){ 
    m_eblgFine = a_eblgFine;
  }
  else{ // Need to define from scratch
    m_eblgFine = EBLevelGrid(a_eblgFine.getDBL(), a_eblgFine.getDomain(), ghostReq, a_eblgFine.getEBIS());
  }

  // Coarse is just a copy.
  
  m_eblgCoar = a_eblgCoar;

  // Define the coarsened fine grids -- needs same number of ghost cells as m_eblgFine.
  coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
  m_eblgCoFi.setMaxRefinementRatio(m_refRat);
}

void EBMultigridInterpolator::defineBuffers(){

  // TLDR: On both the fine and the coarse level we want to have a BoxLayout which is just like the DisjointBoxLayout except that
  //       the boxes should be grown on each level. 
  Timer timer("EBMultigridInterpolator::defineBuffers");

  const int ghostReq = m_ghostCF + m_order;
  
  const DisjointBoxLayout& fineGrids = m_eblgFine.getDBL();
  const DisjointBoxLayout& coFiGrids = m_eblgCoFi.getDBL();

  const EBISLayout& fineEBISL        = m_eblgFine.getEBISL();
  const EBISLayout& coFiEBISL        = m_eblgCoFi.getEBISL();

  const ProblemDomain& fineDomain    = m_eblgFine.getDomain();
  const ProblemDomain& coarDomain    = m_eblgCoFi.getDomain();

  // Coarsen the fine-grid boxes and make a coarse layout, grown by a pretty big number of ghost cells...
  timer.startEvent("make lD");
  LayoutData<Box> grownFineBoxesLayout(fineGrids);
  LayoutData<Box> grownCoarBoxesLayout(coFiGrids);
  timer.stopEvent("make lD");  

  timer.startEvent("dit loop");
  for (DataIterator dit(coFiGrids); dit.ok(); ++dit){
    Box coarBox = grow(coFiGrids[dit()], ghostReq);
    Box fineBox = grow(fineGrids[dit()], ghostReq);

    grownCoarBoxesLayout[dit()] = coarBox & coarDomain;
    grownFineBoxesLayout[dit()] = fineBox & fineDomain;
  }
  timer.stopEvent("dit loop");  

  timer.startEvent("defines 1");
  m_grownFineBoxesLayout.define(grownFineBoxesLayout);
  m_grownCoarBoxesLayout.define(grownCoarBoxesLayout);
  timer.stopEvent("defines 1");  

  timer.startEvent("defines 2");
  m_grownFineData.define(m_grownFineBoxesLayout, m_nComp, EBCellFactory(fineEBISL));
  m_grownCoarData.define(m_grownCoarBoxesLayout, m_nComp, EBCellFactory(coFiEBISL));
  timer.stopEvent("defines 2");

  timer.eventReport();
}

void EBMultigridInterpolator::defineStencils(){
  const int comp  = 0;
  const int nComp = 1;

  const Real dxFine = 1.0;
  const Real dxCoar = dxFine*m_refRat;

  const ProblemDomain& domFine = m_eblgFine.getDomain();
  const ProblemDomain& domCoar = m_eblgCoFi.getDomain();

  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();
  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();

  const EBISLayout& ebislFine = m_eblgFine.getEBISL();
  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();

  m_fineStencils.define(dblFine);
  m_coarStencils.define(dblFine);

  for (DataIterator dit(dblFine); dit.ok(); ++dit){
    const Box origFineBox    = dblFine[dit()];
    const Box ghostedFineBox = grow(origFineBox, m_minGhost);
    const Box grownFineBox   = m_grownFineBoxesLayout[dit()];
    const Box grownCoarBox   = m_grownCoarBoxesLayout[dit()];

    // Define the valid regions such that the interpolation does not include coarse grid cells that fall beneath the fine level,
    // and no fine cells outside the CF.
    DenseIntVectSet validFineCells(origFineBox,  true); // Original patch. We will this patch's ghost cells if those ghost cells overlap with other patches on this level. 
    DenseIntVectSet validCoarCells(grownCoarBox, true); // 

    // Subtract the coarse cells under the fine box. 
    Box cBox = origFineBox;
    cBox.coarsen(m_refRat);
    validCoarCells -= cBox;

    // Same for parts of the current (grown) patch that overlaps with neighboring boxes.
    NeighborIterator nit(dblFine);
    for (nit.begin(dit()); nit.ok(); ++nit){
      const Box fineOverlap  = grownFineBox & dblFine[nit()];
      
      Box coarOverlap = fineOverlap;
      coarOverlap.coarsen(m_refRat);
      
      validFineCells |= DenseIntVectSet(fineOverlap, true);
      validCoarCells -= DenseIntVectSet(coarOverlap, true);
    }

    // Restrict to the ghosted input box
    validFineCells &= ghostedFineBox;

    // Restrict to domain
    validFineCells &= domFine;
    validCoarCells &= domCoar;
    
    // Now go through each ghost cell and get an interpolation stencil to specified order.
    const EBISBox& ebisboxFine = m_eblgFine.getEBISL()[dit()];
    const EBISBox& ebisboxCoar = m_eblgCoFi.getEBISL()[dit()];

    const EBGraph& fineGraph   = ebisboxFine.getEBGraph();
    const EBGraph& coarGraph   = ebisboxCoar.getEBGraph();

    m_fineStencils[dit()].define(m_ghostCells[dit()], fineGraph, nComp);
    m_coarStencils[dit()].define(m_ghostCells[dit()], fineGraph, nComp);

    for (VoFIterator vofit(m_ghostCells[dit()], fineGraph); vofit.ok(); ++vofit){

      const VolIndex& ghostVofFine = vofit();
      const VolIndex& ghostVofCoar = ebislFine.coarsen(ghostVofFine, m_refRat, dit());
      
      VoFStencil& fineSten = m_fineStencils[dit()](ghostVofFine, comp);
      VoFStencil& coarSten = m_coarStencils[dit()](ghostVofFine, comp);

      int order         = m_order;
      bool foundStencil = false;
						       
      while(order > 0 && !foundStencil){
	foundStencil = this->getStencil(fineSten,
					coarSten,
					m_cellLocation,
					ghostVofFine,
					ghostVofCoar,
					ebisboxFine,
					ebisboxCoar,
					validFineCells,
					validCoarCells,
					dxFine,
					dxCoar,
					order,
					m_weight);
	
	order--;

	if(!foundStencil){
	  pout() << "EBMultigridInterpolator -- on domain = " << m_eblgFine.getDomain() << ", dropping order for vof = " << ghostVofFine << endl;
	}
      }

      // Drop to order 0 if we never found a stencil, and issue an error code. 
      if(!foundStencil){
	fineSten.clear();
	coarSten.clear();
	
	coarSten.add(ghostVofCoar, 1.0);
	
	pout() << "EBMultigridInterpolator::defineStencils -- could not find stencil and dropping to order 0" << endl;
      }
    }
  }
}

bool EBMultigridInterpolator::getStencil(VoFStencil&            a_stencilFine,
					 VoFStencil&            a_stencilCoar,
					 const CellLocation&    a_cellLocation,
					 const VolIndex&        a_ghostVofFine,
					 const VolIndex&        a_ghostVofCoar,
					 const EBISBox&         a_ebisboxFine,
					 const EBISBox&         a_ebisboxCoar,
					 const DenseIntVectSet& a_validFineCells,
					 const DenseIntVectSet& a_validCoarCells,
					 const Real&            a_dxFine,
					 const Real&            a_dxCoar,
					 const int&             a_order,
					 const int&             a_weight){
  bool foundStencil = true;
  
  const int fineRadius = a_order;
  const int coarRadius = std::max(2, fineRadius/m_refRat);

  Vector<VolIndex> fineVofs;
  Vector<VolIndex> coarVofs;


  // Get all Vofs in specified radii. Don't use cells that are not in a_validFineCells.
  fineVofs = VofUtils::getVofsInRadius(a_ghostVofFine, a_ebisboxFine, fineRadius, VofUtils::Connectivity::MonotonePath, false);
  coarVofs = VofUtils::getVofsInRadius(a_ghostVofCoar, a_ebisboxCoar, coarRadius, VofUtils::Connectivity::MonotonePath, true );
  
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
      fineDisplacements.push_back(LeastSquares::displacement(a_cellLocation, a_cellLocation, a_ghostVofFine, fineVof, a_ebisboxFine, a_dxFine));
    }

    for (const auto& coarVof : coarVofs.stdVector()){
      coarDisplacements.push_back(LeastSquares::displacement(a_cellLocation, a_cellLocation, a_ghostVofFine, coarVof, a_ebisboxFine, a_ebisboxCoar, a_dxFine, a_dxCoar));
    }

    // LeastSquares computes all unknown terms in a Taylor expansion up to specified order. We want the 0th order term, i.e. the interpolated value,
    // which in multi-index notation is the term (0,0), i.e. IntVect::Zero. The format of the two-level least squares routine is such that the
    // fine stencil lies on the first index. This can be confusing, but the LeastSquares uses a very compact notation. 
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

    foundStencil = true;
  }
  else{
    foundStencil = false;
  }

  return foundStencil;
}

#include <CD_NamespaceFooter.H>
