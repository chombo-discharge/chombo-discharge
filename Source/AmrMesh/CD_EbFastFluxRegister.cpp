/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbFastFluxRegister.cpp
  @brief  Implementation of CD_EbFastFluxRegister.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <NeighborIterator.H>

// Our includes
#include <CD_Timer.H>
#include <CD_EbFastFluxRegister.H>
#include <CD_NamespaceHeader.H>

EbFastFluxRegister::EbFastFluxRegister(){
  CH_TIME("EbFastFluxRegister::EbFastFluxRegister");
}

EbFastFluxRegister::EbFastFluxRegister(const DisjointBoxLayout& a_dblFine,
				       const DisjointBoxLayout& a_dblCoar,
				       const EBISLayout&        a_ebislFine,
				       const EBISLayout&        a_ebislCoar,
				       const Box&               a_domainCoar,
				       const int&               a_nref,
				       const int&               a_nvar,
				       const EBIndexSpace*      a_ebisPtr,
				       const bool               a_forceNoEBCF){
  CH_TIME("EbFastFluxRegister::EbFastFluxRegister");
  
  this->define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar, a_domainCoar, a_nref, a_nvar, a_ebisPtr, a_forceNoEBCF);
}

EbFastFluxRegister::EbFastFluxRegister(const EBLevelGrid& a_eblgFine,
				       const EBLevelGrid& a_eblgCoar,
				       const int&         a_nref,
				       const int&         a_nvar,
				       const bool         a_forceNoEBCF){
  CH_TIME("EbFastFluxRegister::EbFastFluxRegister");
  
  this->define(a_eblgFine, a_eblgCoar, a_nref, a_nvar, a_forceNoEBCF);
}


EbFastFluxRegister::~EbFastFluxRegister(){
  CH_TIME("EbFastFluxRegister::~EbFastFluxRegister");  
}

void EbFastFluxRegister::define(const DisjointBoxLayout& a_gridsFine,
				const DisjointBoxLayout& a_gridsCoar,
				const EBISLayout&        a_ebislFine,
				const EBISLayout&        a_ebislCoar,
				const ProblemDomain&     a_domainCoar,
				const int&               a_nref,
				const int&               a_nvar,
				const EBIndexSpace*      a_ebisPtr,
				const bool               a_forceNoEBCF){
  CH_TIME("EbFastFluxRegister::define");
  
  ProblemDomain domainFine = refine(a_domainCoar, a_nref);
  EBLevelGrid eblgCoar(a_gridsCoar, a_ebislCoar, a_domainCoar);
  EBLevelGrid eblgFine(a_gridsFine, a_ebislFine,   domainFine);

  // Call other version. 
  this->define(eblgFine, eblgCoar, a_nref, a_nvar, a_forceNoEBCF);
}

void EbFastFluxRegister::define(const EBLevelGrid&       a_eblgFine,
				const EBLevelGrid&       a_eblgCoar,
				const int&               a_refRat,
				const int&               a_nvar,
				const bool               a_forceNoEBCF){
  CH_TIME("EbFastFluxRegister::define");

  Timer timer("EbFastFluxRegister::define");
  
  // Regular Chombo stuff
  this->clear();
  
  m_refRat   = a_refRat;
  m_nComp    = a_nvar;
  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;
  if (!m_eblgFine.coarsenable(a_refRat)){
    MayDay::Error("EbFastFluxRegister::define: dbl not coarsenable by refrat");
  }

  // Compute the coarsened fine grids. 
  timer.startEvent("Define coarsened grids");
  coarsen(m_eblgCoFi, a_eblgFine, a_refRat);
  m_eblgCoFi.getEBISL().setMaxRefinementRatio(a_refRat, m_eblgCoFi.getEBIS());
  m_reverseCopier.ghostDefine(m_eblgCoFi.getDBL(), m_eblgCoar.getDBL(),m_eblgCoar.getDomain(),IntVect::Unit);
  timer.stopEvent("Define coarsened grids");  

  m_nrefdmo = pow(m_refRat, SpaceDim-1);;

  // Define the regular flux register. 
  timer.startEvent("Define LevelFluxRegister");
  m_levelFluxReg = new LevelFluxRegister();
  m_levelFluxReg->define(a_eblgFine.getDBL(),
			 a_eblgCoar.getDBL(),
			 a_eblgFine.getDomain(),
			 a_refRat, a_nvar);
  timer.stopEvent("Define LevelFluxRegister");

  // Compute the existence of a coarse-fine boundary
  if (a_forceNoEBCF){
    m_hasEBCF = false;
  }
  else{
    timer.startEvent("Compute EBCF");
    m_hasEBCF = this->computeHasEBCF();
    timer.stopEvent("Compute EBCF");    
  }

  // If no EBCF then nothing happens. Otherwise we need to define
  // a bunch of temporaries that help us bookkeep fluxes over the
  // EBCF interface. 
  if(m_hasEBCF){
    timer.startEvent("Define masks");
    this->defineMasks();
    timer.stopEvent("Define masks");

    timer.startEvent("Define sets");
    this->fastDefineSetsAndIterators();
    timer.stopEvent("Define sets");

    timer.startEvent("Define buffers");
    this->defineBuffers();
    timer.stopEvent("Define buffers");    
  }
  
  //set all registers to zero to start.
  timer.startEvent("Set to zero");
  this->setToZero();
  timer.stopEvent("Set to zero");
  
  m_isDefined = true;

  ParmParse pp("EbFastFluxRegister");
  bool profile = false;
  pp.query("profile", profile);

  if(profile){
    timer.eventReport(pout(), false);
  }
}

void EbFastFluxRegister::fastDefineSetsAndIterators(){
  CH_TIME("EBFsatFluxRegister::fastDefineSetsAndIterators");

  // TLDR: The code for the CoFi things is essentially the same as that in Chombo. I've left that intact because
  //       it hasn't shown up as a bottleneck (yet). For the coarse stuff we've redone the define functions, avoiding
  //       the heinous call to getCFIVSSubset. 
  
  for (int idir = 0; idir < SpaceDim; idir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const int iindex = index(idir, sit());
      
      m_setsCoFi[iindex].define(m_eblgCoFi.getDBL());
      m_vofiCoFi[iindex].define(m_eblgCoFi.getDBL());

      for (DataIterator dit = m_eblgCoFi.getDBL().dataIterator(); dit.ok(); ++dit){
	const Box& grid = m_eblgCoFi.getDBL().get(dit());
	Box boxCoFi;
	if (sit() == Side::Lo){
	  boxCoFi = adjCellLo(grid, idir, 1);}
	else{
	  boxCoFi = adjCellHi(grid, idir, 1);
	}
	const EBGraph& ebgraphCoFi =    m_eblgCoFi.getEBISL()[dit()].getEBGraph();
	const IntVectSet& ivsCF    =  (*m_eblgCoFi.getCFIVS())[dit()];
	IntVectSet ivsEBCF         = ebgraphCoFi.getIrregCells(boxCoFi);

	ivsEBCF &= ivsCF;
	ivsEBCF &= boxCoFi;

	(m_setsCoFi[iindex])[dit()]      = ivsEBCF;
	(m_vofiCoFi[iindex])[dit()].define(ivsEBCF, ebgraphCoFi);
      }

      m_setsCoar[iindex].define(m_eblgCoar.getDBL());
      m_vofiCoar[iindex].define(m_eblgCoar.getDBL());

      // We make some changes here, and we identify the CFIVS using the mask. The coarse sets in EBFluxRegister are a bit silly, but they work. 
      for (DataIterator dit = m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit){
	const Box& boxCoar        = m_eblgCoar.getDBL()[dit()];
	const IntVectSet irregIVS = m_eblgCoar.getEBISL()[dit()].getIrregIVS(boxCoar);
	Vector<IntVectSet>  cfIVSVec;
	cfIVSVec.resize(0);

	// Add CFIVS identification from mask. 
	const IntVectSet& myIVS = m_coarCFIVS[iindex][dit()];
	if(!myIVS.isEmpty()) cfIVSVec.push_back(myIVS);

	const EBGraph& ebgraphCoar =  m_eblgCoar.getEBISL()[dit()].getEBGraph();
	(m_setsCoar[iindex])[dit()].resize(cfIVSVec.size());
	(m_vofiCoar[iindex])[dit()].resize(cfIVSVec.size());
	for (int ivec = 0; ivec < cfIVSVec.size(); ivec++){
	  IntVectSet ivsEBCF = irregIVS;
	  ivsEBCF &= cfIVSVec[ivec];

	  //need to grow by one to get on both sides of
	  //coarse-fine interface
	  (m_setsCoar[iindex])[dit()][ivec]      = ivsEBCF;
	  (m_vofiCoar[iindex])[dit()][ivec].define(ivsEBCF, ebgraphCoar);
	}
      }
    }
  }
}

void EbFastFluxRegister::defineMasks(){
  CH_TIME("EbFastFluxRegister::defineMasks");
  
  // TLDR: This code computes the coarse cells (i.e. mask) on the outside of the CFIVS on the coarse grid.
  //       We do this by using a mask on the CoFi grid. The boxes on that grid are grown by 1 cell in each direction
  //       and we make a BoxLayout<FArrayBox> that holds a value of 1 outside the CFIVS and 0 elsewhere. We set the
  //       mask on the coarse grid to 0, and then add the CoFi BoxLayout to the coarse grid. What we end up with is
  //       a LevelData<FArrayBox> which is 0 on all cells except the cells that are on the coarse side of the CFIVS. 

  constexpr int      iComp  = 0;
  constexpr int      nComp  = 1;
  constexpr Real     zero   = 0.0;
  constexpr Real     one    = 1.0;
  const     Interval interv = Interval(iComp, iComp);

  const DisjointBoxLayout& gridsCoar  = m_eblgCoar.getDBL();
  const DisjointBoxLayout& gridsCoFi  = m_eblgCoFi.getDBL();

  // Make a BoxLayoutData<FArrayBox> on gridsCoFi but where each box is grown by one cell.
  // We will be able to use DataIterator(gridsCoFi) to iterate through this construct. 
  LayoutData<Box> grownBoxes(gridsCoFi);
  for (DataIterator dit(gridsCoFi); dit.ok(); ++dit){
    grownBoxes[dit()] = grow(gridsCoFi[dit()], 1);
  }

  const BoxLayout grownCoFiGrid(grownBoxes);
  BoxLayoutData<FArrayBox> coFiMask(grownCoFiGrid, nComp); 
  
  // Go through and set ghost cells on the outside of the gridsCoFi to 1,
  // all valid cells and ghost cells "inside" the grid are set to zero. 
  for (int idir = 0; idir < SpaceDim; idir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const int iindex = index(idir, sit());

      // Define the coar mask in this direction.
      LevelData<FArrayBox> coarMask;
      LayoutData<IntVectSet>& coarCFIVS = m_coarCFIVS[iindex];
      coarMask.define(gridsCoar, nComp, IntVect::Zero);
      coarCFIVS.define(gridsCoar);

      // Coar mask is 0, we will copy from the coFiMask and into this one. 
      for (DataIterator dit(gridsCoar); dit.ok(); ++dit){
      	coarMask[dit()].setVal(zero);
      }
      
      // Set the coarsened fine mask to zero everywhere except immediately outside the coarsened fine grids. 
      for (DataIterator dit(gridsCoFi); dit.ok(); ++dit){
      	const Box cellBoxCoFi  = gridsCoFi    [dit()];
	const Box grownBoxCoFi = grownCoFiGrid[dit()];

      	// Get CF cells outside the CoFi box. Subtract the ungrown (i.e., valid region) neighbor boxes because they might overlap with
	// the outside of this box. 
      	const Box coarseFineInterfaceBox = adjCellBox(cellBoxCoFi, idir, sit(), 1);
      	DenseIntVectSet cfIvsCoFi(coarseFineInterfaceBox, true);

	NeighborIterator nit(gridsCoFi);	
      	for (nit.begin(dit()); nit.ok(); ++nit){
      	  const Box neighborBox = gridsCoFi[nit()];
	  const Box overlapBox  = cellBoxCoFi & neighborBox;

	  if(!overlapBox.isEmpty()){
	    cfIvsCoFi -= neighborBox;
	  }
      	}

      	// Set the coarase-fine mask. This also sets the ghost cells.
	coFiMask[dit()].setVal(zero);
	for (DenseIntVectSetIterator ivsIt(cfIvsCoFi); ivsIt.ok(); ++ivsIt){
	  coFiMask[dit()](ivsIt(), iComp) = one;
	}
      }

      // Add the coarsened fine mask to the coarse mask. This means that we add the contents in the "ghost regions" on the coarsened fine mask
      // onto the valid region in the coarse mask, leaving us with a view of which valid cells in the coarse grid abut the fine-grid boundary. 
      coFiMask.addTo(Interval(0,0), coarMask, Interval(0,0), m_eblgCoar.getDomain());

      // Go through the mask and set the CFIVS accordingly. 
      for (DataIterator dit = gridsCoar.dataIterator(); dit.ok(); ++dit){
      	const Box coarBox = gridsCoar[dit()];
      	IntVectSet& cfivs = coarCFIVS[dit()];

      	cfivs.makeEmpty();
      	for (BoxIterator bit(coarBox); bit.ok(); ++bit){
      	  if(coarMask[dit()](bit(), iComp) > zero) {
	    cfivs |= bit();
	  }
      	}
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
