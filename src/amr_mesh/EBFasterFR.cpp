/*!
  @brief   EBFasterFR.cpp
  @details Implementation of EBFasterFR.H
  @author  Robert Maskar
  @date    Aug. 2020
*/

#include "EBFasterFR.H"
#include "NeighborIterator.H"

EBFasterFR::EBFasterFR(){

}

EBFasterFR::EBFasterFR(const EBLevelGrid&       a_eblgFine,
		       const EBLevelGrid&       a_eblgCoar,
		       const int&               a_nref,
		       const int&               a_nvar,
		       const bool               a_forceNoEBCF){
  setDefaultValues();
  define(a_eblgFine, a_eblgCoar, a_nref, a_nvar, a_forceNoEBCF);
}

EBFasterFR::~EBFasterFR(){

}

void EBFasterFR::define(const EBLevelGrid&       a_eblgFine,
			const EBLevelGrid&       a_eblgCoar,
			const int&               a_refRat,
			const int&               a_nvar,
			const bool               a_forceNoEBCF){
  CH_TIME("EBFasterFR::define");
  
  clear();
  m_refRat   = a_refRat;
  m_nComp    = a_nvar;
  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;
  if (!m_eblgFine.coarsenable(a_refRat))
    {
      MayDay::Error("EBFastFR::define: dbl not coarsenable by refrat");
    }

  coarsen(m_eblgCoFi, a_eblgFine, a_refRat);
  m_eblgCoFi.getEBISL().setMaxRefinementRatio(a_refRat, m_eblgCoFi.getEBIS());
  m_reverseCopier.ghostDefine(m_eblgCoFi.getDBL(), m_eblgCoar.getDBL(),m_eblgCoar.getDomain(),IntVect::Unit);

#if (CH_SPACEDIM == 2)
  m_nrefdmo = m_refRat;
#elif (CH_SPACEDIM == 3)
  m_nrefdmo = m_refRat*m_refRat;
#else
  bogus spacedim;
#endif

  m_levelFluxReg = new LevelFluxRegister();

  //  pout() << "before regular define" << endl;
  m_levelFluxReg->define(a_eblgFine.getDBL(),
                         a_eblgCoar.getDBL(),
                         a_eblgFine.getDomain(),
                         a_refRat, a_nvar);
  if (a_forceNoEBCF)
    {
      m_hasEBCF = false;
    }
  else
    {
      m_hasEBCF = computeHasEBCF();
    }

  //if no EBCF--nothing happens here but calls to level flux register
  if (m_hasEBCF)
    {
      defineMasks();
      fastDefineSetsAndIterators();
      defineBuffers();
    }
  //set all registers to zero to start.
  setToZero();
  m_isDefined = true;
  //  pout() << "leaving define" << endl;
}

void EBFasterFR::fastDefineSetsAndIterators(){
  CH_TIME("EBFasterFR::fastDefineSetsAndIterators");
  
  for (int idir = 0; idir < SpaceDim; idir++){
    for (SideIterator sit; sit.ok(); ++sit){
      int iindex = index(idir, sit());
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

      for (DataIterator dit = m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit){
	const Box& boxCoar = m_eblgCoar.getDBL()[dit()];
	IntVectSet irregIVS = m_eblgCoar.getEBISL()[dit()].getIrregIVS(boxCoar);
	Vector<IntVectSet>  cfIVSVec;
	int iindex = index(idir, sit());
#if 0 // original code
	getCFIVSSubset(cfIVSVec, boxCoar, m_eblgCoFi.getDBL(), idir, sit());
#else
	IntVectSet myIVS = IntVectSet();
	FArrayBox& coarMask = m_coarMask[iindex][dit()];
	for (BoxIterator bit(boxCoar); bit.ok(); ++bit){
	  if(coarMask(bit()) > 0.0){
	    myIVS |= bit();
	  }
	}
	cfIVSVec.resize(0);
	if(!myIVS.isEmpty()){
	  cfIVSVec.push_back(myIVS);
	}
#endif
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

void EBFasterFR::defineMasks(){
  CH_TIME("EBFasterFR::defineMasks");

  // TLDR: This code computes the coarse cells (i.e. mask) on the outside of the CFIVS on the coarse grid.
  //       We do this by using a mask on the CoFi grid. The boxes on that grid are grown by 1 cell in each direction
  //       and we make a BoxLayout<FArrayBox> that holds a value of 1 outside the CFIVS and 0 elsewhere. We set the
  //       mask on the coarse grid to 0, and then add the CoFi BoxLayout to the coarse grid. What we end up with is
  //       a LevelData<FArrayBox> which is 0 on all cells except the cells that are on the coarse side of the CFIVS. 

  const int icomp = 0;
  const int ncomp = 1;

  const DisjointBoxLayout& gridsCoar = m_eblgCoar.getDBL();
  const DisjointBoxLayout& gridsCoFi = m_eblgCoFi.getDBL();

  const ProblemDomain& coarDomain = m_eblgCoar.getDomain();

  // Make a BoxLayoutData<FArrayBox> on gridsCoFi but where each box is grown by one cell.
  Vector<Box> coFiBoxes = gridsCoFi.boxArray();
  Vector<int> coFiProcs = gridsCoFi.procIDs();

  const int ghostCoFi = 1;
  for (int i = 0; i < coFiBoxes.size(); i++){
    coFiBoxes[i].grow(ghostCoFi);
    coFiBoxes[i] &= coarDomain;
  }

  const DisjointBoxLayout grownCoFiGrid(coFiBoxes, coFiProcs);
  BoxLayoutData<FArrayBox> coFiMask(grownCoFiGrid, ncomp);

  // Go through and set ghost cells on the outside of the gridsCoFi to 1,
  // all valid cells and ghost cells "inside" the grid are set to zero
  for (int idir = 0; idir < SpaceDim; idir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const int iindex = index(idir, sit());

      // Define the coar mask in this direction.
      LevelData<FArrayBox>& coarMask = m_coarMask[iindex];
      coarMask.define(gridsCoar, ncomp, IntVect::Zero);

      // Coar mask is 0, we will copy from the coFiMask and into this one. 
      for (DataIterator dit = gridsCoar.dataIterator(); dit.ok(); ++dit){
	coarMask[dit()].setVal(0.0);
      }

      // Fine mask is also 0, except on the CF region. 
      NeighborIterator nit(gridsCoFi);
      for (DataIterator dit = gridsCoFi.dataIterator(); dit.ok(); ++dit){
	const Box boxCoFi = gridsCoFi[dit()];

	// Get CF cells outside the CoFi box. Subtract the neighbor boxes. 
	Box cfBoxCoFi = adjCellBox(boxCoFi, idir, sit(), 1);
	IntVectSet cfIvsCoFi(cfBoxCoFi);
	for (nit.begin(dit()); nit.ok(); ++nit){
	  const Box neighborBox = gridsCoFi[nit()];

	  cfIvsCoFi -= neighborBox;
	}

	// Set the CoFi mask. This also sets the ghost cells.
	coFiMask[dit()].setVal(0.0);
	for (IVSIterator ivsIt(cfIvsCoFi); ivsIt.ok(); ++ivsIt){
	  coFiMask[dit()](ivsIt(), icomp) = 1.0;
	}
      }

      // Copy CoFi mask to coarse level. But this only copies valid cells, *sigh*
      coFiMask.addTo(Interval(0,0), coarMask, Interval(0,0), m_eblgCoar.getDomain());

#if 0 // Dummy check
      for (DataIterator dit = gridsCoar.dataIterator(); dit.ok(); ++dit){
	const FArrayBox& coarBoxMask = coarMask[dit()];
	const Box bx = gridsCoar[dit()];

	for (BoxIterator bit(bx); bit.ok(); ++bit){
	  const IntVect iv = bit();
	  if(coarBoxMask(iv, 0) > 0.0){
	    std::cout << "hoozah" << std::endl;
	  }
	}
      }
#endif
    }
  }
}

void EBFasterFR::getCFIVSSubset(Vector<IntVectSet>&      a_cfivsCoFi,
				const Box&               a_subreCoar,
				const DisjointBoxLayout& a_layouCoFi,
				const int&               a_idir,
				const Side::LoHiSide&    a_sd){

  a_cfivsCoFi.resize(0);
  for (LayoutIterator lito = a_layouCoFi.layoutIterator(); lito.ok(); ++lito){
    const Box& gridCoFiO = a_layouCoFi[lito()];
    //the - 1 is due to the adjacent cell box thing
    if (gridCoFiO.bigEnd(0) < a_subreCoar.smallEnd(0) - 1 ){
      //can skip rest cuz we haven't gotten to something interesting
      continue;
    }

    Box cfboxCoFi = adjCellBox(gridCoFiO, a_idir, a_sd, 1);

    //only want the cells within the input subregion
    cfboxCoFi &= a_subreCoar;

    IntVectSet ivsCoFi(cfboxCoFi);
    for (LayoutIterator liti = a_layouCoFi.layoutIterator(); liti.ok(); ++liti){
      const Box& gridCoFiI = a_layouCoFi[liti()];
      if (gridCoFiI.bigEnd(0) < cfboxCoFi.smallEnd(0)){
	//can skip rest cuz we haven't gotten to something interesting
	continue;
      }

      //remove from the set boxes that overlap on the own grid
      //(fine-fine interface points)
      ivsCoFi -= gridCoFiI;

      if ((gridCoFiI.smallEnd(0) > cfboxCoFi.bigEnd(0)) || ivsCoFi.isEmpty()){
	//can break out of loop, since we know that the smallEnd
	// of all the remaining boxes are lexigraphically beyond this ghosted box.
	break;
      }
    }
    if (!ivsCoFi.isEmpty()){
      a_cfivsCoFi.push_back(ivsCoFi);
    }

    //the + 1 is due to the adjacent cell box thing
    if (gridCoFiO.smallEnd(0) > a_subreCoar.bigEnd(0) + 1 ){
      //can break out of loop, since we know that the smallEnd
      // of all the remaining boxes are lexigraphically beyond this ghosted box.
      break;
    }
  }
}
