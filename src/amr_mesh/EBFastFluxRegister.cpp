/*!
  @file   EBFastFluxRegister.cpp
  @brief  Implementation of EBFastFluxRegister.H
  @author Robert Marskar
  @date   Aug. 2020
*/

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#include <iomanip>
#include "EBArith.H"
#include "EBFastFluxRegister.H"
#include "LayoutIterator.H"
#include "BaseIVFactory.H"
#include "EBCoarToCoarRedist.H"
#include "EBCoarToFineRedist.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "EBIndexSpace.H"
#include "parstream.H"
#include "EBLevelDataOps.H"

#include "CH_Timer.H"
#include "NeighborIterator.H"

#include "CD_NamespaceHeader.H"

EBFastFluxRegister::EBFastFluxRegister(){

}

EBFastFluxRegister::EBFastFluxRegister(const DisjointBoxLayout& a_dblFine,
				       const DisjointBoxLayout& a_dblCoar,
				       const EBISLayout&        a_ebislFine,
				       const EBISLayout&        a_ebislCoar,
				       const Box&               a_domainCoar,
				       const int&               a_nref,
				       const int&               a_nvar,
				       const EBIndexSpace*      ebisPtr,
				       const bool               a_forceNoEBCF){
  define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
	 a_domainCoar, a_nref, a_nvar, ebisPtr, a_forceNoEBCF);
}

EBFastFluxRegister::EBFastFluxRegister(const EBLevelGrid& a_eblgFine,
				       const EBLevelGrid& a_eblgCoar,
				       const int&         a_nref,
				       const int&         a_nvar,
				       const bool         a_forceNoEBCF){
  define(a_eblgFine, a_eblgCoar, a_nref, a_nvar, a_forceNoEBCF);
}


EBFastFluxRegister::~EBFastFluxRegister(){

}

void EBFastFluxRegister::define(const DisjointBoxLayout& a_gridsFine,
				const DisjointBoxLayout& a_gridsCoar,
				const EBISLayout&        a_ebislFine,
				const EBISLayout&        a_ebislCoar,
				const ProblemDomain&     a_domainCoar,
				const int&               a_nref,
				const int&               a_nvar,
				const EBIndexSpace*      ebisPtr,
				const bool               a_forceNoEBCF){
  ProblemDomain domainFine = refine(a_domainCoar, a_nref);
  EBLevelGrid eblgCoar(a_gridsCoar, a_ebislCoar, a_domainCoar);
  EBLevelGrid eblgFine(a_gridsFine, a_ebislFine,   domainFine);

  // Call other version. 
  define(eblgFine, eblgCoar, a_nref, a_nvar, a_forceNoEBCF);

}

void EBFastFluxRegister::define(const EBLevelGrid&       a_eblgFine,
				const EBLevelGrid&       a_eblgCoar,
				const int&               a_refRat,
				const int&               a_nvar,
				const bool               a_forceNoEBCF){
  CH_TIME("EBFastFluxRegister::define");
  
  clear();
  m_refRat   = a_refRat;
  m_nComp    = a_nvar;
  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;
  if (!m_eblgFine.coarsenable(a_refRat)){
    MayDay::Error("EBFastFluxRegister::define: dbl not coarsenable by refrat");
  }

  coarsen(m_eblgCoFi, a_eblgFine, a_refRat);
  m_eblgCoFi.getEBISL().setMaxRefinementRatio(a_refRat, m_eblgCoFi.getEBIS());
  m_reverseCopier.ghostDefine(m_eblgCoFi.getDBL(), m_eblgCoar.getDBL(),m_eblgCoar.getDomain(),IntVect::Unit);

  m_nrefdmo = pow(m_refRat, SpaceDim-1);;

  m_levelFluxReg = new LevelFluxRegister();

  //  pout() << "before regular define" << endl;
  m_levelFluxReg->define(a_eblgFine.getDBL(),
			 a_eblgCoar.getDBL(),
			 a_eblgFine.getDomain(),
			 a_refRat, a_nvar);
  if (a_forceNoEBCF){
    m_hasEBCF = false;
  }
  else{
    m_hasEBCF = computeHasEBCF();
  }

  //if no EBCF--nothing happens here but calls to level flux register
  if (m_hasEBCF){
    defineMasks();
    fastDefineSetsAndIterators();
    defineBuffers();
  }
  
  //set all registers to zero to start.
  setToZero();
  m_isDefined = true;
}

void EBFastFluxRegister::fastDefineSetsAndIterators(){
  CH_TIME("EBFsatFluxRegister::fastDefineSetsAndIterators");
  
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

void EBFastFluxRegister::defineMasks(){
  CH_TIME("EBFastFluxRegister::defineMasks");

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
    //    coFiBoxes[i] &= coarDomain;
  }

  const DisjointBoxLayout grownCoFiGrid(coFiBoxes, coFiProcs);
  BoxLayoutData<FArrayBox> coFiMask(grownCoFiGrid, ncomp);

  // Go through and set ghost cells on the outside of the gridsCoFi to 1,
  // all valid cells and ghost cells "inside" the grid are set to zero
  for (int idir = 0; idir < SpaceDim; idir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const int iindex = index(idir, sit());

      // Define the coar mask in this direction.
      LevelData<FArrayBox> coarMask;
      LayoutData<IntVectSet>& coarCFIVS = m_coarCFIVS[iindex];
      coarMask.define(gridsCoar, ncomp, IntVect::Zero);
      coarCFIVS.define(gridsCoar);

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

      // Go through the mask and set the CFIVS
      for (DataIterator dit = gridsCoar.dataIterator(); dit.ok(); ++dit){
	const Box coarBox = gridsCoar[dit()];
	IntVectSet& cfivs = coarCFIVS[dit()];

	cfivs.makeEmpty();
	for (BoxIterator bit(coarBox); bit.ok(); ++bit){
	  if(coarMask[dit()](bit(), icomp) > 0.0) cfivs |= bit();
	}
      }
    }
  }
}
#include "CD_NamespaceFooter.H"
