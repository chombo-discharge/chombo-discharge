/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbFastFineToCoarRedist.cpp
  @brief  Implementation of EbFastFineToCoarRedist
  @author Robert Marskar
*/

// Chombo includes
#include <NeighborIterator.H>
#include <ParmParse.H>

// Our includes
#include <CD_Timer.H>
#include <CD_EbFastFineToCoarRedist.H>
#include <CD_NamespaceHeader.H>

#define EBFASTF2C_DEBUG 0

EbFastFineToCoarRedist::EbFastFineToCoarRedist() : EBFineToCoarRedist(){

}

EbFastFineToCoarRedist::~EbFastFineToCoarRedist(){

}

void EbFastFineToCoarRedist::fastDefine(const EBLevelGrid&                      a_eblgFine,
					const EBLevelGrid&                      a_eblgCoar,
					const int&                              a_nRef,
					const int&                              a_nVar,
					const int&                              a_redistRad){
  CH_TIME("EbFastFineToCoarRedist::define");

  const EBIndexSpace* ebisPtr = a_eblgFine.getEBIS();
  CH_assert(ebisPtr->isDefined());  

  Timer timer("EbFastFineToCoarRedist");
  
  m_isDefined = true;
  m_nComp     = a_nVar;
  m_refRat    = a_nRef;
  m_redistRad = a_redistRad;
  m_domainCoar= a_eblgCoar.getDomain().domainBox();
  m_gridsFine = a_eblgFine.getDBL();
  m_gridsCoar = a_eblgCoar.getDBL();
  m_ebislFine = a_eblgFine.getEBISL();
  m_ebislCoar = a_eblgCoar.getEBISL();
  
  //created the coarsened fine layout
  timer.startEvent("Make refined grids");
  m_gridsRefCoar = DisjointBoxLayout();
  refine(m_gridsRefCoar, m_gridsCoar, m_refRat);
  int nghost = 3*m_redistRad;
  Box domainFine = refine(m_domainCoar, m_refRat);
  ebisPtr->fillEBISLayout(m_ebislRefCoar, m_gridsRefCoar, domainFine, nghost);
  m_ebislRefCoar.setMaxCoarseningRatio(m_refRat,ebisPtr);
  timer.stopEvent("Make refined grids");    

  //define the intvectsets over which the objects live
  timer.startEvent("Define sets");
  m_setsFine.   define(m_gridsFine);
  m_setsRefCoar.define(m_gridsCoar);
  timer.stopEvent("Define sets");  

  // The new define function creates a local view of the inside CFIVS in each fine box, and then
  // copies that mask to the refined coarse grids. The refCoar mask also holds updated ghost cells
  // so there is no need for extra communications beyond the mask copy and exchange
  timer.startEvent("Define masks");
  constexpr int numMaskComp = 1;
  
  LevelData<BaseFab<bool> > fineShellMask;
  LevelData<BaseFab<bool> > refCoarShellMask;
  
  fineShellMask.   define(m_gridsFine,    numMaskComp,             IntVect::Zero);
  refCoarShellMask.define(m_gridsRefCoar, numMaskComp, m_redistRad*IntVect::Unit);
  timer.stopEvent("Define masks");  

  // Make the fine mask and copy it to the mask on the refined coarse grids
  timer.startEvent("Make fine mask");
  this->makeFineMask(fineShellMask);
  timer.stopEvent("Make fine mask");
  timer.startEvent("Make coar mask");
  for (DataIterator dit = m_gridsRefCoar.dataIterator(); dit.ok(); ++dit){
    refCoarShellMask[dit()].setVal(false);
  }
  fineShellMask.copyTo(refCoarShellMask);
  refCoarShellMask.exchange();
  timer.stopEvent("Make coar mask");    

  timer.startEvent("Define fine set");
  this->makeFineSets(fineShellMask);
  timer.stopEvent("Define fine set");
  timer.startEvent("Define coar set");  
  this->makeCoarSets(refCoarShellMask);
  timer.stopEvent("Define coar set");    

  timer.startEvent("Define data");
  this->defineDataHolders();
  this->setToZero();
  timer.stopEvent("Define data");  

  ParmParse pp("EbFastFineToCoarRedist");

  bool profile = false;
  bool debug   = false;

  pp.get("profile", profile);
  pp.get("debug",   debug);  

  if(profile){
    timer.eventReport(pout(), false);
  }

  if(debug){
    IntVectSet newFineSet,    oldFineSet;
    IntVectSet newRefCoarSet, oldRefCoarSet;
  
    this->gatherSetsFine(newFineSet);
    this->gatherSetsRefCoar(newRefCoarSet);

    // Call old define
    EBFineToCoarRedist::define(a_eblgFine, a_eblgCoar, a_nRef, a_nVar, a_redistRad);
    this->gatherSetsFine(oldFineSet);
    this->gatherSetsRefCoar(oldRefCoarSet);

    const IntVectSet diffSet1 = newFineSet    - oldFineSet;
    const IntVectSet diffSet2 = oldFineSet    - newFineSet;
    const IntVectSet diffSet3 = newRefCoarSet - oldRefCoarSet;
    const IntVectSet diffSet4 = oldRefCoarSet - newRefCoarSet;

    if(diffSet1.numPts() != 0) MayDay::Error("EbFastFineToCoarRedist::define - diffSet1 not empty");
    if(diffSet2.numPts() != 0) MayDay::Error("EbFastFineToCoarRedist::define - diffSet2 not empty");
    if(diffSet3.numPts() != 0) MayDay::Error("EbFastFineToCoarRedist::define - diffSet3 not empty");
    if(diffSet4.numPts() != 0) MayDay::Error("EbFastFineToCoarRedist::define - diffSet4 not empty");
  }
}

void EbFastFineToCoarRedist::makeFineMask(LevelData<BaseFab<bool> >& a_fineShellMask){
  CH_TIME("EbFastFineToCoarRedist::makeFineMask");

  // TLDR: We are trying to figure out which cells in each fine-grid patch redistributes over the coarse-fine interface. To do that
  //       we need to know about all cells on the fine level that are a distance m_redistRad away from the interface. We simply compute
  //       the cells
  constexpr int maskComp = 0;
  
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
    const Box& fineBox = m_gridsFine.get(dit());
    
    // Make the coarse-fine interface. I think it's safe to use a NeighborIterator here because we are only dealing with valid cells.
    Box grownBox = grow(fineBox, 1);
    grownBox &= m_ebislFine.getDomain();
    
    IntVectSet coarseFineInterface(grownBox);
    coarseFineInterface -= fineBox;

    NeighborIterator nit(m_gridsFine);
    for (nit.begin(dit()); nit.ok(); ++nit){
      coarseFineInterface -= m_gridsFine[nit()];
    }
    
    // Make the local shell on the inside of each box by growing the CFISV and then restricting to the fine-grid box. 
    IntVectSet fineShell;
    for (IVSIterator ivsIt(coarseFineInterface); ivsIt.ok(); ++ivsIt){
      const IntVect iv = ivsIt();
      Box b(iv,iv);
      b.grow(m_redistRad);
      fineShell |= IntVectSet(b);
    }
    fineShell &= fineBox;


    a_fineShellMask[dit()].setVal(false);
    for (IVSIterator ivsIt(fineShell); ivsIt.ok(); ++ivsIt){
      const IntVect iv = ivsIt();
      a_fineShellMask[dit()](iv, maskComp) = true;
    }
  }
}

void EbFastFineToCoarRedist::makeFineSets(const LevelData<BaseFab<bool> >& a_fineMask){
  CH_TIME("EbFastFineToCoarRedist::makeFineSets");

  // TLDR: m_setsFine are the set of cut-cells that redistribute to the coarse level. We have made a mask (a_fineMask) which
  //       contains all cells that fall within m_redistRad on the inside of refinement boundary. So all we need to do is check
  //       which of the cut-cells in the fine-grid patches intersect with that mask. 

  constexpr int maskComp = 0;  
  
  for (DataIterator dit(m_gridsFine); dit.ok(); ++dit){
    const Box&           fineBox  = m_gridsFine[dit()];
    const EBISBox&       ebisBox  = m_ebislFine[dit()];
    const BaseFab<bool>& fineMask = a_fineMask [dit()];

    const bool isAllRegular = ebisBox.isAllRegular();
    const bool isAllCovered = ebisBox.isAllCovered();
    const bool isIrregular  = !isAllCovered && !isAllRegular;

    m_setsFine[dit()].makeEmpty();

    if(isIrregular){
      const IntVectSet irregIVS = ebisBox.getIrregIVS(fineBox);
      for (IVSIterator ivsIt(irregIVS); ivsIt.ok(); ++ivsIt){
	const IntVect iv = ivsIt();
	
	if(fineMask(iv, maskComp)){
	  m_setsFine[dit()] |= iv;
	}
      }
    }
  }
}

void EbFastFineToCoarRedist::makeCoarSets(const LevelData<BaseFab<bool> >& a_refCoarMask){
  CH_TIME("EbFastFineToCoarRedist::makeCoarSets");

  constexpr int maskComp = 0;

  const Box domainFine = refine(m_domainCoar, m_refRat);

  for (DataIterator dit(m_gridsRefCoar); dit.ok(); ++dit){
    const Box&           coarBox  = m_gridsRefCoar[dit()];
    const EBISBox&       ebisBox  = m_ebislRefCoar[dit()];
    const BaseFab<bool>& maskCoar = a_refCoarMask [dit()];

    const bool isAllRegular = ebisBox.isAllRegular();
    const bool isAllCovered = ebisBox.isAllCovered();
    const bool isIrregular  = !isAllCovered && !isAllRegular;

    m_setsRefCoar[dit()].makeEmpty();

    if(isIrregular){
      Box grownBox = grow(coarBox, m_redistRad);
      grownBox &= domainFine;
      
      const IntVectSet irregIVS = ebisBox.getIrregIVS(grownBox);

      for (IVSIterator ivsIt(irregIVS); ivsIt.ok(); ++ivsIt){
	const IntVect iv = ivsIt();
	if(maskCoar(iv, maskComp)){
	  m_setsRefCoar[dit()] |= iv;
	}
      }
    }
  }
}

void EbFastFineToCoarRedist::gatherBroadcast(IntVectSet& a_set){
  CH_TIME("EbFastFineToCoarRedist::gatherBroadcast");
#ifdef CH_MPI
  Vector<IntVectSet> procSet;
  const int destProc = uniqueProc(SerialTask::compute);
  gather(procSet, a_set, destProc);
  a_set.makeEmpty();
  if(procID() == destProc){
    for (int i = 0; i < procSet.size(); i++){
      a_set |= procSet[i];
    }
  }
  broadcast(a_set, destProc);
#endif
}

void EbFastFineToCoarRedist::gatherSetsFine(IntVectSet& a_setsFine){
  CH_TIME("EbFastFineToCoarRedist::gatherSetsFine");

  a_setsFine.makeEmpty();
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
    a_setsFine |= m_setsFine[dit()];
  }
  
  this->gatherBroadcast(a_setsFine);
}

void EbFastFineToCoarRedist::gatherSetsRefCoar(IntVectSet& a_setsRefCoar){
  CH_TIME("EbFastFineToCoarRedist::gatherSetsRefCoar");

  a_setsRefCoar.makeEmpty();
  for (DataIterator dit = m_gridsRefCoar.dataIterator(); dit.ok(); ++dit){
    a_setsRefCoar |= m_setsRefCoar[dit()];
  }
  
  this->gatherBroadcast(a_setsRefCoar);
}

#include <CD_NamespaceFooter.H>
