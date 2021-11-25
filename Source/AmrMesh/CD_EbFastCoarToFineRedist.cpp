/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbFastCoarToFineRedist.cpp
  @brief  Implementation of CD_EbFastCoarToFineRedist.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <NeighborIterator.H>
#include <ParmParse.H>

// Our includes
#include <CD_Timer.H>
#include <CD_EbFastCoarToFineRedist.H>
#include <CD_NamespaceHeader.H>

EbFastCoarToFineRedist::EbFastCoarToFineRedist() : EBCoarToFineRedist(){

}

EbFastCoarToFineRedist::~EbFastCoarToFineRedist(){

}

void EbFastCoarToFineRedist::fastDefine(const EBLevelGrid& a_eblgFine,
					const EBLevelGrid& a_eblgCoar,
					const int&         a_refRat,
					const int&         a_nComp,
					const int&         a_redistRad){
  CH_TIME("EbFastCoarToFineRedist::fastDefine");

  const EBIndexSpace* const ebisPtr = a_eblgFine.getEBIS();
  
  CH_assert(ebisPtr->isDefined());

  Timer timer("EbFastCoarToFineRedist::fastDefine");

  // Regular stuff for base class. From here we can assume the redistRad == 1
  m_isDefined  = true;
  m_nComp      = a_nComp;
  m_refRat     = a_refRat;
  m_domainCoar = a_eblgCoar.getDomain().domainBox();
  m_gridsFine  = a_eblgFine.getDBL();
  m_gridsCoar  = a_eblgCoar.getDBL();
  m_ebislCoar  = a_eblgCoar.getEBISL();
  m_redistRad  = a_redistRad;

  // First, create the coarsened fine layout.
  // Comment, RM: Why does chombo need a buffer with 3*m_redistRad ghost cells in the EBISBoxes? That seems awfully big. 
  //  m_gridsCedFine = DisjointBoxLayout();
  timer.startEvent("Make coarsened grids");
  coarsen(m_gridsCedFine, m_gridsFine, m_refRat);
  ebisPtr->fillEBISLayout(m_ebislCedFine, m_gridsCedFine, m_domainCoar, 3*m_redistRad);
  m_ebislCedFine.setMaxRefinementRatio(m_refRat, ebisPtr);
  timer.stopEvent("Make coarsened grids");

  // Define the sets over which the objects live.
  timer.startEvent("Define sets");
  m_setsCoar.   define(m_gridsCoar);
  m_setsCedFine.define(m_gridsCedFine);
  timer.stopEvent("Define sets");

  // Make sets
  timer.startEvent("Make fine set");
  this->makeCedFineSets();
  timer.stopEvent("Make fine set");    
  timer.startEvent("Make coar set");  
  this->makeCoarSets();  
  timer.stopEvent("Make coar set");  

  // Define data holders and reset buffers
  timer.startEvent("Define data");
  this->defineDataHolders();
  this->setToZero();
  timer.stopEvent("Define data");

  ParmParse pp("EbFastCoarToFineRedist");
  bool profile = false;
  bool debug   = false;
  pp.query("debug",   debug);
  pp.query("profile", profile);

  if(profile){
    timer.eventReport(pout(), false);
  }

  if(debug){
    IntVectSet newCedFineSet, oldCedFineSet;
    IntVectSet newCoarSet,    oldCoarSet;
  
    this->gatherSetsCedFine(newCedFineSet);
    this->gatherSetsCoar(   newCoarSet);

    // Call the old define function. 
    EBCoarToFineRedist::define(a_eblgFine, a_eblgCoar, a_refRat, a_nComp, a_redistRad);
    this->gatherSetsCedFine(oldCedFineSet);
    this->gatherSetsCoar(   oldCoarSet);

    // Compare the sets. If they're not empty, the old and new define functions don't do the same thing
    // and we have a bug. 
    const IntVectSet diffSet1 = newCedFineSet - oldCedFineSet;
    const IntVectSet diffSet2 = oldCedFineSet - newCedFineSet;
    const IntVectSet diffSet3 = newCoarSet    - oldCoarSet;
    const IntVectSet diffSet4 = oldCoarSet    - newCoarSet;

    if(procID() == uniqueProc(SerialTask::compute)){
      if(diffSet1.numPts() != 0) MayDay::Error("EbFastCoarToFineRedist::define - diffSet1 not empty");
      if(diffSet2.numPts() != 0) MayDay::Error("EbFastCoarToFineRedist::define - diffSet2 not empty");
      if(diffSet3.numPts() != 0) MayDay::Error("EbFastCoarToFineRedist::define - diffSet3 not empty");
      if(diffSet4.numPts() != 0) MayDay::Error("EbFastCoarToFineRedist::define - diffSet4 not empty");
    }
  }
}

void EbFastCoarToFineRedist::makeCedFineSets() {
  CH_TIME("EbFastCoarToFineRedist::makeCedFineSets");

  for (DataIterator dit(m_gridsCedFine); dit.ok(); ++dit){
    const Box&     cellBox = m_gridsCedFine[dit()];
    const EBISBox& ebisBox = m_ebislCedFine[dit()];

    const bool isAllRegular = ebisBox.isAllRegular();
    const bool isAllCovered = ebisBox.isAllCovered();
    const bool isIrregular  = !isAllRegular && !isAllCovered;

    m_setsCedFine[dit()].makeEmpty();

    if(isIrregular){
      // Make the grown box
      Box grownBox = grow(cellBox, m_redistRad);
      grownBox    &= m_domainCoar;

      const IntVectSet irregIVS = ebisBox.getIrregIVS(grownBox);      

      // We want all cut-cells that will redistribute to this patch. This includes every coarse cut-cell in a radius
      // m_redistRad outside the coarse-fine interface. Run with a DenseIntVectSet here because it will have faster
      // intersection functions. When we make the setsCedFine, convert to a regular IntVectSet. 
      DenseIntVectSet coarCellsRedist(grownBox, true);

      // Subtract of the cells in this box, and then subtract off cells in the neighbor boxes. 
      coarCellsRedist -= cellBox;
    
      for (LayoutIterator lit = m_gridsCedFine.layoutIterator(); lit.ok(); ++lit){
	const Box neighborBox = m_gridsCedFine[lit()];
	const Box overlapBox  = neighborBox & grownBox;

	if(!overlapBox.isEmpty()){
	  coarCellsRedist -= overlapBox;
	}
      }

      // coarCellsRedist now contain the cells on the coarse level that are NOT covered by a finer level, and that are within range
      // m_redistRadius from the current grid box. These are the only cells that redistribute over the coarse-fine interface and
      // into this (fine-level) grid box. 
      for (IVSIterator ivsIt(irregIVS); ivsIt.ok(); ++ivsIt){
	const IntVect iv = ivsIt();
	if(coarCellsRedist[iv]){
	  m_setsCedFine[dit()] |= iv;
	}
      }
    }
  }
}

void EbFastCoarToFineRedist::makeCoarSets(){
  CH_TIME("EbFastCoarToFineRedist::makeCoarSets");

  // The coar sets is the set of cells on each coarse-grid patch that can redistribute to the fine grid. Unfortunately,
  // this is slightly more complicated than the other set because we need to intersect the coarse-grid cut-cells with
  // the entire refinement boundary (appropriately grown by m_redistRad). The original Chombo version did a local reconstruction
  // of the refinement boundary, which killed performance. Instead, we coarsen the fine grid and create a mask with value = 1
  // in cells that lie in the refinment bounary (again, grown by m_redistRad). This mask is simply copied to the coarse-level,
  // which provides us a view of which cells on the coarse level that lie on the refinement boundary. Note that we use a BoxLayoutData
  // to achieve this (I would have used BaseFab<bool>, but it doesn't have an incrementation operator). So, rather than copyTo we set
  // the mask to 1 in the refinement boundary and use addTo to increment the mask on the coarse level.
  //
  // NOTE: We could, alternatively, simply gather m_setsCedFine globally and compute the intersection of that sets with the irregular cells
  //       in each coarse-grid box. But I prefer this way to avoid a global reduction of something whose size is unknown. Nonetheless, I'm leaving
  //       the alternative method in place, should we ever see that the default method becomes a bottleneck. 

#if 1 // Default method
  constexpr Real     zero   = 0.0;
  constexpr Real     one    = 1.0;    
  constexpr int      comp   = 0;
  constexpr int      nComp  = 1;
  const     Interval interv = Interval(comp, comp);

  const IntVect dstGhost =             IntVect::Zero;
  const IntVect srcGhost = m_redistRad*IntVect::Unit;
  
  LevelData<FArrayBox> coarMask(m_gridsCoar,    nComp, dstGhost);
  LevelData<FArrayBox> coFiMask(m_gridsCedFine, nComp, srcGhost);

  // Set the coarse mask to be zero everywhere. 
  for (DataIterator dit(m_gridsCoar); dit.ok(); ++dit){
    coarMask[dit()].setVal(zero);
  }

  // On the coarsened fine mask we iterate through all the cells and set the mask to be = 1 outside the refinement boundary (up to m_redistRad away).
  for (DataIterator dit(m_gridsCedFine); dit.ok(); ++dit){
    FArrayBox& mask       = coFiMask[dit()];
    const Box& cellBox    = m_gridsCedFine[dit()]; // This is the unghosted box.
    const Box& ghostedBox = mask.box();            // This is the ghosted box, i.e. cellBox but grown by srcGhost.

    mask.setVal(one);
    mask.setVal(zero, cellBox, comp, nComp);

    // This will 
    NeighborIterator nit(m_gridsCedFine);
    for (nit.begin(dit()); nit.ok(); ++nit){
      const Box overlapBox  = m_gridsCedFine[nit()] & ghostedBox;

      if(!overlapBox.isEmpty()){
	mask.setVal(zero, overlapBox, comp, nComp);
      }
    }
  }

  // Add the result to the coarse grid. This should add from ghost cells and into the valid region.
  Copier copier;
  copier.ghostDefine(m_gridsCedFine, m_gridsCoar, m_domainCoar, srcGhost);
  coFiMask.copyTo(interv, coarMask, interv, copier, LDaddOp<FArrayBox>());

  // Create the sets. 
  for (DataIterator dit(m_gridsCoar); dit.ok(); ++dit){
    const Box        cellBox = m_gridsCoar[dit()];
    const EBISBox&   ebisBox = m_ebislCoar[dit()];
    const FArrayBox& mask    = coarMask   [dit()];

    const bool isAllRegular = ebisBox.isAllRegular();
    const bool isAllCovered = ebisBox.isAllCovered();
    const bool isIrregular  = !isAllRegular && !isAllCovered;
    
    m_setsCoar[dit()].makeEmpty();

    // Recall: The mask has a value > 0 in all cells that are a distance m_redistRad cells from the refinement boundary. This includes regular
    //         cells but those don't redistribute. So, here we fetch all the valid coarse-grid cells that redistribute over the refinement boundary
    //         and into the fine grid.
    if(isIrregular){
      const IntVectSet irregIVS = ebisBox.getIrregIVS(cellBox);

      for (IVSIterator ivsIt(irregIVS); ivsIt.ok(); ++ivsIt){
	const IntVect iv = ivsIt();

	if(mask(iv, comp) > zero){
	  m_setsCoar[dit()] |= iv;
	}
      }
    }
  }  
#else
  IntVectSet globalIVS;

  for (DataIterator dit(m_gridsCedFine); dit.ok(); ++dit){
    globalIVS |= m_setsCedFine[dit()];
  }

  this->gatherBroadcast(globalIVS);

  for (DataIterator dit(m_gridsCoar); dit.ok(); ++dit){
    m_setsCoar[dit()] = globalIVS & m_ebislCoar[dit()].getIrregIVS(m_gridsCoar[dit()]);
  }
#endif
}

void EbFastCoarToFineRedist::gatherBroadcast(IntVectSet& a_set){
  CH_TIME("EBFastCoarToFine::gatherBroadcast");
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

void EbFastCoarToFineRedist::gatherSetsCedFine(IntVectSet& a_setsCedFine){
  CH_TIME("EbFastCoarToFineRedist::gatherSetsCedFine");

  a_setsCedFine.makeEmpty();
  for (DataIterator dit = m_gridsCedFine.dataIterator(); dit.ok(); ++dit){
    a_setsCedFine |= m_setsCedFine[dit()];
  }
  this->gatherBroadcast(a_setsCedFine);
}

void EbFastCoarToFineRedist::gatherSetsCoar(IntVectSet& a_setsCoar){
  CH_TIME("EbFastCoarToFineRedist::gatherSetsCoar");

  a_setsCoar.makeEmpty();
  for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit){
    a_setsCoar |= m_setsCoar[dit()];
  }
  
  this->gatherBroadcast(a_setsCoar);
}

#include <CD_NamespaceFooter.H>
