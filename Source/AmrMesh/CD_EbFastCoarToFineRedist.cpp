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
#include <ParmParse.H>

// Our includes
#include <CD_Timer.H>
#include <CD_EbFastCoarToFineRedist.H>
#include <CD_NamespaceHeader.H>

#define EBFASTC2F_DEBUG 0

EbFastCoarToFineRedist::EbFastCoarToFineRedist() : EBCoarToFineRedist(){

}

EbFastCoarToFineRedist::~EbFastCoarToFineRedist(){

}

void EbFastCoarToFineRedist::fastDefine(const EBLevelGrid&                      a_eblgFine,
					const EBLevelGrid&                      a_eblgCoar,
					const LayoutData<Vector<LayoutIndex> >& a_neighborsFine,
					const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar,
					const int&                              a_nRef,
					const int&                              a_nVar,
					const int&                              a_redistRad){
  CH_TIME("EbFastCoarToFineRedist::fastDefine");

  const EBIndexSpace* const ebisPtr = a_eblgFine.getEBIS();
  
  CH_assert(ebisPtr->isDefined());

  Timer timer("EbFastCoarToFineRedist::fastDefine");

  // Regular stuff for base class. From here we can assume the redistRad == 1
  m_isDefined  = true;
  m_nComp      = a_nVar;
  m_refRat     = a_nRef;
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

  // Make sets.
  timer.startEvent("Make sets");  
  IntVectSet globalCFIVS;
  this->makeCedFineSets(globalCFIVS, a_neighborsFine);
  this->makeCoarSets(   globalCFIVS, a_neighborsCoar);
  timer.stopEvent("Make sets");    

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

    // Call old define
    EBCoarToFineRedist::define(a_eblgFine, a_eblgCoar, a_nRef, a_nVar, a_redistRad);
    this->gatherSetsCedFine(oldCedFineSet);
    this->gatherSetsCoar(   oldCoarSet);

    const IntVectSet diffSet1 = newCedFineSet - oldCedFineSet;
    const IntVectSet diffSet2 = oldCedFineSet - newCedFineSet;
    const IntVectSet diffSet3 = newCoarSet    - oldCoarSet;
    const IntVectSet diffSet4 = oldCoarSet    - newCoarSet;

    if(procID() == uniqueProc(SerialTask::compute)){
      if(diffSet1.numPts() != 0) MayDay::Abort("EbFastCoarToFineRedist::define - diffSet1 not empty");
      if(diffSet2.numPts() != 0) MayDay::Abort("EbFastCoarToFineRedist::define - diffSet2 not empty");
      if(diffSet3.numPts() != 0) MayDay::Abort("EbFastCoarToFineRedist::define - diffSet3 not empty");
      if(diffSet4.numPts() != 0) MayDay::Abort("EbFastCoarToFineRedist::define - diffSet4 not empty");
    }
  }
}

void EbFastCoarToFineRedist::makeCedFineSets(IntVectSet& a_cfivs, const LayoutData<Vector<LayoutIndex> >& a_neighborsFine){
  CH_TIME("EbFastCoarToFineRedist::makeCedFineSets");

  a_cfivs.makeEmpty();
  for (DataIterator dit = m_gridsCedFine.dataIterator(); dit.ok(); ++dit){
    const Box& box = m_gridsCedFine.get(dit());
    Box grownBox   = grow(box, m_redistRad);
    grownBox      &= m_domainCoar;

    IntVectSet boxCFIVS(grownBox);
    boxCFIVS -= box;

    m_setsCedFine[dit()]  = m_ebislCedFine[dit()].getIrregIVS(grownBox);
    m_setsCedFine[dit()] -= box;

    // Subtract off cells from the coarsened fine grid
    const Vector<LayoutIndex>& neighbors = a_neighborsFine[dit()];
    for (int i = 0; i < neighbors.size(); i++){
      Box cedBox = m_gridsFine.get(neighbors[i]);
      cedBox.coarsen(m_refRat);

      boxCFIVS -= cedBox;
      m_setsCedFine[dit()] -= cedBox;
    }

    a_cfivs |= boxCFIVS;
  }
  
  this->gatherBroadcast(a_cfivs);
}


void EbFastCoarToFineRedist::makeCoarSets(const IntVectSet& a_cfivs, const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar){
  CH_TIME("EbFastCoarToFineRedist::makeCoarSets");

  // We have a global view of the CFIVS so there's not much to this...
  for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit){
    const Box& coarBox = m_gridsCoar.get(dit());
    m_setsCoar[dit()]  = m_ebislCoar[dit()].getIrregIVS(coarBox);
    m_setsCoar[dit()] &= a_cfivs;
  }
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
