/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   EbFastCoarToFineRedist.cpp
  @brief  Implementation of EbFastCoarToFineRedist.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EbFastCoarToFineRedist.H>
#include <CD_NamespaceHeader.H>

#define EBFASTC2F_DEBUG 0

EbFastCoarToFineRedist::EbFastCoarToFineRedist() : EBCoarToFineRedist(){

}

EbFastCoarToFineRedist::~EbFastCoarToFineRedist(){

}

void EbFastCoarToFineRedist::define(const EBLevelGrid&                      a_eblgFine,
				    const EBLevelGrid&                      a_eblgCoar,
				    const LayoutData<Vector<LayoutIndex> >& a_neighborsFine,
				    const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar,
				    const int&                              a_nref,
				    const int&                              a_nvar,
				    const int&                              a_redistRad){
  CH_TIME("EbFastCoarToFineRedist::define");

  //from here we can assume the redistRad == 1
  m_isDefined  = true;
  m_nComp      = a_nvar;
  m_refRat     = a_nref;
  m_domainCoar = a_eblgCoar.getDomain().domainBox();
  m_gridsFine  = a_eblgFine.getDBL();
  m_gridsCoar  = a_eblgCoar.getDBL();

  m_ebislCoar  = a_eblgCoar.getEBISL();
  m_redistRad  = a_redistRad;

  //created the coarsened fine layout
  m_gridsCedFine = DisjointBoxLayout();
  coarsen(m_gridsCedFine, m_gridsFine, m_refRat);
  const EBIndexSpace* const ebisPtr = a_eblgFine.getEBIS();
  CH_assert(ebisPtr->isDefined());
  int nghost = 3*m_redistRad;
  ebisPtr->fillEBISLayout(m_ebislCedFine, m_gridsCedFine, m_domainCoar, nghost);
  m_ebislCedFine.setMaxRefinementRatio(m_refRat, ebisPtr);

  //define the intvectsets over which the objects live
  m_setsCoar.define(m_gridsCoar);
  m_setsCedFine.define(m_gridsCedFine);

  // Make sets
  IntVectSet globalCFIVS;
  this->makeCedFineSets(globalCFIVS, a_neighborsFine);
  this->makeCoarSets(   globalCFIVS, a_neighborsCoar);

  // Define data holders and reset buffers
  this->defineDataHolders();
  this->setToZero();

#if EBFASTFC2F_DEBUG // This debugging hook calls the original function and checks that the sets completely overlap. 
  IntVectSet newCedFineSet, oldCedFineSet;
  IntVectSet newCoarSet,    oldCoarSet;
  
  this->gatherSetsCedFine(newCedFineSet);
  this->gatherSetsCoar(   newCoarSet);

  // Call old define
  EBCoarToFineRedist::define(a_eblgFine, a_eblgCoar, a_nref, a_nvar, a_redistRad);
  this->gatherSetsCedFine(oldCedFineSet);
  this->gatherSetsCoar(   oldCoarSet);

  const IntVectSet diffSet1 = newCedFineSet - oldCedFineSet;
  const IntVectSet diffSet2 = oldCedFineSet - newCedFineSet;
  const IntVectSet diffSet3 = newCoarSet    - oldCoarSet;
  const IntVectSet diffSet4 = oldCoarSet    - newCoarSet;

  if(diffSet1.numPts() != 0) MayDay::Abort("EbFastCoarToFineRedist::define - diffSet1 not empty");
  if(diffSet2.numPts() != 0) MayDay::Abort("EbFastCoarToFineRedist::define - diffSet2 not empty");
  if(diffSet3.numPts() != 0) MayDay::Abort("EbFastCoarToFineRedist::define - diffSet3 not empty");
  if(diffSet4.numPts() != 0) MayDay::Abort("EbFastCoarToFineRedist::define - diffSet4 not empty");
#endif
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
  gatherBroadcast(a_cfivs);
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
