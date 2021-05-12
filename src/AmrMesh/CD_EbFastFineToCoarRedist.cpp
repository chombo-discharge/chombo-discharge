/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_EbFastFineToCoarRedist.cpp
  @brief  Implementation of EbFastFineToCoarRedist
  @author Robert Marskar
*/

// Our includes
#include <CD_EbFastFineToCoarRedist.H>
#include <CD_NamespaceHeader.H>

#define EBFASTF2C_DEBUG 0

EbFastFineToCoarRedist::EbFastFineToCoarRedist() : EBFineToCoarRedist(){

}

EbFastFineToCoarRedist::~EbFastFineToCoarRedist(){

}

void EbFastFineToCoarRedist::define(const EBLevelGrid&                      a_eblgFine,
				    const EBLevelGrid&                      a_eblgCoar,
				    const LayoutData<Vector<LayoutIndex> >& a_neighborsFine,
				    const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar,
				    const int&                              a_nRef,
				    const int&                              a_nVar,
				    const int&                              a_redistRad){
  CH_TIME("EbFastFineToCoarRedist::define");

  const int comp  = 0;
  const int ncomp = 1;
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
  m_gridsRefCoar = DisjointBoxLayout();
  refine(m_gridsRefCoar, m_gridsCoar, m_refRat);

  const EBIndexSpace* ebisPtr = a_eblgFine.getEBIS();
  CH_assert(ebisPtr->isDefined());
  int nghost = 3*m_redistRad;
  Box domainFine = refine(m_domainCoar, m_refRat);
  ebisPtr->fillEBISLayout(m_ebislRefCoar, m_gridsRefCoar, domainFine, nghost);

  m_ebislRefCoar.setMaxCoarseningRatio(m_refRat,ebisPtr);

  //define the intvectsets over which the objects live
  m_setsFine.define(m_gridsFine);
  m_setsRefCoar.define(m_gridsCoar);

  // The new define function creates a local view of the inside CFIVS in each fine box, and then
  // copies that mask to the refined coarse grids. The refCoar mask also holds updated ghost cells
  // so there is no need for extra communications beyond the mask copy and exchange
  LevelData<BaseFab<bool> > fineShellMask;
  LevelData<BaseFab<bool> > refCoarShellMask;

  fineShellMask.define(   m_gridsFine,    ncomp, IntVect::Zero);
  refCoarShellMask.define(m_gridsRefCoar, ncomp, m_redistRad*IntVect::Unit);

  // Make the fine mask and copy it to the mask on the refined coarse grids
  makeFineMask(fineShellMask, a_neighborsFine);
  for (DataIterator dit = m_gridsRefCoar.dataIterator(); dit.ok(); ++dit){
    refCoarShellMask[dit()].setVal(false);
  }
  fineShellMask.copyTo(refCoarShellMask);
  refCoarShellMask.exchange();

  this->makeFineSets(fineShellMask);
  this->makeCoarSets(refCoarShellMask);

  this->defineDataHolders();
  this->setToZero();

#if EBFASTF2C_DEBUG // This debugging hook calls the original function and checks that the sets completely overlap. 
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

  if(diffSet1.numPts() != 0) MayDay::Abort("EbFastFineToCoarRedist::define - diffSet1 not empty");
  if(diffSet2.numPts() != 0) MayDay::Abort("EbFastFineToCoarRedist::define - diffSet2 not empty");
  if(diffSet3.numPts() != 0) MayDay::Abort("EbFastFineToCoarRedist::define - diffSet3 not empty");
  if(diffSet4.numPts() != 0) MayDay::Abort("EbFastFineToCoarRedist::define - diffSet4 not empty");
#endif
}

void EbFastFineToCoarRedist::makeFineMask(LevelData<BaseFab<bool> >&              a_fineShellMask,
					  const LayoutData<Vector<LayoutIndex> >& a_neighborsFine){
  CH_TIME("EbFastFineToCoarRedist::makeFineMask");

  const int comp = 0;
  
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
    const Box& fineBox = m_gridsFine.get(dit());
    
    // Make the inside shell
    Box grownBox = grow(fineBox, 1);
    IntVectSet cfivs(grownBox);
    cfivs -= fineBox;
    const Vector<LayoutIndex>& neighbors = a_neighborsFine[dit()];
    for (int i = 0; i < neighbors.size(); i++){
      cfivs -= m_gridsFine[neighbors[i]];
    }
    
    // Make the local shell by growing the CFISV and restricting to box
    IntVectSet shell;
    for (IVSIterator ivsIt(cfivs); ivsIt.ok(); ++ivsIt){
      const IntVect iv = ivsIt();
      Box b(iv,iv);
      b.grow(m_redistRad);
      shell |= IntVectSet(b);
    }
    shell &= fineBox;


    a_fineShellMask[dit()].setVal(false);
    for (IVSIterator ivsIt(shell); ivsIt.ok(); ++ivsIt){
      const IntVect iv = ivsIt();
      a_fineShellMask[dit()](iv, comp) = true;
    }
  }
}

void EbFastFineToCoarRedist::makeFineSets(const LevelData<BaseFab<bool> >& a_fineMask){
  CH_TIME("EbFastFineToCoarRedist::makeFineSets");

  const int comp = 0;
  
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
    const Box& fineBox = m_gridsFine.get(dit());
    
    // Make mask into IntVect
    IntVectSet localShell;
    const BaseFab<bool>& mask = a_fineMask[dit()];
    for(BoxIterator bit(mask.box()); bit.ok(); ++bit){
      const IntVect iv = bit();
      if(mask(iv, comp)) localShell |= iv;
    }

    m_setsFine[dit()]  = m_ebislFine[dit()].getIrregIVS(fineBox);
    m_setsFine[dit()] &= localShell;
  }
}

void EbFastFineToCoarRedist::makeCoarSets(const LevelData<BaseFab<bool> >& a_refCoarMask){
  CH_TIME("EbFastFineToCoarRedist::makeCoarSets");

  const int comp = 0;

  const Box domainFine = refine(m_domainCoar, m_refRat);

  for (DataIterator dit = m_gridsRefCoar.dataIterator(); dit.ok(); ++dit){
    const Box& refCoarBox = m_gridsRefCoar.get(dit());
    
    // Make mask into IntVect shell
    IntVectSet localShell;
    const BaseFab<bool>& mask = a_refCoarMask[dit()];
    for(BoxIterator bit(mask.box()); bit.ok(); ++bit){
      const IntVect iv = bit();
      if(mask(iv, comp)) localShell |= iv;
    }

    Box grownBox = grow(refCoarBox, m_redistRad);
    grownBox &= domainFine;

    m_setsRefCoar[dit()]  = m_ebislRefCoar[dit()].getIrregIVS(grownBox);
    m_setsRefCoar[dit()] &= localShell;
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
