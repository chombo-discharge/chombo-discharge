/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbFastCoarToCoarRedist.cpp
  @brief  Implementation of EbFastCoarToCoarRedist
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_EbFastCoarToCoarRedist.H>
#include <CD_NamespaceHeader.H>

#define EBFASTC2C_DEBUG 0

EbFastCoarToCoarRedist::EbFastCoarToCoarRedist() : EBCoarToCoarRedist(){

}

EbFastCoarToCoarRedist::~EbFastCoarToCoarRedist(){

}

void EbFastCoarToCoarRedist::fastDefine(const EBLevelGrid& a_eblgFine,
					const EBLevelGrid& a_eblgCoar,
					const int&         a_refRat,
					const int&         a_nComp,
					const int&         a_redistRad){
  CH_TIME("EbFastCoarToCoarRedist::define");

  // TLDR: The coar-to-coar redistribution utility removes mass that was redistributed to the coarse level from invalid regions, i.e.
  //       regions on the coarse level that is covered by the fine grid. This happens if there is a coarse-fine interface near the
  //       embedded boundary.
  //
  //       The default Chombo version 



  // Usual Chombo stuff
  m_isDefined  = true;
  m_nComp      =  a_nComp;
  m_refRat     =  a_refRat;
  m_redistRad  =  a_redistRad;
  m_domainCoar =  a_eblgCoar.getDomain().domainBox();
  m_gridsCoar  =  a_eblgCoar.getDBL();
  m_ebislCoar  =  a_eblgCoar.getEBISL();

  // Make a mask == true on the coarsened fine grid and copy it to the coarse grid. Include ghost cells of radius
  // m_redistRad so we can build a local, grown view of the covering IVS
  constexpr int maskComp    = 0;
  constexpr int numMaskComp = 1;
  
  DisjointBoxLayout gridsCoFi;
  coarsen(gridsCoFi, a_eblgFine.getDBL(), m_refRat);
  
  LevelData<BaseFab<bool> > coFiMask(gridsCoFi,   numMaskComp, IntVect::Zero);
  LevelData<BaseFab<bool> > coarMask(m_gridsCoar, numMaskComp, m_redistRad*IntVect::Unit);
  
  for (DataIterator dit(gridsCoFi); dit.ok(); ++dit){
    coFiMask[dit()].setVal(true);
  }
  for (DataIterator dit(m_gridsCoar); dit.ok(); ++dit){
    coarMask[dit()].setVal(false);
  }
  coFiMask.copyTo(coarMask);
  coarMask.exchange();

  // Make sets. Recall that m_setsCoar are essentially the cut-cells on the coarse level that are covered by a finer level. 
  m_setsCoar.define(m_gridsCoar);
  for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit){
    Box gridBox = m_gridsCoar.get(dit());
    gridBox.grow(m_redistRad);
    gridBox &= m_domainCoar;

    // coarMask holds a local view of the fine level grid, including ghost cells. Make that view into
    // something that IntVectSet can use
    DenseIntVectSet coveredIVS(gridBox, false);
    const BaseFab<bool>& mask = coarMask[dit()];
    for (BoxIterator bit(mask.box()); bit.ok(); ++bit){
      const IntVect iv = bit();
      if(mask(iv, maskComp)) coveredIVS |= iv;
    }
    
    m_setsCoar[dit()]  = m_ebislCoar[dit()].getIrregIVS(gridBox);
    m_setsCoar[dit()] &= IntVectSet(coveredIVS);    
  }

  // Define data
  EBCoarToCoarRedist::defineDataHolders();
  EBCoarToCoarRedist::setToZero();

  // Below is a debugging hook, intended to make sure the sets definition is the same in both cases. 
  ParmParse pp("EbFastCoarToCoarRedist");
  bool debug = false;
  pp.query("debug", debug);

  if(debug) {

    // Gather m_setsCoar on all ranks using both the new and old define function. 
    IntVectSet newDefineSet;
    IntVectSet oldDefineSet;
    
    this->gatherCoarSet(newDefineSet);
    EBCoarToCoarRedist::define(a_eblgFine, a_eblgCoar, a_refRat, a_nComp, a_redistRad);
    this->gatherCoarSet(oldDefineSet);

    // Compute the intersection of the sets and print the difference between them. If the difference is non-zero then
    // we have a bug. 
    const IntVectSet diffSet1 = newDefineSet - oldDefineSet;
    const IntVectSet diffSet2 = oldDefineSet - newDefineSet;

    if(diffSet1.numPts() > 0 || diffSet2.numPts() > 0){
      MayDay::Warning("EbFastCoarToCoarRedist::fastDefine -- sets do not match");
      
      if(procID() == uniqueProc(SerialTask::compute)) std::cout << "diffSet1 = " << diffSet1 << std::endl;
      if(procID() == uniqueProc(SerialTask::compute)) std::cout << "diffSet2 = " << diffSet2 << std::endl;
    }
  }
}

void EbFastCoarToCoarRedist::gatherBroadcast(IntVectSet& a_set){
  CH_TIME("EbFastCoarToCoarRedist::gatherBroadcast");
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

void EbFastCoarToCoarRedist::gatherCoarSet(IntVectSet& a_coarSet){
  CH_TIME("EbFastCoarToCoarRedist::gatherCoarSet");

  a_coarSet.makeEmpty();
  for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit){
    a_coarSet |= m_setsCoar[dit()];
  }
  this->gatherBroadcast(a_coarSet);
}

#include <CD_NamespaceFooter.H>
