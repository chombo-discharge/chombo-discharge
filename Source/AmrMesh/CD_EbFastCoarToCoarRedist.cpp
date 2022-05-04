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
#include <CD_Timer.H>
#include <CD_EbFastCoarToCoarRedist.H>
#include <CD_NamespaceHeader.H>

EbFastCoarToCoarRedist::EbFastCoarToCoarRedist() : EBCoarToCoarRedist() {}

EbFastCoarToCoarRedist::~EbFastCoarToCoarRedist() {}

void
EbFastCoarToCoarRedist::fastDefine(const EBLevelGrid& a_eblgFine,
                                   const EBLevelGrid& a_eblgCoar,
                                   const int&         a_refRat,
                                   const int&         a_nComp,
                                   const int&         a_redistRad)
{
  CH_TIME("EbFastCoarToCoarRedist::define");

  Timer timer("EbFastCoarToCoarRedist::fastDefine");

  // TLDR: The coar-to-coar redistribution utility removes mass that was redistributed to the coarse level from invalid regions, i.e.
  //       regions on the coarse level that is covered by the fine grid. This happens if there is a coarse-fine interface near the
  //       embedded boundary.
  //
  //       The default Chombo version is not usable at scale because it computes a global version of the fine-grid cells (which can
  //       range in the billions). This version avoids that and will therefore be much faster.

  // Usual Chombo stuff
  m_isDefined  = true;
  m_nComp      = a_nComp;
  m_refRat     = a_refRat;
  m_redistRad  = a_redistRad;
  m_domainCoar = a_eblgCoar.getDomain().domainBox();
  m_gridsCoar  = a_eblgCoar.getDBL();
  m_ebislCoar  = a_eblgCoar.getEBISL();

  // Make a mask == true on the coarsened fine grid and copy it to the coarse grid. This creates a BaseFab<bool> mask
  // on the coarse grid where bool=true means that the coarse cell is covered by a finer grid. Note that we include ghost cells
  // in the mask, permitting us to avoid having anything "global". We only include the covering IVS up to a certain range (the
  // redistribution radius).
  timer.startEvent("Define mask");
  constexpr int maskComp    = 0;
  constexpr int numMaskComp = 1;

  DisjointBoxLayout gridsCoFi;
  coarsen(gridsCoFi, a_eblgFine.getDBL(), m_refRat);

  LevelData<BaseFab<bool>> coFiMask(gridsCoFi, numMaskComp, IntVect::Zero);
  LevelData<BaseFab<bool>> coarMask(m_gridsCoar, numMaskComp, m_redistRad * IntVect::Unit);

  for (DataIterator dit(gridsCoFi); dit.ok(); ++dit) {
    coFiMask[dit()].setVal(true);
  }
  for (DataIterator dit(m_gridsCoar); dit.ok(); ++dit) {
    coarMask[dit()].setVal(false);
  }
  timer.stopEvent("Define mask");

  timer.startEvent("Mask exchange");
  coFiMask.copyTo(coarMask);
  coarMask.exchange();
  timer.stopEvent("Mask exchange");

  // Make sets. Recall that m_setsCoar are is the set of cells on the coarse level that will redistribute to this box. So,
  // we are looking for cut-cells within a slightly larger box (grown by the redistribution radius) that are also covered
  // by a finer grid level.
  timer.startEvent("Define sets");
  m_setsCoar.define(m_gridsCoar);
  for (DataIterator dit(m_gridsCoar); dit.ok(); ++dit) {
    const EBISBox& ebisBox = m_ebislCoar[dit()];
    const Box      cellBox = m_gridsCoar[dit()];

    const bool isAllRegular = ebisBox.isAllRegular();
    const bool isAllCovered = ebisBox.isAllCovered();
    const bool isIrregular  = !isAllRegular && !isAllCovered;

    m_setsCoar[dit()].makeEmpty();

    if (isIrregular) {

      // Make the slightly larger box.
      Box gridBox = cellBox;
      gridBox.grow(m_redistRad);
      gridBox &= m_domainCoar;

      // coarMask holds a local view of the fine level grid. If the mask is true, it means that we found a cut-cell
      // that will redistribute from an invalid region and into this box.
#if 0 // Old kernel. Should be fast, but the below kernel should be even faster. 
      DenseIntVectSet coveredIVS(gridBox, false);
      const BaseFab<bool>& mask = coarMask[dit()];
      for (BoxIterator bit(mask.box()); bit.ok(); ++bit){
	const IntVect iv = bit();
	if(mask(iv, maskComp)) coveredIVS |= iv;
      }
    
      m_setsCoar[dit()]  = m_ebislCoar[dit()].getIrregIVS(gridBox);
      m_setsCoar[dit()] &= IntVectSet(coveredIVS);
#else // This version SHOULD be faster than the one above.
      const IntVectSet     irregIVS = ebisBox.getIrregIVS(gridBox);
      const BaseFab<bool>& mask     = coarMask[dit()];

      for (IVSIterator ivsIt(irregIVS); ivsIt.ok(); ++ivsIt) {
        const IntVect iv = ivsIt();
        if (mask(iv, maskComp)) {
          m_setsCoar[dit()] |= iv;
        }
      }
#endif
    }
  }
  timer.stopEvent("Define sets");

  // Define data
  timer.startEvent("Define data");
  EBCoarToCoarRedist::defineDataHolders();
  EBCoarToCoarRedist::setToZero();
  timer.stopEvent("Define data");

  // Below is a debugging hook, intended to make sure the sets definition is the same in both cases.
  ParmParse pp("EbFastCoarToCoarRedist");
  bool      profile = false;
  bool      debug   = false;

  pp.query("profile", profile);
  pp.query("debug", debug);

  if (profile) {
    timer.eventReport(pout(), false);
  }

  if (debug) {
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

    if (diffSet1.numPts() > 0 || diffSet2.numPts() > 0) {
      MayDay::Warning("EbFastCoarToCoarRedist::fastDefine -- sets do not match");

      if (procID() == uniqueProc(SerialTask::compute))
        std::cout << "diffSet1 = " << diffSet1 << std::endl;
      if (procID() == uniqueProc(SerialTask::compute))
        std::cout << "diffSet2 = " << diffSet2 << std::endl;
    }
  }
}

void
EbFastCoarToCoarRedist::gatherBroadcast(IntVectSet& a_set)
{
  CH_TIME("EbFastCoarToCoarRedist::gatherBroadcast");
#ifdef CH_MPI
  Vector<IntVectSet> procSet;
  const int          destProc = uniqueProc(SerialTask::compute);
  gather(procSet, a_set, destProc);
  a_set.makeEmpty();
  if (procID() == destProc) {
    for (int i = 0; i < procSet.size(); i++) {
      a_set |= procSet[i];
    }
  }
  broadcast(a_set, destProc);
#endif
}

void
EbFastCoarToCoarRedist::gatherCoarSet(IntVectSet& a_coarSet)
{
  CH_TIME("EbFastCoarToCoarRedist::gatherCoarSet");

  a_coarSet.makeEmpty();
  for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit) {
    a_coarSet |= m_setsCoar[dit()];
  }
  this->gatherBroadcast(a_coarSet);
}

#include <CD_NamespaceFooter.H>
