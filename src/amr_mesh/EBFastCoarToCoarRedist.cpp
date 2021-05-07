/*!
  @file   EBFastCoarToCoarRedist.cpp
  @brief  Implementation of EBFastCoarToCoarRedist
  @author Robert Marskar
*/

#include "EBFastCoarToCoarRedist.H"

#define EBFASTC2C_DEBUG 0

namespace ChomboDischarge {

  EBFastCoarToCoarRedist::EBFastCoarToCoarRedist() : EBCoarToCoarRedist(){

  }

  EBFastCoarToCoarRedist::~EBFastCoarToCoarRedist(){

  }

  void EBFastCoarToCoarRedist::define(const EBLevelGrid&                      a_eblgFine,
				      const EBLevelGrid&                      a_eblgCoar,
				      const LayoutData<Vector<LayoutIndex> >& a_neighborsFine,
				      const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar,
				      const int&                              a_nref,
				      const int&                              a_nvar,
				      const int&                              a_redistRad){
    CH_TIME("EBFastCoarToCoarRedist::define");

    const int comp  = 0;
    const int ncomp = 1;

    m_isDefined  = true;
    m_nComp      =  a_nvar;
    m_refRat     =  a_nref;
    m_redistRad  =  a_redistRad;
    m_domainCoar =  a_eblgCoar.getDomain().domainBox();
    m_gridsCoar  =  a_eblgCoar.getDBL();
    m_ebislCoar  =  a_eblgCoar.getEBISL();


    // Make a mask == true on the coarsened fine grid and copy it to the coarse grid. Include ghost cells of radius
    // m_redistRad so we can build a local, grown view of the covering IVS
    DisjointBoxLayout gridsCedFine;
    coarsen(gridsCedFine, a_eblgFine.getDBL(), m_refRat);
    LevelData<BaseFab<bool> > cedFineMask(gridsCedFine, ncomp, IntVect::Zero);
    LevelData<BaseFab<bool> > coarMask(   m_gridsCoar,  ncomp, m_redistRad*IntVect::Unit);
    for (DataIterator dit = gridsCedFine.dataIterator(); dit.ok(); ++dit){
      cedFineMask[dit()].setVal(true);
    }
    for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit){
      coarMask[dit()].setVal(false);
    }
    cedFineMask.copyTo(coarMask);
    coarMask.exchange();


    // Make sets
    m_setsCoar.define(m_gridsCoar);
    for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit){
      Box gridBox = m_gridsCoar.get(dit());
      gridBox.grow(m_redistRad);
      gridBox &= m_domainCoar;

      // coarMask holds a local view of the fine level grid, including ghost cells. Make that view into
      // something that IntVectSet can use
      IntVectSet coveringIVS;
      const BaseFab<bool>& mask = coarMask[dit()];
      for (BoxIterator bit(mask.box()); bit.ok(); ++bit){
	const IntVect iv = bit();
	if(mask(iv, comp)) coveringIVS |= iv;
      }
    
      m_setsCoar[dit()]  = m_ebislCoar[dit()].getIrregIVS(gridBox);
      m_setsCoar[dit()] &= coveringIVS;
    }

    // Define data
    defineDataHolders();
    setToZero();

#if EBFASTC2C_DEBUG // Original function
    IntVectSet newDefineSet, oldDefineSet;
    gatherCoarSet(newDefineSet);
    EBCoarToCoarRedist::define(a_eblgFine, a_eblgCoar, a_nref, a_nvar, a_redistRad);
    gatherCoarSet(oldDefineSet);

    const IntVectSet diffSet1 = newDefineSet - oldDefineSet;
    const IntVectSet diffSet2 = oldDefineSet - newDefineSet;

    if(procID() == uniqueProc(SerialTask::compute)) std::cout << diffSet1 << std::endl;
    if(procID() == uniqueProc(SerialTask::compute)) std::cout << diffSet2 << std::endl;
#endif
  }

  void EBFastCoarToCoarRedist::gatherBroadcast(IntVectSet& a_set){
    CH_TIME("EBFastCoarToCoarRedist::gatherBroadcast");
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

  void EBFastCoarToCoarRedist::gatherCoarSet(IntVectSet& a_coarSet){
    CH_TIME("EBFastCoarToCoarRedist::gatherCoarSet");

    a_coarSet.makeEmpty();
    for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit){
      a_coarSet |= m_setsCoar[dit()];
    }
    gatherBroadcast(a_coarSet);
  }
}
