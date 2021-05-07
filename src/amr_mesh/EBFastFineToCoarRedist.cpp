/*!
  @file   EBFastFineToCoarRedist.cpp
  @brief  Implementation of EBFastFineToCoarRedist
  @author Robert Marskar
*/

#include "EBFastFineToCoarRedist.H"

#include <EBArith.H>

#define EBFASTF2C_DEBUG 0

namespace ChomboDischarge {

  EBFastFineToCoarRedist::EBFastFineToCoarRedist() : EBFineToCoarRedist(){

  }

  EBFastFineToCoarRedist::~EBFastFineToCoarRedist(){

  }

  void EBFastFineToCoarRedist::define(const EBLevelGrid&                      a_eblgFine,
				      const EBLevelGrid&                      a_eblgCoar,
				      const LayoutData<Vector<LayoutIndex> >& a_neighborsFine,
				      const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar,
				      const int&                              a_nref,
				      const int&                              a_nvar,
				      const int&                              a_redistRad){
    CH_TIME("EBFastFineToCoarRedist::define");

    const int comp  = 0;
    const int ncomp = 1;
    m_isDefined = true;
    m_nComp     = a_nvar;
    m_refRat    = a_nref;
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

    makeFineSets(fineShellMask);
    makeCoarSets(refCoarShellMask);

    defineDataHolders();
    setToZero();

#if EBFASTF2C_DEBUG // This debugging hook calls the original function and checks that the sets completely overlap. 
    IntVectSet newFineSet,    oldFineSet;
    IntVectSet newRefCoarSet, oldRefCoarSet;
  
    gatherSetsFine(newFineSet);
    gatherSetsRefCoar(newRefCoarSet);

    // Call old define
    EBFineToCoarRedist::define(a_eblgFine, a_eblgCoar, a_nref, a_nvar, a_redistRad);
    gatherSetsFine(oldFineSet);
    gatherSetsRefCoar(oldRefCoarSet);

    const IntVectSet diffSet1 = newFineSet    - oldFineSet;
    const IntVectSet diffSet2 = oldFineSet    - newFineSet;
    const IntVectSet diffSet3 = newRefCoarSet - oldRefCoarSet;
    const IntVectSet diffSet4 = oldRefCoarSet - newRefCoarSet;

    if(diffSet1.numPts() != 0) MayDay::Abort("EBFastFineToCoarRedist::define - diffSet1 not empty");
    if(diffSet2.numPts() != 0) MayDay::Abort("EBFastFineToCoarRedist::define - diffSet2 not empty");
    if(diffSet3.numPts() != 0) MayDay::Abort("EBFastFineToCoarRedist::define - diffSet3 not empty");
    if(diffSet4.numPts() != 0) MayDay::Abort("EBFastFineToCoarRedist::define - diffSet4 not empty");
#endif
  }

  void EBFastFineToCoarRedist::makeFineMask(LevelData<BaseFab<bool> >&              a_fineShellMask,
					    const LayoutData<Vector<LayoutIndex> >& a_neighborsFine){
    CH_TIME("EBFastFineToCoarRedist::makeFineMask");

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

  void EBFastFineToCoarRedist::makeFineSets(const LevelData<BaseFab<bool> >& a_fineMask){
    CH_TIME("EBFastFineToCoarRedist::makeFineSets");

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

  void EBFastFineToCoarRedist::makeCoarSets(const LevelData<BaseFab<bool> >& a_refCoarMask){
    CH_TIME("EBFastFineToCoarRedist::makeCoarSets");

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

  void EBFastFineToCoarRedist::gatherBroadcast(IntVectSet& a_set){
    CH_TIME("EBFastFineToCoarRedist::gatherBroadcast");
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

  void EBFastFineToCoarRedist::gatherSetsFine(IntVectSet& a_setsFine){
    CH_TIME("EBFastFineToCoarRedist::gatherSetsFine");

    a_setsFine.makeEmpty();
    for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
      a_setsFine |= m_setsFine[dit()];
    }
    gatherBroadcast(a_setsFine);
  }

  void EBFastFineToCoarRedist::gatherSetsRefCoar(IntVectSet& a_setsRefCoar){
    CH_TIME("EBFastFineToCoarRedist::gatherSetsRefCoar");

    a_setsRefCoar.makeEmpty();
    for (DataIterator dit = m_gridsRefCoar.dataIterator(); dit.ok(); ++dit){
      a_setsRefCoar |= m_setsRefCoar[dit()];
    }
    gatherBroadcast(a_setsRefCoar);
  }
}
