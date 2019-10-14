/*!
  @file   EBFastFineToCoarRedist.cpp
  @brief  Implementation of EBFastFineToCoarRedist
  @author Robert Marskar
*/

#include "EBFastFineToCoarRedist.H"

#include <EBArith.H>


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
  ebisPtr->fillEBISLayout(m_ebislRefCoar, m_gridsRefCoar,
                          domainFine, nghost);
  m_ebislRefCoar.setMaxCoarseningRatio(m_refRat,ebisPtr);

  //define the intvectsets over which the objects live
  m_setsFine.define(m_gridsFine);
  m_setsRefCoar.define(m_gridsCoar);

  IntVectSet globalShell;
  globalShell.makeEmpty();

  makeFineSets(globalShell, a_neighborsFine); // a_neighborsFine is used for rapidly building the inside CFIVS 
  makeCoarSets(globalShell, a_neighborsCoar); // a_neighborsCoar is not used. 

  // Define data as usual
  defineDataHolders();
  setToZero();
}

void EBFastFineToCoarRedist::makeFineSets(IntVectSet& a_globalShell, const LayoutData<Vector<LayoutIndex> >& a_neighborsFine){
  CH_TIME("EBFastFineToCoarRedist::makeFineSets");

  IntVectSet localShell;
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
    const Box& fineBox = m_gridsFine.get(dit());
    
    // Make a local view of the CFIVS
    Box grownBox = grow(fineBox, 1);
    IntVectSet localCFIVS = IntVectSet(grownBox);
    localCFIVS -= fineBox;
    const Vector<LayoutIndex>& neighbors = a_neighborsFine[dit()];
    for (int i = 0; i < neighbors.size(); i++){
      localCFIVS -= m_gridsFine[neighbors[i]];
    }
    
    // Make the local shell by growing the CFISV and restricting to box
    IntVectSet curShell;
    for (IVSIterator ivsIt(localCFIVS); ivsIt.ok(); ++ivsIt){
      const IntVect iv = ivsIt();
      Box b(iv,iv);
      b.grow(m_redistRad);
      curShell |= IntVectSet(b);
    }
    curShell &= fineBox;
    localShell |= curShell;
    
    // Define the fine set
    m_setsFine[dit()]  = m_ebislFine[dit()].getIrregIVS(fineBox);
    m_setsFine[dit()] &= localShell;
  }

#ifdef CH_MPI  // Build the global view of the shell
  Vector<IntVectSet> procShells;
  const int destProc = uniqueProc(SerialTask::compute);
  gather(procShells, localShell, destProc);
  if(procID() == destProc){
    for (int i = 0; i < procShells.size(); i++){
      a_globalShell |= procShells[i];
    }
  }
  broadcast(a_globalShell, destProc);
#endif
}

void EBFastFineToCoarRedist::makeCoarSets(const IntVectSet& a_globalShell,
					  const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar){
  CH_TIME("EBFastFineToCoarRedist::makeCoarSets");

  const Box domainFine = refine(m_domainCoar, m_refRat);
  for (DataIterator dit = m_gridsRefCoar.dataIterator(); dit.ok(); ++dit){
    const Box& refCoarBox = m_gridsRefCoar.get(dit());

    // Make the coar set
    Box grownBox = grow(refCoarBox, m_redistRad);
    grownBox &= domainFine;
    m_setsRefCoar[dit()]  = m_ebislRefCoar[dit()].getIrregIVS(grownBox);
    m_setsRefCoar[dit()] &= a_globalShell;
  }
}
