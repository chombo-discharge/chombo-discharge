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

  //make sets
  m_fineShell    = new LevelData<BaseFab<bool> >(m_gridsFine,    1, IntVect::Zero);
  m_refCoarShell = new LevelData<BaseFab<bool> >(m_gridsRefCoar, 1, IntVect::Zero);
  
  makeEmptyFineShell();    // Init mask to zero
  makeEmptyRefCoarShell(); // Init mask to zero
  
  makeFineSets(a_neighborsFine); // a_neighborsFine is used for rapidly building the inside CFIVS 
  makeCoarSets(a_neighborsCoar); // a_neighborsCoar is not used. 

  delete m_fineShell;
  delete m_refCoarShell;

  // Define data as usual
  defineDataHolders();
  setToZero();
}

void EBFastFineToCoarRedist::makeEmptyFineShell(){
  CH_TIME("EBFastFineToCoarRedist::makeEmptyFineShell");
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
    (*m_fineShell)[dit()].setVal(false);
  }
}

void EBFastFineToCoarRedist::makeEmptyRefCoarShell(){
  CH_TIME("EBFastFineToCoarRedist::makeEmptyCoarShell");
  for (DataIterator dit = m_gridsRefCoar.dataIterator(); dit.ok(); ++dit){
    (*m_refCoarShell)[dit()].setVal(false);
  }
}

void EBFastFineToCoarRedist::makeFineSets(const LayoutData<Vector<LayoutIndex> >& a_neighborsFine){
  CH_TIME("EBFastFineToCoarRedist::makeFineSets");

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
    IntVectSet localShell;
    for (IVSIterator ivsIt(localCFIVS); ivsIt.ok(); ++ivsIt){
      const IntVect iv = ivsIt();
      Box b(iv,iv);
      b.grow(m_redistRad);
      localShell |= IntVectSet(b);
    }
    localShell &= fineBox;

    // Cache the fineShell
    for (IVSIterator ivsIt(localShell); ivsIt.ok(); ++ivsIt){
      (*m_fineShell)[dit()](ivsIt(), 0) = true;
    }
    
    // Define the fine set
    m_setsFine[dit()]  = m_ebislFine[dit()].getIrregIVS(fineBox);
    m_setsFine[dit()] &= localShell;
  }
}

void EBFastFineToCoarRedist::makeCoarSets(const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar){
  CH_TIME("EBFastFineToCoarRedist::makeCoarSets");

  const Box domainFine = refine(m_domainCoar, m_refRat);

  // Copy mask onto coarsened fine grids
  m_fineShell->copyTo(*m_refCoarShell);
  
  for (DataIterator dit = m_gridsRefCoar.dataIterator(); dit.ok(); ++dit){
    const Box& refCoarBox = m_gridsRefCoar.get(dit());

    // Gather the shell here
    IntVectSet localShell;
    for (BoxIterator bit(refCoarBox); bit.ok(); ++bit){
      const IntVect iv = bit();
      if((*m_refCoarShell)[dit()](iv,0)) localShell |= iv;
    }

    // Make the coar set
    Box grownBox = grow(m_gridsRefCoar.get(dit()), m_redistRad);
    grownBox &= domainFine;
    m_setsRefCoar[dit()] = m_ebislRefCoar[dit()].getIrregIVS(grownBox);
    m_setsRefCoar[dit()] &= localShell;
  }
}
