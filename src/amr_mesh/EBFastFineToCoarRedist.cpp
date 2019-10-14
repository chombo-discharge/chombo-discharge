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

void EBFastFineToCoarRedist::new_define(const EBLevelGrid&                      a_eblgFine,
				     const EBLevelGrid&                      a_eblgCoar,
				     const LayoutData<Vector<LayoutIndex> >& a_neighborsFine,
				     const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar,
				     const int&                              a_nref,
				     const int&                              a_nvar,
				     const int&                              a_redistRad){
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
  //global set consists of redistrad on the fine side of the coarse-fine
  //interface
  //the fine set is that within one fine box.
  //the refcoar set is that set within one refined coarse box.
  IntVectSet fineShell;    
  
  {
    CH_TIME("make_fine_sets");

    // Get the strip of cells immediately on the outside of this DBL. No restriction to ProblemDomain here! 
    LayoutData<IntVectSet> outsideCFIVS;
    outsideCFIVS.define(m_gridsFine);
    for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
      Box grownBox = grow(m_gridsFine.get(dit()), 1);
      outsideCFIVS[dit()] = IntVectSet(grownBox);
      for (LayoutIterator lit = m_gridsFine.layoutIterator(); lit.ok(); ++lit){
	outsideCFIVS[dit()] -= m_gridsFine[lit()];
      }
    }

    IntVectSet localShell;
    for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
      const IntVectSet& cfivs = outsideCFIVS[dit()];

      // Loop through the CFIVS and grow each point into a cube. The length of the cube
      // in each direction is 1+2*m_redistRad. 
      IntVectSet grown_cfivs;
      for (IVSIterator ivsIt(cfivs); ivsIt.ok(); ++ivsIt){
	const IntVect iv = ivsIt();
	
	Box grownBox(iv, iv);
	grownBox.grow(m_redistRad);
	
	grown_cfivs |= IntVectSet(grownBox);
      }
      

      // Intersect the cube-grown CFIVS with the grid. This will discard everything on the outside of the CFIVS but
      // we keep everything on the inside. Importantly, we also grew the corner cells so there shouldn't be any issues
      // with internal corners on the fine shell. 
      for (LayoutIterator lit = m_gridsFine.layoutIterator(); lit.ok(); ++lit){
	const Box fineBox = m_gridsFine.get(dit());
	grown_cfivs &= fineBox;
      }

      localShell |= grown_cfivs;
    }

    // Gather-broadcast the local shells to a global shell. 
    Vector<IntVectSet> all_shells;
    const int dest_proc = uniqueProc(SerialTask::compute);
    gather(all_shells, localShell, dest_proc);
    if(procID() == dest_proc){
      for (int i = 0; i < all_shells.size(); i++){
	fineShell |= all_shells[i];
      }
      fineShell.compact();
    }
    broadcast(fineShell, dest_proc);

    // Define fine sets
    for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
      const Box& fineBox = m_gridsFine.get(dit());
      m_setsFine[dit()]  = m_ebislFine[dit()].getIrregIVS(fineBox);
      m_setsFine[dit()] &= fineShell;
    }
  }
  {
    CH_TIME("make_coar_sets");
    for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit)
      {
        Box grownBox = grow(m_gridsRefCoar.get(dit()), m_redistRad);
        grownBox &= domainFine;
        m_setsRefCoar[dit()] = m_ebislRefCoar[dit()].getIrregIVS(grownBox);
        m_setsRefCoar[dit()] &= fineShell;
      }
  }

  // Call old stuff
  //  EBFineToCoarRedist::define(a_eblgFine, a_eblgCoar, a_nref, a_nvar, a_redistRad);

#if 0 // Debugging hook for coarse set. Make sure we have the correct amount of cells
  IntVectSet localSetsCoar, globalSetsCoar;

  // Gather fine sets
  for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit){
    localSetsCoar |= m_setsRefCoar[dit()];
  }
  const int dest_proc = uniqueProc(SerialTask::compute);
  Vector<IntVectSet> procSetsCoar;
  gather(procSetsCoar, localSetsCoar, dest_proc);
  if(procID() == dest_proc){
    for (int i = 0; i < procSetsCoar.size(); i++){
      globalSetsCoar |= procSetsCoar[i];
      }
    globalSetsCoar.compact();
  }
  broadcast(globalSetsCoar, dest_proc);

  if(procID() == 0){
    std::cout << "old define: \n" << globalSetsCoar << std::endl;
  }
#endif

#if 0 // Debugging hook for fine set. Make sure we have the correct amount of cells
  IntVectSet localSetsFine, globalSetsFine;

  // Gather fine sets
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
    localSetsFine |= m_setsFine[dit()];
  }
  const int dest_proc = uniqueProc(SerialTask::compute);
  Vector<IntVectSet> procSetsFine;
  gather(procSetsFine, localSetsFine, dest_proc);
  if(procID() == dest_proc){
    for (int i = 0; i < procSetsFine.size(); i++){
      globalSetsFine |= procSetsFine[i];
      }
    globalSetsFine.compact();
  }
  broadcast(globalSetsFine, dest_proc);

  if(procID() == 0){
    std::cout << "old define: \n" << globalSetsFine << std::endl;
  }
#endif

  

  // Define data as usual
  defineDataHolders();
  setToZero();
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
  //global set consists of redistrad on the fine side of the coarse-fine interface interface, but is not computed
  //the fine set is that within one fine box, using local information about the CFIVS 
  //the refcoar set is that set same within one refined coarse box. This needs to be computed using a slower computation
  m_fineShell    = new LevelData<BaseFab<bool> >(m_gridsFine,    1, IntVect::Zero);
  m_refCoarShell = new LevelData<BaseFab<bool> >(m_gridsRefCoar, 1, IntVect::Zero);
  
  makeEmptyFineShell();    // Init bits to zero
  makeEmptyRefCoarShell(); // Init bits to zero
  
  makeFineSets(a_neighborsFine);
  makeCoarSets(a_neighborsCoar); // a_neighborsCoar is not used. 

  delete m_fineShell;
  delete m_refCoarShell;

  // Define data as usual
  defineDataHolders();
  setToZero();


#if 0 // Debugging hook. Make sure we have the correct amount of cells
  IntVectSet localSetsCoar, globalSetsCoar;

  // Gather fine sets
  for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit){
    localSetsCoar |= m_setsRefCoar[dit()];
  }
  const int dest_proc = uniqueProc(SerialTask::compute);
  Vector<IntVectSet> procSetsCoar;
  gather(procSetsCoar, localSetsCoar, dest_proc);
  if(procID() == dest_proc){
    for (int i = 0; i < procSetsCoar.size(); i++){
      globalSetsCoar |= procSetsCoar[i];
      }
    globalSetsCoar.compact();
  }
  broadcast(globalSetsCoar, dest_proc);

  if(procID() == 0){
    std::cout << "new define = \n" << globalSetsCoar << std::endl;
  }
#endif

#if 0 // Debugging hook for fine set. Make sure we have the correct amount of cells
  IntVectSet localSetsFine, globalSetsFine;

  // Gather fine sets
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
    localSetsFine |= m_setsFine[dit()];
  }
  const int dest_proc = uniqueProc(SerialTask::compute);
  Vector<IntVectSet> procSetsFine;
  gather(procSetsFine, localSetsFine, dest_proc);
  if(procID() == dest_proc){
    for (int i = 0; i < procSetsFine.size(); i++){
      globalSetsFine |= procSetsFine[i];
      }
    globalSetsFine.compact();
  }
  broadcast(globalSetsFine, dest_proc);

  if(procID() == 0){
    std::cout << "new define: \n" << globalSetsFine << std::endl;
  }
#endif
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

    // Gather shell shell here
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
