/*!
  @file   ebfinetocoar_redist.cpp
  @brief  Implementation of ebfinetocoar_redist
  @author Robert Marskar
*/

#include "ebfinetocoar_redist.H"

#include <EBArith.H>


ebfinetocoar_redist::ebfinetocoar_redist() : EBFineToCoarRedist(){

}

ebfinetocoar_redist::~ebfinetocoar_redist(){

}

void ebfinetocoar_redist::new_define(const EBLevelGrid&  a_eblgFine,
				     const EBLevelGrid&  a_eblgCoar,
				     const int&          a_nref,
				     const int&          a_nvar,
				     const int&          a_redistRad){

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
#if 1 // Debug code
  IntVectSet fineShell2;
#endif
  
  {
    CH_TIME("make_fine_sets");
#if 0 // Old way of building the shell
    IntVectSet fineInterior = a_eblgFine.getCoveringIVS();
    EBArith::shrinkIVS(fineInterior, m_redistRad);
    fineShell = a_eblgFine.getCoveringIVS();
    fineShell -= fineInterior;
#else // New way of doing this


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

      // Loop through the CFIVS and grow point into a cube in each direction. The length of the cube
      // in each direction is 1+2*m_redistRad. 
      IntVectSet grown_cfivs;
      for (IVSIterator ivsIt(cfivs); ivsIt.ok(); ++ivsIt){
	const IntVect iv = ivsIt();

	Box grownBox(iv, iv);
	grownBox.grow(m_redistRad);

	grown_cfivs |= IntVectSet(grownBox);
      }
      

      // Intersect the cube-grown CFIVS with the grid. This will discard everything on the outside of the CFIVS but
      // we keep everything on the inside. Importantly, we also grew the corner cells. 
      for (LayoutIterator lit = m_gridsFine.layoutIterator(); lit.ok(); ++lit){
	const Box fineBox = m_gridsFine.get(dit());
	grown_cfivs &= fineBox;
      }

      localShell |= grown_cfivs;
    }

    // Gather-broadcast the local shells
    Vector<IntVectSet> all_shells;
    const int dest_proc = uniqueProc(SerialTask::compute);
    gather(all_shells, localShell, dest_proc);
    if(procID() == dest_proc){
      for (int i = 0; i < all_shells.size(); i++){
	fineShell2 |= all_shells[i];
      }

      fineShell2.compact();
    }

    broadcast(fineShell2, dest_proc);

#if 0 // Debug
    if(procID() == 0){
      const int numShell  = fineShell.numPts();
      const int numShell2 = fineShell2.numPts();
      std::cout << "fineShell.numPts() = << " << numShell << "\t fineShell2.numPts() = " << numShell2 << std::endl;
    }

    IntVectSet diffShell  = fineShell;
    diffShell -= fineShell2;

    IntVectSet diffShell2 = fineShell2;
    diffShell2 -= fineShell;
    if(diffShell.numPts() != 0 || diffShell2.numPts() != 0){
      MayDay::Abort("stop");
    }
#endif
#endif
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

  // Define data as usual
  defineDataHolders();
  setToZero();
}
