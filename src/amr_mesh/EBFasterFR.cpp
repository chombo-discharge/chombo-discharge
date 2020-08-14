/*!
  @brief   EBFasterFR.cpp
  @details Implementation of EBFasterFR.H
  @author  Robert Maskar
  @date    Aug. 2020
*/

#include "EBFasterFR.H"

EBFasterFR::EBFasterFR(){

}

EBFasterFR::EBFasterFR(const EBLevelGrid&       a_eblgFine,
		       const EBLevelGrid&       a_eblgCoar,
		       const int&               a_nref,
		       const int&               a_nvar,
		       const bool               a_forceNoEBCF){
  setDefaultValues();
  define(a_eblgFine, a_eblgCoar, a_nref, a_nvar, a_forceNoEBCF);
}

EBFasterFR::~EBFasterFR(){

}

void EBFasterFR::define(const EBLevelGrid&       a_eblgFine,
			const EBLevelGrid&       a_eblgCoar,
			const int&               a_refRat,
			const int&               a_nvar,
			const bool               a_forceNoEBCF){
  CH_TIME("EBFasterFR::define");
  
  clear();
  m_refRat   = a_refRat;
  m_nComp    = a_nvar;
  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;
  if (!m_eblgFine.coarsenable(a_refRat))
    {
      MayDay::Error("EBFastFR::define: dbl not coarsenable by refrat");
    }

  coarsen(m_eblgCoFi, a_eblgFine, a_refRat);
  m_eblgCoFi.getEBISL().setMaxRefinementRatio(a_refRat, m_eblgCoFi.getEBIS());
  m_reverseCopier.ghostDefine(m_eblgCoFi.getDBL(), m_eblgCoar.getDBL(),m_eblgCoar.getDomain(),IntVect::Unit);

#if (CH_SPACEDIM == 2)
  m_nrefdmo = m_refRat;
#elif (CH_SPACEDIM == 3)
  m_nrefdmo = m_refRat*m_refRat;
#else
  bogus spacedim;
#endif

  m_levelFluxReg = new LevelFluxRegister();

  //  pout() << "before regular define" << endl;
  m_levelFluxReg->define(a_eblgFine.getDBL(),
                         a_eblgCoar.getDBL(),
                         a_eblgFine.getDomain(),
                         a_refRat, a_nvar);
  if (a_forceNoEBCF)
    {
      m_hasEBCF = false;
    }
  else
    {
      m_hasEBCF = computeHasEBCF();
    }

  //if no EBCF--nothing happens here but calls to level flux register
  if (m_hasEBCF)
    {
      fastDefineSetsAndIterators();
      defineBuffers();
    }
  //set all registers to zero to start.
  setToZero();
  m_isDefined = true;
  //  pout() << "leaving define" << endl;
}

void EBFasterFR::fastDefineSetsAndIterators(){
  CH_TIME("EBFasterFR::fastDefineSetsAndIterators");
  
  for (int idir = 0; idir < SpaceDim; idir++){
    for (SideIterator sit; sit.ok(); ++sit){
      int iindex = index(idir, sit());
      m_setsCoFi[iindex].define(m_eblgCoFi.getDBL());
      m_vofiCoFi[iindex].define(m_eblgCoFi.getDBL());

      for (DataIterator dit = m_eblgCoFi.getDBL().dataIterator(); dit.ok(); ++dit){
	const Box& grid = m_eblgCoFi.getDBL().get(dit());
	Box boxCoFi;
	if (sit() == Side::Lo){
	  boxCoFi = adjCellLo(grid, idir, 1);}
	else{
	  boxCoFi = adjCellHi(grid, idir, 1);
	}
	const EBGraph& ebgraphCoFi =    m_eblgCoFi.getEBISL()[dit()].getEBGraph();
	const IntVectSet& ivsCF    =  (*m_eblgCoFi.getCFIVS())[dit()];
	IntVectSet ivsEBCF         = ebgraphCoFi.getIrregCells(boxCoFi);

	ivsEBCF &= ivsCF;
	ivsEBCF &= boxCoFi;

	(m_setsCoFi[iindex])[dit()]      = ivsEBCF;
	(m_vofiCoFi[iindex])[dit()].define(ivsEBCF, ebgraphCoFi);
      }

      m_setsCoar[iindex].define(m_eblgCoar.getDBL());
      m_vofiCoar[iindex].define(m_eblgCoar.getDBL());

      for (DataIterator dit = m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit){
	const Box& boxCoar = m_eblgCoar.getDBL()[dit()];
	IntVectSet irregIVS = m_eblgCoar.getEBISL()[dit()].getIrregIVS(boxCoar);
	Vector<IntVectSet>  cfIVSVec;
	int iindex = index(idir, sit());
	getCFIVSSubset(cfIVSVec, boxCoar, m_eblgCoFi.getDBL(), idir, sit());
	const EBGraph& ebgraphCoar =  m_eblgCoar.getEBISL()[dit()].getEBGraph();
	(m_setsCoar[iindex])[dit()].resize(cfIVSVec.size());
	(m_vofiCoar[iindex])[dit()].resize(cfIVSVec.size());
	for (int ivec = 0; ivec < cfIVSVec.size(); ivec++){
	  IntVectSet ivsEBCF = irregIVS;
	  ivsEBCF &= cfIVSVec[ivec];

	  //need to grow by one to get on both sides of
	  //coarse-fine interface
	  (m_setsCoar[iindex])[dit()][ivec]      = ivsEBCF;
	  (m_vofiCoar[iindex])[dit()][ivec].define(ivsEBCF, ebgraphCoar);
	}
      }
    }
  }
}

void EBFasterFR::getCFIVSSubset(Vector<IntVectSet>&      a_cfivsCoFi,
				const Box&               a_subreCoar,
				const DisjointBoxLayout& a_layouCoFi,
				const int&               a_idir,
				const Side::LoHiSide&    a_sd){

  a_cfivsCoFi.resize(0);
  for (LayoutIterator lito = a_layouCoFi.layoutIterator(); lito.ok(); ++lito){
    const Box& gridCoFiO = a_layouCoFi[lito()];
    //the - 1 is due to the adjacent cell box thing
    if (gridCoFiO.bigEnd(0) < a_subreCoar.smallEnd(0) - 1 ){
      //can skip rest cuz we haven't gotten to something interesting
      continue;
    }

    Box cfboxCoFi = adjCellBox(gridCoFiO, a_idir, a_sd, 1);

    //only want the cells within the input subregion
    cfboxCoFi &= a_subreCoar;

    IntVectSet ivsCoFi(cfboxCoFi);
    for (LayoutIterator liti = a_layouCoFi.layoutIterator(); liti.ok(); ++liti){
      const Box& gridCoFiI = a_layouCoFi[liti()];
      if (gridCoFiI.bigEnd(0) < cfboxCoFi.smallEnd(0)){
	//can skip rest cuz we haven't gotten to something interesting
	continue;
      }

      //remove from the set boxes that overlap on the own grid
      //(fine-fine interface points)
      ivsCoFi -= gridCoFiI;

      if ((gridCoFiI.smallEnd(0) > cfboxCoFi.bigEnd(0)) || ivsCoFi.isEmpty()){
	//can break out of loop, since we know that the smallEnd
	// of all the remaining boxes are lexigraphically beyond this ghosted box.
	break;
      }
    }
    if (!ivsCoFi.isEmpty()){
      a_cfivsCoFi.push_back(ivsCoFi);
    }

    //the + 1 is due to the adjacent cell box thing
    if (gridCoFiO.smallEnd(0) > a_subreCoar.bigEnd(0) + 1 ){
      //can break out of loop, since we know that the smallEnd
      // of all the remaining boxes are lexigraphically beyond this ghosted box.
      break;
    }
  }
}
