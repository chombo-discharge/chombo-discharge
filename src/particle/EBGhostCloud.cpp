/*!
  @file   EBGhostCloud.cpp
  @brief  Implementation of EBGhostCloud.H
  @author Robert Marskar
  @date   April 2020
*/

#include "EBGhostCloud.H"
#include "EBGhostCloudF_F.H"
#include "EBCellFactory.H"
#include "EBAverageF_F.H"
#include "EBAlias.H"

EBGhostCloud::EBGhostCloud(){

}

EBGhostCloud::EBGhostCloud(const DisjointBoxLayout& a_gridsCoar,
			   const DisjointBoxLayout& a_gridsFine,
			   const EBLevelGrid&       a_eblgCoar,
			   const EBLevelGrid&       a_eblgFine,
			   const ProblemDomain&     a_domainCoar,
			   const ProblemDomain&     a_domainFine,
			   const int                a_refRat,
			   const int                a_nComp,
			   const int                a_ghostFine){

  this->define(a_gridsCoar, a_gridsFine, a_eblgCoar, a_eblgFine, a_domainCoar, a_domainFine, a_refRat, a_nComp, a_ghostFine);
}

EBGhostCloud::~EBGhostCloud(){

}

void EBGhostCloud::define(const DisjointBoxLayout& a_gridsCoar,
			  const DisjointBoxLayout& a_gridsFine,
			  const EBLevelGrid&       a_eblgCoar,
			  const EBLevelGrid&       a_eblgFine,
			  const ProblemDomain&     a_domainCoar,
			  const ProblemDomain&     a_domainFine,
			  const int                a_refRat,
			  const int                a_nComp,
			  const int                a_ghostFine){

  m_refRat     = a_refRat;
  m_nComp      = a_nComp;
  m_ghost      = a_ghostFine;
  m_gridsFine  = a_gridsFine;
  m_domainCoar = a_domainCoar;

  // Coarsen the grids and make a BoxLayout grid
  DisjointBoxLayout dblCoFi;
  
  coarsen(dblCoFi,    a_gridsFine, m_refRat);
  coarsen(m_eblgCoFi, a_eblgFine,  m_refRat);
  
  Vector<Box> coFiBoxes = dblCoFi.boxArray();
  Vector<int> coFiProcs = dblCoFi.procIDs();

  const int ghostCoFi = (a_ghostFine + m_refRat - 1)/m_refRat;
  for (int i = 0; i < coFiBoxes.size(); i++){
    coFiBoxes[i].grow(ghostCoFi);
    coFiBoxes[i] &= a_domainCoar;
  }
  
  // Create scratch data. The fine data should just be a raw copy. The coarse data is the coarse data. 
  const EBCellFactory factFine(a_eblgFine.getEBISL());
  const EBCellFactory factCoFi(m_eblgCoFi.getEBISL());

  m_scratchFine.define(a_gridsFine, a_nComp, m_ghost*IntVect::Unit, factFine);
  
  m_gridsCoFi.define(coFiBoxes, coFiProcs);
  m_dataCoFi.define(m_gridsCoFi, a_nComp);

  m_isDefined = true;
}

void EBGhostCloud::addFineGhostsToCoarse(LevelData<EBCellFAB>& a_coarData, const LevelData<EBCellFAB>& a_fineData){
  CH_assert(m_isDefined);
  
  // Copy the fine data to scratch and reset interior cells
  a_fineData.localCopyTo(m_scratchFine);
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
    const Box& box = m_gridsFine.get(dit());
    m_scratchFine[dit()].setVal(0.0, box, 0, m_nComp);
  }
  m_scratchFine.exchange();

  const Real factor = 1./pow(m_refRat, SpaceDim);

  // Coarsen the fine grid data
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
    
    FArrayBox& coFiFab        = m_dataCoFi[dit()];//.getFArrayBox();
    const FArrayBox& fineData = m_scratchFine[dit()].getFArrayBox();

    const Box coarBox         = m_gridsCoFi[dit()];
    const Box fineBox         = fineData.box();

    const Box refbox(IntVect::Zero, (m_refRat-1)*IntVect::Unit);

    coFiFab.setVal(0.0);

    // Do all the regular cells. This also does irregular cells, but not multi-valued ones.
    // If you have a refinement boundary that also crosses a multi-valued cell, feel free to fuck off. 
    for (int comp = 0; comp < m_nComp; comp++){
      FORT_AVERAGE_GHOSTS(CHF_FRA1(coFiFab, comp),
			  CHF_CONST_FRA1(fineData, comp),
			  CHF_BOX(fineBox),
			  CHF_CONST_INT(m_refRat),
			  CHF_CONST_REAL(factor));
    }
  }

  // Add the data
  const Interval interv(0, m_nComp-1);
  LevelData<FArrayBox> coarAlias;
  aliasEB(coarAlias, a_coarData);
  m_dataCoFi.addTo(interv, coarAlias, interv, m_domainCoar.domainBox());
}
