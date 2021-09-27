/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbGhostCloud.cpp
  @brief  Implementation of CD_EbGhostCloud.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <EBCellFactory.H>
#include <EBAverageF_F.H>
#include <EBAlias.H>

// Our includes
#include <CD_EbGhostCloud.H>
#include <CD_EbGhostCloudF_F.H>
#include <CD_NamespaceHeader.H>
  
EbGhostCloud::EbGhostCloud(){
  CH_TIME("EbGhostCloud::EbGhostCloud");
  
  m_isDefined = false;
}

EbGhostCloud::EbGhostCloud(const EBLevelGrid&       a_eblgCoar,
			   const EBLevelGrid&       a_eblgFine,
			   const int                a_refRat,
			   const int                a_nComp,
			   const int                a_ghostFine){
  CH_TIME("EbGhostCloud::EbGhostCloud");
  
  this->define(a_eblgCoar, a_eblgFine, a_refRat, a_nComp, a_ghostFine);
}

EbGhostCloud::~EbGhostCloud(){
  CH_TIME("EbGhostCloud::~EbGhostCloud");
}

void EbGhostCloud::define(const EBLevelGrid&       a_eblgCoar,
			  const EBLevelGrid&       a_eblgFine,
			  const int                a_refRat,
			  const int                a_nComp,
			  const int                a_ghostFine){
  CH_TIME("EbGhostCloud::define");
  
  m_refRat     = a_refRat;
  m_nComp      = a_nComp;
  m_ghost      = a_ghostFine;

  m_eblgCoar   = a_eblgCoar;
  m_eblgFine   = a_eblgFine;

  m_gridsCoar  = m_eblgCoar.getDBL();
  m_gridsFine  = m_eblgFine.getDBL();
  
  m_domainCoar = m_eblgCoar.getDomain();
  m_domainFine = m_eblgFine.getDomain();

  // Make coarsened/refined stuff
  this->makeCoFiStuff();
  this->makeFiCoStuff();

  // Create fine scratch data. The fine data should just be a raw copy. 
  m_scratchFine.define(m_gridsFine, m_nComp, m_ghost*IntVect::Unit, EBCellFactory(m_eblgFine.getEBISL()));

  m_isDefined = true;
}

void EbGhostCloud::makeCoFiStuff(){
  CH_TIME("EbGhostCloud::makeCoFiStuff");
  
  // TLDR: Coarsen the grids and make a BoxLayout data
  
  DisjointBoxLayout dblCoFi;
  coarsen(dblCoFi, m_gridsFine, m_refRat);
  
  Vector<Box> coFiBoxes = dblCoFi.boxArray();
  Vector<int> coFiProcs = dblCoFi.procIDs();

  const int ghostCoFi = (m_ghost + m_refRat - 1)/m_refRat; // Rounds upwards. 
  for (int i = 0; i < coFiBoxes.size(); i++){
    coFiBoxes[i].grow(ghostCoFi);
    coFiBoxes[i] &= m_domainCoar;
  }
  
  m_gridsCoFi.define(coFiBoxes, coFiProcs);
  m_dataCoFi. define(m_gridsCoFi, m_nComp);
}

void EbGhostCloud::makeFiCoStuff(){
  CH_TIME("EbGhostCloud::makeFiCoStuff");
  
  DisjointBoxLayout dblFiCo;
  refine(dblFiCo, m_gridsCoar, m_refRat);

  Vector<Box> fiCoBoxes = dblFiCo.boxArray();
  Vector<int> fiCoProcs = dblFiCo.procIDs();

  const int ghostFiCo = m_ghost;
  for (int i = 0; i < fiCoBoxes.size(); i++){
    fiCoBoxes[i].grow(ghostFiCo);
  }

  // Define the grids and the data holder
  m_gridsFiCo.define(fiCoBoxes, fiCoProcs);
  m_dataFiCo. define(m_gridsFiCo, m_nComp);

  m_eblgFiCo.define(dblFiCo, m_domainFine, ghostFiCo, m_eblgFine.getEBIS());
}

void EbGhostCloud::addFineGhostsToCoarse(LevelData<EBCellFAB>& a_coarData, const LevelData<EBCellFAB>& a_fineData){
  CH_TIME("EbGhostCloud::addFineGhostsToCoarse");
  
  CH_assert(m_isDefined);

  // TLDR: This routine will take the ghost cells that are in a_fineData and add them to the coarse level. We use buffers for this,
  //       copying the data from a_fineData to our nifty m_scratchFine buffer holder. The ghost cells in that scratch data are
  //       coarsened onto yet another buffer (m_dataCoFi). Finally, we add the contents in that buffer to the coarse data. 
  
  // Copy the fine data to scratch and reset the interior cells. We do this by copying everything to the fine scratch data,
  // and then set all the valid data in each box to zero. The exchange() operation will then take care of ghost cells that. 
  // overlap with valid regions in different boxes (we could use a NeighborIterator to achieve the same). 
  a_fineData.localCopyTo(m_scratchFine);
  for (DataIterator dit(m_gridsFine); dit.ok(); ++dit){
    const Box& box = m_gridsFine[dit()];
    m_scratchFine[dit()].setVal(0.0, box, 0, m_nComp);
  }
  m_scratchFine.exchange();


  // Coarsen the fine grid data. We do this by AVERAGING the fine-grid data onto the coarse grid. This actually involves entire patches
  // and not just individual ghost cells regions, but that's ok because the rest of the data (in the valid fine regions) are set to zero above,
  // so the coarse data underneath those regions will be zero, also. 
  const Real factor = 1./pow(m_refRat, SpaceDim);  
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit){
    
    FArrayBox& coFiFab        = m_dataCoFi[dit()];//.getFArrayBox();
    const FArrayBox& fineData = m_scratchFine[dit()].getFArrayBox();

    const Box coarBox         = m_gridsCoFi[dit()];
    const Box fineBox         = fineData.box() & m_domainFine;

    const Box refbox(IntVect::Zero, (m_refRat-1)*IntVect::Unit);

    coFiFab.setVal(0.0);
    
    // Do all the regular cells. This also does irregular cells, but not multi-valued ones.
    // This will not work if you have a refinement boundary that also crosses a multi-valued cell, but I really don't know how to
    // handle deposition in that case. 
    for (int comp = 0; comp < m_nComp; comp++){
      FORT_AVERAGE_GHOSTS(CHF_FRA1(coFiFab, comp),
			  CHF_CONST_FRA1(fineData, comp),
			  CHF_BOX(fineBox),
			  CHF_CONST_INT(m_refRat),
			  CHF_CONST_REAL(factor));
    }
  }

  // Add the data from teh buffer to the coarse data. Again, this will only do single-valued cells!
  const Interval interv(0, m_nComp-1);
  LevelData<FArrayBox> coarAlias;
  aliasEB(coarAlias, a_coarData);
  m_dataCoFi.addTo(interv, coarAlias, interv, m_domainCoar.domainBox());
}

void EbGhostCloud::addFiCoDataToFine(LevelData<EBCellFAB>& a_fineData, const BoxLayoutData<FArrayBox>& a_fiCoData){
  CH_TIME("EbGhostCloud::addFiCoDataToFine");

  // This is really simple because the user will already have used the buffer to deposit particles onto the "fine level". We simply add the data directly here. Again,
  // there might be better ways of doing this because a_fiCoData actually contains a lot of mesh!
  
  const Interval interv(0, m_nComp-1);
  LevelData<FArrayBox> fineAlias;
  aliasEB(fineAlias, a_fineData);
  a_fiCoData.addTo(interv, fineAlias, interv, m_domainFine.domainBox());
}

const BoxLayoutData<FArrayBox>& EbGhostCloud::getFiCoBuffer() const {
  CH_TIME("EbGhostCloud::getFiCoBuffer");
  
  return m_dataFiCo;
}

BoxLayoutData<FArrayBox>& EbGhostCloud::getFiCoBuffer() {
  CH_TIME("EbGhostCloud::getFiCoBuffer");
  
  return m_dataFiCo;
}

const EBLevelGrid& EbGhostCloud::getEblgFiCo() const {
  CH_TIME("EbGhostCloud::getEblgFiCo");
  
  return m_eblgFiCo;
}

#include <CD_NamespaceFooter.H>
