/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBCoarseFineParticleMesh.cpp
  @brief  Implementation of CD_EBCoarseFineParticleMesh.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <EBCellFactory.H>
#include <EBAverageF_F.H>
#include <EBAlias.H>

// Our includes
#include <CD_EBCoarseFineParticleMesh.H>
#include <CD_EBCoarseFineParticleMeshF_F.H>
#include <CD_EBAddOp.H>
#include <CD_NamespaceHeader.H>
  
EBCoarseFineParticleMesh::EBCoarseFineParticleMesh(){
  CH_TIME("EBCoarseFineParticleMesh::EBCoarseFineParticleMesh");
  
  m_isDefined = false;
}

EBCoarseFineParticleMesh::EBCoarseFineParticleMesh(const EBLevelGrid& a_eblgCoar,
						   const EBLevelGrid& a_eblgFine,
						   const int          a_refRat,
						   const int          a_nComp,
						   const int          a_ghostFine){
  CH_TIME("EBCoarseFineParticleMesh::EBCoarseFineParticleMesh");
  
  this->define(a_eblgCoar, a_eblgFine, a_refRat, a_nComp, a_ghostFine);
}

EBCoarseFineParticleMesh::~EBCoarseFineParticleMesh(){
  CH_TIME("EBCoarseFineParticleMesh::~EBCoarseFineParticleMesh");
}

void EBCoarseFineParticleMesh::define(const EBLevelGrid& a_eblgCoar,
				      const EBLevelGrid& a_eblgFine,
				      const int          a_refRat,
				      const int          a_nComp,
				      const int          a_ghostFine) {
  CH_TIME("EBCoarseFineParticleMesh::define");
  
  m_refRat     = a_refRat;
  m_nComp      = a_nComp;
  m_ghost      = a_ghostFine;

  m_eblgCoar   = a_eblgCoar;
  m_eblgFine   = a_eblgFine;

  // Make the coarsened and refine grids.
  coarsen(m_eblgCoFi, m_eblgFine, m_refRat);  
  refine (m_eblgFiCo, m_eblgCoar, m_refRat);
  
  // Create various buffers. The fine data should just be a raw copy.
  m_bufferCoFi.define(m_eblgCoFi.getDBL(), m_nComp, m_ghost*IntVect::Unit, EBCellFactory(m_eblgCoFi.getEBISL()));
  m_bufferFiCo.define(m_eblgFiCo.getDBL(), m_nComp, m_ghost*IntVect::Unit, EBCellFactory(m_eblgFiCo.getEBISL()));  
  m_bufferFine.define(m_eblgFine.getDBL(), m_nComp, m_ghost*IntVect::Unit, EBCellFactory(m_eblgFine.getEBISL()));

  // Define Copiers
  m_copierFiCoToFine.ghostDefine(m_eblgFiCo.getDBL(), m_eblgFine.getDBL(), m_eblgFine.getDomain(), m_ghost*IntVect::Unit);
  m_copierCoFiToCoar.ghostDefine(m_eblgCoFi.getDBL(), m_eblgCoar.getDBL(), m_eblgCoar.getDomain(), m_ghost*IntVect::Unit);

  // Define VoF iterators
  this->defineVoFIterators();

  m_isDefined = true;
}

void EBCoarseFineParticleMesh::defineVoFIterators(){
  CH_TIME("EBCoarseFineParticleMesh::defineVoFIterators");

  const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL   ();
  const ProblemDomain&     domainFine = m_eblgFine.getDomain();
  const EBISLayout&        ebislFine  = m_eblgFine.getEBISL ();

  const DisjointBoxLayout& dblCoFi    = m_eblgCoFi.getDBL   ();
  const ProblemDomain&     domainCoFi = m_eblgCoFi.getDomain();
  const EBISLayout&        ebislCoFi  = m_eblgCoFi.getEBISL ();  

  m_vofIterFineGhosts.define(dblFine);
  m_vofIterCoFiGhosts.define(dblCoFi);  
  
  for (DataIterator dit(dblFine); dit.ok(); ++dit){
    const Box&     cellBoxFine = dblFine  [dit()];    
    const EBISBox& ebisBoxFine = ebislFine[dit()];
    const EBGraph& ebgraphFine = ebisBoxFine.getEBGraph();

    const Box&     cellBoxCoFi = dblCoFi  [dit()];    
    const EBISBox& ebisBoxCoFi = ebislCoFi[dit()];
    const EBGraph& ebgraphCoFi = ebisBoxCoFi.getEBGraph();    
    
    Box grownBoxFine = grow(cellBoxFine, m_ghost);
    grownBoxFine &= domainFine;

    // On the fine grid I only want the ghost cells that are irregular. 
    IntVectSet irregIVSFine = ebisBoxFine.getIrregIVS(grownBoxFine);
    irregIVSFine -= cellBoxFine;
    m_vofIterFineGhosts[dit()].define(irregIVSFine, ebgraphFine);

    // On the coarse grid I want coarsenings of the irregular ghost cells on the fine level.
    IntVectSet ivsCoFi = coarsen(irregIVSFine, m_refRat);
    m_vofIterCoFiGhosts[dit()].define(ivsCoFi, ebgraphCoFi);
  }
}

void EBCoarseFineParticleMesh::addFineGhostsToCoarse(LevelData<EBCellFAB>& a_coarData, const LevelData<EBCellFAB>& a_fineData) const {
  CH_TIME("EBCoarseFineParticleMesh::addFineGhostsToCoarse");
  
  CH_assert(m_isDefined);

  // TLDR: This routine will take the ghost cells that are in a_fineData and add them to the coarse level. We use buffers for this,
  //       copying the data from a_fineData to our nifty m_bufferFine buffer holder. The ghost cells in that scratch data are
  //       coarsened onto yet another buffer (m_bufferCoFi). Finally, we add the contents in that buffer to the coarse data.

  const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL   ();
  const ProblemDomain&     domainFine = m_eblgFine.getDomain();
  const EBISLayout&        ebislFine  = m_eblgFine.getEBISL ();

  const DisjointBoxLayout& dblCoFi    = m_eblgCoFi.getDBL   ();
  const ProblemDomain&     domainCoFi = m_eblgCoFi.getDomain();
  const EBISLayout&        ebislCoFi  = m_eblgCoFi.getEBISL ();  
  
  // Copy the fine data to scratch and reset the interior cells. We do this by copying everything to the fine scratch data,
  // and then set all the valid data in each box to zero. The exchange() operation will then take care of ghost cells that. 
  // overlap with valid regions in different boxes (we could use a NeighborIterator to achieve the same). Note that this is critical
  // for the result because we essentially end up setting all valid data to zero so we don't accidentally add mass to invalid region
  // of the coarse grid (i.e., the region underneath the fine grid).
  a_fineData.localCopyTo(m_bufferFine);
  for (DataIterator dit(dblFine); dit.ok(); ++dit){
    const Box& box = dblFine[dit()];
    m_bufferFine[dit()].setVal(0.0, box, 0, m_nComp);
  }
  m_bufferFine.exchange();


  // Coarsen the fine grid data. We do this by AVERAGING the fine-grid data onto the coarse grid. This actually involves entire patches
  // and not just individual ghost cells regions, but that's ok because the rest of the data (in the valid fine regions) are set to zero above,
  // so the coarse data underneath those regions will be zero, also. 
  const Real factor = 1./pow(m_refRat, SpaceDim);  
  for (DataIterator dit(dblFine); dit.ok(); ++dit){
    const EBISBox&   ebisBoxFine = ebislFine   [dit()];
    const EBISBox&   ebisBoxCoFi = ebislCoFi   [dit()];    
    EBCellFAB&       coFiData    = m_bufferCoFi[dit()];
    const EBCellFAB& fineData    = m_bufferFine[dit()];


    coFiData.setVal(0.0);

    // Do all the regular cells. 
    FArrayBox&       coFiDataReg = coFiData.getFArrayBox();
    const FArrayBox& fineDataReg = fineData.getFArrayBox();
    
    const Box fineBox = fineDataReg.box() & domainFine;    
    for (int comp = 0; comp < m_nComp; comp++){
      FORT_AVERAGE_GHOSTS(CHF_FRA1(coFiDataReg, comp),
			  CHF_CONST_FRA1(fineDataReg, comp),
			  CHF_BOX(fineBox),
			  CHF_CONST_INT(m_refRat),
			  CHF_CONST_REAL(factor));
    }

    // Now do the irregular cells. First reset the coarse cell values and then increment with the fine-grid values.
    VoFIterator& vofitCoar = m_vofIterCoFiGhosts[dit()];

    for (vofitCoar.reset(); vofitCoar.ok(); ++vofitCoar){
      const VolIndex&        coarVof  = vofitCoar();
      const Vector<VolIndex> fineVofs = ebislCoFi.refine(coarVof, m_refRat, dit());

      // Take the average of the fine-vof values. 
      for (int comp = 0; comp < m_nComp; comp++){
	coFiData(coarVof, comp) = 0.0;

	for (int ivof = 0; ivof < fineVofs.size(); ivof++){
	  coFiData(coarVof, comp) += fineData(fineVofs[ivof], comp);
	}

	coFiData(coarVof, comp) *= factor;
      }
    }
  }

  // At this point we have coarsened the data in the ghost regions around the fine grid onto m_bufferCoFi. We should be table to take that data and add
  // it directly to 
  const Interval interv(0, m_nComp-1);
  
  m_bufferCoFi.copyTo(interv, a_coarData, interv, m_copierCoFiToCoar, EBAddOp());
}

void EBCoarseFineParticleMesh::addFiCoDataToFine(LevelData<EBCellFAB>& a_fineData, const LevelData<EBCellFAB>& a_fiCoData) const {
  CH_TIME("EBCoarseFineParticleMesh::addFiCoDataToFine");

  const Interval interv(0, m_nComp-1);
  a_fiCoData.copyTo(interv, a_fineData, interv, m_copierFiCoToFine, EBAddOp());
}

LevelData<EBCellFAB>& EBCoarseFineParticleMesh::getFiCoBuffer() const {
  CH_TIME("EBCoarseFineParticleMesh::getFiCoBuffer");
  
  return m_bufferFiCo;
}

const EBLevelGrid& EBCoarseFineParticleMesh::getEblgFiCo() const {
  CH_TIME("EBCoarseFineParticleMesh::getEblgFiCo");
  
  return m_eblgFiCo;
}

#include <CD_NamespaceFooter.H>
