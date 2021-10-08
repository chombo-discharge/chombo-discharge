/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbCoarAve.cpp
  @brief  Implementation of CD_EbCoarAve.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBFluxFactory.H>
#include <EBAverageF_F.H>

// Our includes
#include <CD_EbCoarAve.H>
#include <CD_NamespaceHeader.H>

EbCoarAve::EbCoarAve(){
  EBCoarseAverage::setDefaultValues();
}

EbCoarAve::~EbCoarAve(){
}

EbCoarAve::EbCoarAve(const DisjointBoxLayout& a_dblFine,
		     const DisjointBoxLayout& a_dblCoar,
		     const EBISLayout&        a_ebislFine,
		     const EBISLayout&        a_ebislCoar,
		     const ProblemDomain&     a_domainCoar,
		     const int&               a_refRat,
		     const int&               a_nComp,
		     const EBIndexSpace*      a_ebisPtr){
  CH_TIME("EbCoarAve::EbCoarAve");
  
  EBCoarseAverage::setDefaultValues();
  EBCoarseAverage::define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
			  a_domainCoar, a_refRat, a_nComp, a_ebisPtr);
}

EbCoarAve::EbCoarAve(const EBLevelGrid& a_eblgFine,
		     const EBLevelGrid& a_eblgCoar,
		     const EBLevelGrid& a_eblgCoFi,
		     const int&         a_refRat,
		     const int&         a_nComp){
  CH_TIME("EbCoarAve::EbCoarAve");
  
  EBCoarseAverage::setDefaultValues();
  EBCoarseAverage::define(a_eblgFine, a_eblgCoar, a_eblgCoFi, a_refRat, a_nComp);
}

void EbCoarAve::conservativeAverage(LevelData<BaseIVFAB<Real> >&        a_coarData,
				    const LevelData<BaseIVFAB<Real> >&  a_fineData,
				    const Interval&                     a_variables){
  CH_TIME("EbCoarAve::conservativeAverage(LD<BaseIVFAB>)");
  
  LevelData<BaseIVFAB<Real> > coarFiData;
  LevelData<BaseIVFAB<Real> > fineBuffer;
  
  CH_assert(isDefined());
  {
    CH_TIME("buffer allocation");
    BaseIVFactory<Real> factCoFi(m_eblgCoFi.getEBISL(), m_irregSetsCoFi);
    coarFiData.define(m_eblgCoFi.getDBL(), m_nComp, IntVect::Zero, factCoFi);
    if (m_useFineBuffer){
      BaseIVFactory<Real> factFine(m_eblgFine.getEBISL(), m_irregSetsFine);
      coarFiData.define(m_eblgFine.getDBL(), m_nComp, IntVect::Zero, factFine);
    }
  }

  if (m_useFineBuffer){
    CH_TIME("fine_copy");
    a_fineData.copyTo(a_variables, fineBuffer, a_variables);
  }

  {
    CH_TIME("averaging");
    for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit){
      const BaseIVFAB<Real>* fineFABPtr = NULL;
      if (m_useFineBuffer){
	fineFABPtr = &fineBuffer[dit()];
      }
      else{
	fineFABPtr = &a_fineData[dit()];
      }
      BaseIVFAB<Real>&       cofiFAB = coarFiData[dit()];
      const BaseIVFAB<Real>& fineFAB = *fineFABPtr;
      
      this->conservativeAverageFAB(cofiFAB,
				   fineFAB,
				   dit(),
				   a_variables);
    }
  }
  {
    CH_TIME("copy_coar");
    coarFiData.copyTo(a_variables, a_coarData, a_variables);
  }
}

void EbCoarAve::conservativeAverageFAB(BaseIVFAB<Real>&       a_coar,
				       const BaseIVFAB<Real>& a_fine,
				       const DataIndex&       a_datInd,
				       const Interval&        a_variables) const{
  CH_assert(isDefined());

  const EBISBox& ebisBoxCoar     = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine     = m_eblgFine.getEBISL()[a_datInd];
  
  const IntVectSet& coarIrregIVS = a_coar.getIVS();
  const IntVectSet& fineIrregIVS = a_fine.getIVS();

  const int nref2 = std::pow(m_refRat, SpaceDim-1);

  // This loop computes phiCoar = sum(phiFine*areaFine)/areaCoar. If areaCoar == 0
  // then we use a safety factor here. 
  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph()); vofitCoar.ok(); ++vofitCoar) {
    const VolIndex& coarVoF         = vofitCoar();
    const Vector<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

    for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++){
      int  numVoFs = 0;
      Real areaTot = 0;
      Real dataVal = 0;
	
      for (int ifine = 0; ifine < fineVoFs.size(); ifine++){
	const VolIndex& fineVoF = fineVoFs[ifine];
	
	if (fineIrregIVS.contains(fineVoF.gridIndex())){
	  Real bndryArea = ebisBoxFine.bndryArea(fineVoF);
	  if (bndryArea > 0){
	    areaTot += bndryArea;
	    numVoFs++;
	    dataVal += bndryArea*a_fine(fineVoF, ivar);
	  }
	}
      }
      
      constexpr Real safety = 1.E-8;
      
      a_coar(coarVoF, ivar) = dataVal/(nref2*(safety + ebisBoxCoar.bndryArea(coarVoF)));
    }
  }
}

void EbCoarAve::averageFaceData(LevelData<EBFluxFAB>&       a_coarData,
				const LevelData<EBFluxFAB>& a_fineData,
				const Interval&             a_variables){
  CH_TIME("EbCoarAve::averageFaceData");

  // TLDR: This routine computes the arithmetic average of the fine faces. This means that the coarse face value is
  //
  //          phi(coarFace) = sum[phi(fineFace)]/num(fineFaces)  

  CH_assert(isDefined());

  LevelData<EBFluxFAB> coarFiData;
  LevelData<EBFluxFAB> fineBuffer;

  EBFluxFactory factCoFi(m_eblgCoFi.getEBISL());
  EBFluxFactory factFine(m_eblgFine.getEBISL());

  // Define buffer data.
  const int numComp = a_variables.end() + 1;  
  coarFiData.define(m_eblgCoFi.getDBL(), numComp, IntVect::Zero, factCoFi);
  if (m_useFineBuffer){
    fineBuffer.define(m_eblgFine.getDBL(), numComp, IntVect::Zero, factFine);
    a_fineData.copyTo(a_variables, fineBuffer, a_variables);    
  }

  // Grid loop. Switch between buffer data or not. 
  for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit) {
    const EBFluxFAB* fineFABPtr = nullptr;
      
    if (m_useFineBuffer){
      fineFABPtr = &fineBuffer[dit()];
    }
    else{
      fineFABPtr = &a_fineData[dit()];
    }
      
    EBFluxFAB&       cofiFAB = coarFiData[dit()];
    const EBFluxFAB& fineFAB = *fineFABPtr;
      
    for (int idir = 0; idir < SpaceDim; idir++){
      this->averageFaces(cofiFAB[idir],
			 fineFAB[idir],
			 dit(),
			 a_variables,
			 a_variables,
			 idir);
    }
  }

  // Copy back to input data. 
  coarFiData.copyTo(a_variables, a_coarData, a_variables);
}

void EbCoarAve::averageFaces(EBFaceFAB&       a_coarData,
			     const EBFaceFAB& a_fineData,
			     const DataIndex& a_datInd,
			     const Interval&  a_fineInterv,
			     const Interval&  a_coarInterv,
			     const int&       a_dir) {
  CH_TIME("EBCoarseAverage::averageFaces");
  
  CH_assert(isDefined());
  CH_assert(a_fineInterv.size() == a_coarInterv.size());

  // TLDR: This routine computes the arithmetic average of the fine faces. This means that the coarse face value is
  //
  //          phi(coarFace) = sum[phi(fineFace)]/num(fineFaces)


  //do all cells as if they were regular
  BaseFab<Real>&       coarRegFAB = a_coarData.getSingleValuedFAB();
  const BaseFab<Real>& fineRegFAB = a_fineData.getSingleValuedFAB();

  const Box& coarDBLBox = m_eblgCoFi.getDBL().get(a_datInd);
  
  Box refbox(IntVect::Zero,(m_refRat-1)*IntVect::Unit);
  Box coarFaceBox = coarDBLBox;
  
  refbox.     surroundingNodes(a_dir);
  coarFaceBox.surroundingNodes(a_dir);

  for (int ioff = 0; ioff < a_fineInterv.size(); ioff++){
    const int ivarf = a_fineInterv.begin() + ioff;
    const int ivarc = a_coarInterv.begin() + ioff;
    FORT_EBAVERAGEFACE(CHF_FRA1(coarRegFAB,ivarc),
		       CHF_CONST_FRA1(fineRegFAB,ivarf),
		       CHF_BOX(coarFaceBox),
		       CHF_CONST_INT(m_refRat),
		       CHF_BOX(refbox),
		       CHF_CONST_INT(a_dir));
  }

  // Now do the irregular faces. 
  const EBISBox&   ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const EBGraph&   ebgraphCoar  = ebisBoxCoar.getEBGraph();
  const IntVectSet coarIrregIVS = ebisBoxCoar.getIrregIVS(coarDBLBox);

  for (FaceIterator faceItCoar(coarIrregIVS, ebgraphCoar, a_dir, FaceStop::SurroundingWithBoundary); faceItCoar.ok(); ++faceItCoar){
    const FaceIndex&  coarFace  = faceItCoar();
    
    const std::vector<FaceIndex> fineFaces = m_eblgCoFi.getEBISL().refine(coarFace, m_refRat, a_datInd).stdVector();

    const int numFineFaces = fineFaces.size();

    for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
      const int ivarf = a_fineInterv.begin() + ioff;
      const int ivarc = a_coarInterv.begin() + ioff;

      Real dataVal = 0.0;
      if (numFineFaces > 0){
	for (const auto& fineFace : fineFaces){
	  dataVal += a_fineData(fineFace, ivarf);
	}

	dataVal /= Real(fineFaces.size());
      }

      a_coarData(coarFace, ivarc) = dataVal;
    }
  }
}

#include <CD_NamespaceFooter.H>
