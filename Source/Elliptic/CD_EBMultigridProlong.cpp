/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBMultigridProlong.cpp
  @brief  Implementation of CD_EBMultigridProlong.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBCellFactory.H>

// Our includes
#include <CD_EBMultigridProlong.H>
#include <CD_EBMultigridDanceF_F.H>
#include <CD_NamespaceHeader.H>

constexpr int EBMultigridProlong::m_comp;
constexpr int EBMultigridProlong::m_nComp;

EBMultigridProlong::EBMultigridProlong(){
  m_isDefined = false;
}

EBMultigridProlong::EBMultigridProlong(const EBLevelGrid& a_eblgFine,
				       const EBLevelGrid& a_eblgCoar,
				       const int&         a_refRat,
				       const bool&        a_volumeWeighted){
  this->define(a_eblgFine, a_eblgCoar, a_refRat, a_volumeWeighted);

  MayDay::Abort("EBMultigridProlong -- this code breaks!!!");
}

void EBMultigridProlong::define(const EBLevelGrid& a_eblgFine,
				const EBLevelGrid& a_eblgCoar,
				const int&         a_refRat,
				const bool&        a_volumeWeighted){
  m_eblgFine       = a_eblgFine;
  m_eblgCoar       = a_eblgCoar;
  m_refRat         = a_refRat;
  m_volumeWeighted = a_volumeWeighted;

  // Refine the coarse grid. This way we have an iterator which we can use on the coarse grid but which can access fine grid data. 
  refine(m_eblgFiCo, m_eblgCoar, m_refRat);

  m_fineData.define(m_eblgFiCo.getDBL(), m_nComp, IntVect::Zero, EBCellFactory(m_eblgFiCo.getEBISL()));

  // Define stencils
  this->defineStencils();

  m_isDefined = true;
}


EBMultigridProlong::~EBMultigridProlong(){

}

void EBMultigridProlong::defineStencils(){
  const DisjointBoxLayout& coarDBL = m_eblgCoar.getDBL();  
  const DisjointBoxLayout& fineDBL = m_eblgFiCo.getDBL();
  
  const EBISLayout& coarEBISL      = m_eblgCoar.getEBISL();
  const EBISLayout& fineEBISL      = m_eblgFiCo.getEBISL();
  
  m_prolongStencils.define(fineDBL);
  m_vofitIrregFine. define(fineDBL);

  const int numFinePerCoar = std::pow(m_refRat, SpaceDim);

  for (DataIterator dit(coarDBL); dit.ok(); ++dit){
    const Box cellBoxCoar      = coarDBL[dit()];
    const Box cellBoxFine      = fineDBL[dit()];

    const EBISBox& ebisboxCoar = coarEBISL[dit()];
    const EBISBox& ebisboxFine = fineEBISL[dit()];
    
    const EBGraph& ebgraphCoar = ebisboxCoar.getEBGraph();
    const EBGraph& ebgraphFine = ebisboxFine.getEBGraph();

    const IntVectSet irregCoar = ebisboxCoar.getIrregIVS(cellBoxCoar);
    const IntVectSet irregFine = refine(irregCoar, m_refRat);

    BaseIVFAB<VoFStencil>& stencils = m_prolongStencils[dit()];
    VoFIterator& vofIter            = m_vofitIrregFine[dit()];
    
    stencils.define(irregFine, ebgraphFine, m_nComp);
    vofIter. define(irregFine, ebgraphFine);

    // Iterate through the coarse VoFs and get all the refined vofs
    for (VoFIterator vofit(irregCoar, ebgraphCoar); vofit.ok(); ++vofit){
      const VolIndex& vofCoar = vofit();

      const Real kappaCoar    = ebisboxCoar.volFrac(vofCoar);

      const Vector<VolIndex>& fineVofs = coarEBISL.refine(vofCoar, m_refRat, dit());
#if 0
      Real sumKappaFine = 0.0;
      for (int i = 0; i < fineVofs.size(); i++){
	sumKappaFine += ebisboxFine.volFrac(fineVofs[i]);
      }
      
      const Real volFactor = kappaCoar/(numFinePerCoar*sumKappaFine);
      
      for (int i = 0; i < fineVofs.size(); i++){
	VoFStencil& sten = stencils(fineVofs[i], m_comp);
	sten.clear();

	if(m_volumeWeighted && volFactor > 0.0){
	  sten.add(vofCoar, volFactor);
	}
	else{
	  sten.add(vofCoar, 1.0);
	}
      }
#endif
    }
  }
}

void EBMultigridProlong::prolong(LevelData<EBCellFAB>& a_fineData, const LevelData<EBCellFAB>& a_coarData, const Interval& a_variables) const {
  CH_assert(m_isDefined);

  const DisjointBoxLayout& fineDBL = m_eblgFiCo.getDBL();

  Box refBox(IntVect::Zero, IntVect::Zero);
  refBox.refine(m_refRat);

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++){
    const Interval srcInterv(m_comp, m_comp);
    const Interval dstInterv(ivar  , ivar  );
    
    for (DataIterator dit(fineDBL); dit.ok(); ++dit){
      const Box coarBox = fineDBL[dit()];

      EBCellFAB&       fineData = m_fineData[dit()];
      const EBCellFAB& coarData = a_coarData[dit()];

      BaseFab<Real>&       fineDataReg = fineData.getSingleValuedFAB();
      const BaseFab<Real>& coarDataReg = coarData.getSingleValuedFAB();

      //      fineData.setVal(0.0);//
      FORT_PROLONGREGULAR(CHF_FRA1(fineDataReg,ivar),
			  CHF_CONST_FRA1(coarDataReg,ivar),
			  CHF_BOX(coarBox),
			  CHF_BOX(refBox),
			  CHF_CONST_INT(m_refRat));
      
      VoFIterator& vofit = m_vofitIrregFine[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex& fineVof = vofit();

	const VoFStencil& sten = m_prolongStencils[dit()](fineVof, m_comp);

	fineData(fineVof, m_comp) = 0.0;
	for (int i = 0; i < sten.size(); i++){
	  const Real     weight  = sten.weight(i);
	  const VolIndex coarVoF = sten.vof(i);
	  
	  fineData(fineVof, m_comp) += weight * coarData(coarVoF, ivar);
	}
      }
    }

    m_fineData.copyTo(srcInterv, a_fineData, dstInterv);
  }
}

#include <CD_NamespaceFooter.H>
