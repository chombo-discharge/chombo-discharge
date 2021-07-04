/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBMultigridRestrict.cpp
  @brief  Implementation of CD_EBMultigridRestrict.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBCellFactory.H>

// Our includes
#include <CD_EBMultigridRestrict.H>
#include <CD_EBMultigridRestrictF_F.H>
#include <CD_NamespaceHeader.H>

constexpr int EBMultigridRestrict::m_comp;
constexpr int EBMultigridRestrict::m_nComp;

EBMultigridRestrict::EBMultigridRestrict(){
  m_isDefined = false;
}

EBMultigridRestrict::EBMultigridRestrict(const EBLevelGrid& a_eblgFine,
					 const EBLevelGrid& a_eblgCoar,
					 const int&         a_refRat,
					 const bool&        a_volumeWeighted){
  this->define(a_eblgFine, a_eblgCoar, a_refRat, a_volumeWeighted);
}

void EBMultigridRestrict::define(const EBLevelGrid& a_eblgFine,
				 const EBLevelGrid& a_eblgCoar,
				 const int&         a_refRat,
				 const bool&        a_volumeWeighted){
  m_eblgFine       = a_eblgFine;
  m_eblgCoar       = a_eblgCoar;
  m_refRat         = a_refRat;
  m_volumeWeighted = a_volumeWeighted;

  // Define the caorsened fine grid. Need to be careful here because a_eblgCoar can consist of 1x1 grids and those can't be coarsened further. In that
  // case we run the opposite procedure of what would otherwise appear to be normal...
  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

  if(!m_eblgFine.coarsenable(m_refRat)) MayDay::Error("EBMultigridRestrict -- fine layout is not coarsenable!");

  coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
  m_eblgCoFi.setMaxRefinementRatio(m_refRat);
  
  m_coFiData.define(m_eblgCoFi.getDBL(), m_nComp, IntVect::Zero, EBCellFactory(m_eblgCoFi.getEBISL()));

  // Define stencils
  this->defineStencils();

  m_isDefined = true;
}


EBMultigridRestrict::~EBMultigridRestrict(){

}

void EBMultigridRestrict::defineStencils(){
  const DisjointBoxLayout& coarDBL = m_eblgCoFi.getDBL();
  const DisjointBoxLayout& fineDBL = m_eblgFine.getDBL();  
  
  const EBISLayout& coarEBISL      = m_eblgCoFi.getEBISL();
  const EBISLayout& fineEBISL      = m_eblgFine.getEBISL();
  
  m_avgStencils.define(coarDBL);
  m_vofitIrreg. define(coarDBL);

  const int numFinePerCoar = std::pow(m_refRat, SpaceDim);

  for (DataIterator dit(coarDBL); dit.ok(); ++dit){
    const Box cellBoxCoar      = coarDBL[dit()];
    const Box cellBoxFine      = fineDBL[dit()];
    
    const EBISBox& ebisboxCoar = coarEBISL[dit()];
    const EBISBox& ebisboxFine = fineEBISL[dit()];
    
    const EBGraph& ebgraphCoar = ebisboxCoar.getEBGraph();
    const EBGraph& ebgraphFine = ebisboxFine.getEBGraph();
    
    const IntVectSet irregCoar = ebisboxCoar.getIrregIVS(cellBoxCoar);
    const IntVectSet irregFine = ebisboxFine.getIrregIVS(cellBoxFine);

    BaseIVFAB<VoFStencil>& stencils = m_avgStencils[dit()];
    VoFIterator& vofIter            = m_vofitIrreg[dit()];
    
    stencils.define(irregCoar, ebgraphCoar, m_nComp);
    vofIter. define(irregCoar, ebgraphCoar);

    for (vofIter.reset(); vofIter.ok(); ++vofIter){
      const VolIndex& vofCoar = vofIter();
      const Real kappaCoar    = ebisboxCoar.volFrac(vofCoar);
      
      const Vector<VolIndex>& fineVofs = coarEBISL.refine(vofCoar, m_refRat, dit());

      VoFStencil& sten = stencils(vofCoar, m_comp);
      sten.clear();

      if(m_volumeWeighted){
	if(kappaCoar > 0.0){
	  for (int i = 0; i < fineVofs.size(); i++){
	    const Real kappaFine = ebisboxFine.volFrac(fineVofs[i]);
	    sten.add(fineVofs[i], kappaFine/numFinePerCoar);
	  }
	  sten *= 1./kappaCoar;
	}
	else{ // Just set it as the average if there is no real volume on the coarse grid. 
	  for (int i = 0; i < fineVofs.size(); i++){
	    sten.add(fineVofs[i], 1.0);
	  }
	  sten *= 1./fineVofs.size();
	}
      }
      else{ 
	for (int i = 0; i < fineVofs.size(); i++){
	  sten.add(fineVofs[i], 1./numFinePerCoar);
	}
      }
    }
  }
}

void EBMultigridRestrict::restrict(LevelData<EBCellFAB>& a_coarData, const LevelData<EBCellFAB>& a_fineData, const Interval& a_variables) const {
  CH_assert(m_isDefined);

  const DisjointBoxLayout& coarDBL = m_eblgCoFi.getDBL();

  Box refBox(IntVect::Zero, IntVect::Zero);
  refBox.refine(m_refRat);
  const int numFinePerCoar = refBox.numPts();

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++){
    const Interval srcInterv(m_comp, m_comp);
    const Interval dstInterv(ivar  , ivar  );
    
    for (DataIterator dit(coarDBL); dit.ok(); ++dit){
      const Box coarBox = coarDBL[dit()];
      
      EBCellFAB&       coFiData = m_coFiData[dit()];
      const EBCellFAB& fineData = a_fineData[dit()];

      BaseFab<Real>&       coFiDataReg = coFiData.getSingleValuedFAB();
      const BaseFab<Real>& fineDataReg = fineData.getSingleValuedFAB();

      coFiDataReg.setVal(0.0);
      FORT_RESTRICTREGULAR(CHF_FRA1(coFiDataReg,ivar),
			   CHF_CONST_FRA1(fineDataReg,ivar),
			   CHF_BOX(coarBox),
			   CHF_BOX(refBox),
			   CHF_CONST_INT(numFinePerCoar),
			   CHF_CONST_INT(m_refRat));
      
      VoFIterator& vofit = m_vofitIrreg[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();

	const VoFStencil& sten = m_avgStencils[dit()](vof, m_comp);

	coFiData(vof, m_comp) = 0.0;
	for (int i = 0; i < sten.size(); i++){
	  coFiData(vof, m_comp) += sten.weight(i)*fineData(sten.vof(i), ivar);
	}
      }
    }

    m_coFiData.copyTo(srcInterv, a_coarData, dstInterv);
  }
}

#include <CD_NamespaceFooter.H>
