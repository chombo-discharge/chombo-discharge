#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves fri aug 17, 2001
//major rework 4/2009 dtg

#include "ebcoarseaverage.H"
#include "EBAverageF_F.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "DisjointBoxLayout.H"
#include "EBCellFactory.H"
#include "EBFluxFactory.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "CH_Timer.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "EBLoadBalance.H"
//#include "NamespaceHeader.H"

#include "CD_NamespaceHeader.H"

/************************************/
ebcoarseaverage::ebcoarseaverage()
{
  EBCoarseAverage::setDefaultValues();
}
/************************************/
ebcoarseaverage::~ebcoarseaverage()
{
}
/************************************/
ebcoarseaverage::ebcoarseaverage(const DisjointBoxLayout& a_dblFine,
				 const DisjointBoxLayout& a_dblCoar,
				 const EBISLayout& a_ebislFine,
				 const EBISLayout& a_ebislCoar,
				 const ProblemDomain& a_domainCoar,
				 const int& a_nref,
				 const int& a_nvar,
				 const EBIndexSpace* ebisPtr)
{
  EBCoarseAverage::setDefaultValues();

  EBCoarseAverage::define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
	 a_domainCoar, a_nref, a_nvar, ebisPtr);
}
/************************************/
ebcoarseaverage::ebcoarseaverage(const EBLevelGrid& a_eblgFine,
				 const EBLevelGrid& a_eblgCoar,
				 const EBLevelGrid& a_eblgCoFi,
				 const int& a_nref,
				 const int& a_nvar)
{
  EBCoarseAverage::setDefaultValues();

  EBCoarseAverage::define(a_eblgFine, a_eblgCoar, a_eblgCoFi, a_nref, a_nvar);
}

void ebcoarseaverage::conservativeAverage(LevelData<BaseIVFAB<Real> >&        a_coarData,
					  const LevelData<BaseIVFAB<Real> >&  a_fineData,
					  const Interval&                     a_variables){
  CH_TIME("ebcoarseaverage::average(LD<BaseIVFAB>)");
  LevelData<BaseIVFAB<Real> > coarFiData;
  LevelData<BaseIVFAB<Real> > fineBuffer;
  CH_assert(isDefined());
  {
    CH_TIME("buffer allocation");
    BaseIVFactory<Real> factCoFi(m_eblgCoFi.getEBISL(), m_irregSetsCoFi);
    coarFiData.define(m_eblgCoFi.getDBL(), m_nComp, IntVect::Zero, factCoFi);
    if (m_useFineBuffer)
      {
	BaseIVFactory<Real> factFine(m_eblgFine.getEBISL(), m_irregSetsFine);
	coarFiData.define(m_eblgFine.getDBL(), m_nComp, IntVect::Zero, factFine);
      }
  }

  if (m_useFineBuffer)
    {
      CH_TIME("fine_copy");
      a_fineData.copyTo(a_variables, fineBuffer, a_variables);
    }

  {
    CH_TIME("averaging");
    int ifnerg = 0;
    for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit)
      {
	const BaseIVFAB<Real>* fineFABPtr = NULL;
	if (m_useFineBuffer)
	  {
	    fineFABPtr = &fineBuffer[dit()];
	  }
	else
	  {
	    fineFABPtr = &a_fineData[dit()];
	  }
	BaseIVFAB<Real>&       cofiFAB = coarFiData[dit()];
	const BaseIVFAB<Real>& fineFAB = *fineFABPtr;
	conservativeAverageFAB(cofiFAB,
			       fineFAB,
			       dit(),
			       a_variables);

	ifnerg++;
      }
  }
  {
    CH_TIME("copy_coar");
    coarFiData.copyTo(a_variables, a_coarData, a_variables);
  }
}

void ebcoarseaverage::conservativeAverageFAB(BaseIVFAB<Real>&       a_coar,
					     const BaseIVFAB<Real>& a_fine,
					     const DataIndex&       a_datInd,
					     const Interval&        a_variables) const
{
  CH_assert(isDefined());
  //recall that datInd is from the fine layout.
  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_datInd];
  const IntVectSet& coarIrregIVS = a_coar.getIVS();
  const IntVectSet& fineIrregIVS = a_fine.getIVS();

  const int nref2 = pow(m_refRat, SpaceDim-1);

  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph());
       vofitCoar.ok(); ++vofitCoar)
    {
      const VolIndex& coarVoF = vofitCoar();
      Vector<VolIndex> fineVoFs =
	m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

      for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
	{
	  int  numVoFs = 0;
	  Real areaTot = 0;
	  Real dataVal = 0;
	  for (int ifine = 0; ifine < fineVoFs.size(); ifine++)
	    {
	      const VolIndex& fineVoF = fineVoFs[ifine];
	      if (fineIrregIVS.contains(fineVoF.gridIndex()))
		{
		  Real bndryArea = ebisBoxFine.bndryArea(fineVoF);
		  if (bndryArea > 0)
		    {
		      areaTot += bndryArea;
		      numVoFs++;
		      dataVal += bndryArea*a_fine(fineVoF, ivar);
		    }
		}
	    }
	  const Real safety = 1.E-8;
	  a_coar(coarVoF, ivar) = dataVal/(nref2*(safety + ebisBoxCoar.bndryArea(coarVoF)));
	}
    }
}


#include "CD_NamespaceFooter.H"
