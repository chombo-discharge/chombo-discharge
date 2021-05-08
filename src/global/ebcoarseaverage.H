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

#ifndef _ebcoarseaverage_H_
#define _ebcoarseaverage_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBLevelGrid.H"
#include "EBCellFAB.H"
#include "EBFluxFAB.H"
#include "BaseIVFAB.H"
#include "BaseIVFactory.H"
#include "EBFaceFAB.H"
#include "Interval.H"
#include "NamespaceHeader.H"
class EBIndexSpace;

#include "CD_NamespaceHeader.H"
/// replaces coarse level data with an average of fine level data.
/**
   This class replaces data at a coarse level of refinement with an
   average of data at a finer level of refinement, in areas where fine
   data is present.  Coarse level data is not modified where fine
   level data is not present.
*/
class ebcoarseaverage
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  ebcoarseaverage();

  ///
  ~ebcoarseaverage();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().
  */
  ebcoarseaverage(const DisjointBoxLayout& dblFine,
		  const DisjointBoxLayout& dblCoar,
		  const EBISLayout& ebislFine,
		  const EBISLayout& ebislCoar,
		  const ProblemDomain& domainCoar,
		  const int& nref,
		  const int& nvar,
		  const EBIndexSpace* ebisPtr);

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().
  */
  ebcoarseaverage(const EBLevelGrid& a_eblgFine,
		  const EBLevelGrid& a_eblgCoar,
		  const EBLevelGrid& a_eblgCoFi,
		  const int& nref,
		  const int& nvar);
  ///
  /**
     Defines this object.  Existing information is overriden.

  */
  void define(const DisjointBoxLayout& dblFine,
	      const DisjointBoxLayout& dblCoar,
	      const EBISLayout& ebislFine,
	      const EBISLayout& ebislCoar,
	      const ProblemDomain& domainCoar,
	      const int& nref,
	      const int& nvar,
	      const EBIndexSpace* ebisPtr);

  ///
  /**
     Defines this object.  Existing information is overriden.

  */
  void define(const EBLevelGrid& a_eblgFine,
	      const EBLevelGrid& a_eblgCoar,
	      const EBLevelGrid& a_eblgCoFi,
	      const int& nref,
	      const int& nvar);
  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.
  */
  bool isDefined() const;

  ///
  /**
     Replaces a_coarse_data with the average of a_fine_data within the
     coarsening of the domain of the fine level.  Elsewhere, the
     coarse data is unchanged.  It is an error to call if not
     this->isDefined(). The domain of a_fine_data should be
     the same as the fine domain specified in the most recent call to
     define().  It is expected that the coarse and fine level's
     domains are properly nested.  Both a_coarse_data and a_fine_data
     should have the same number of components specified in the most
     recent call to define().

  */
  void
  average(LevelData<EBCellFAB>&       a_coarData,
	  const LevelData<EBCellFAB>& a_fineData,
	  const Interval&             a_variables);

  void
  average(LevelData<EBFluxFAB>&       a_coarData,
	  const LevelData<EBFluxFAB>& a_fineData,
	  const Interval&             a_variables);


  void
  average(LevelData<BaseIVFAB<Real> >&       a_coarData,
	  const LevelData<BaseIVFAB<Real> >& a_fineData,
	  const Interval&                    a_variables);

  void
  conservative_average(LevelData<BaseIVFAB<Real> >&       a_coarData,
		       const LevelData<BaseIVFAB<Real> >& a_fineData,
		       const Interval&                    a_variables);

protected:
  void
  averageFAB(EBCellFAB&       a_coar,
	     const EBCellFAB& a_fine,
	     const DataIndex& a_datInd,
	     const Interval&  a_variables) const;

  void
  averageFAB(EBFaceFAB&       a_coar,
	     const EBFaceFAB& a_fine,
	     const DataIndex& a_datInd,
	     const Interval&  a_variables,
	     const int&       a_dir) const;

  void
  averageFAB(BaseIVFAB<Real>&       a_coar,
	     const BaseIVFAB<Real>& a_fine,
	     const DataIndex&       a_datInd,
	     const Interval&        a_variables) const;

  void
  conservative_averageFAB(BaseIVFAB<Real>&       a_coar,
			  const BaseIVFAB<Real>& a_fine,
			  const DataIndex&       a_datInd,
			  const Interval&        a_variables) const;


  void
  setDefaultValues();

  bool m_isDefined;
  //true only if the layout is not coarsenable
  //and the grid covers the domain.   This
  //happens in multigrid sometimes
  bool m_useFineBuffer;
  EBLevelGrid m_eblgCoar;
  EBLevelGrid m_eblgFine;
  EBLevelGrid m_eblgCoFi;
  LayoutData<IntVectSet> m_irregSetsCoFi;

  //only defined if using fine buffer
  LayoutData<IntVectSet> m_irregSetsFine;

  int m_refRat;
  int m_nComp;

private:
  //disallowed for all the usual reasons
  ebcoarseaverage(const ebcoarseaverage& ebcin)
  {
    MayDay::Error("ebca 2 invalid operator");
  }
  void operator=(const ebcoarseaverage& fabin)
  {
    MayDay::Error("ebca 3 invalid operator");
  }

};
#include "CD_NamespaceFooter.H"

#include "NamespaceFooter.H"
#endif
