/*!
  @file   EBGhostCloud
  @brief  A class for adding coarsened ghost cells to a coarse level
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _EBGHOSTCLOUD_
#define _EBGHOSTCLOUD_

#include "DisjointBoxLayout.H"
#include "EBLevelGrid.H"
#include "EBCellFAB.H"
#include "ProblemDomain.H"

class EBGhostCloud {
public:

  /*!
    @brief Default constructor. Needs to call define
  */
  EBGhostCloud();

  /*!
    @brief Full constructor. Calls define. 
  */
  EBGhostCloud(const DisjointBoxLayout& a_gridsCoar,
	       const DisjointBoxLayout& a_gridsFine,
	       const EBLevelGrid&       a_eblgCoar,
	       const EBLevelGrid&       a_eblgFine,
	       const ProblemDomain&     a_domainCoar,
	       const ProblemDomain&     a_domainFine,
	       const int                a_ref_rat,
	       const int                a_nComp,
	       const int                a_ghost);

  ~EBGhostCloud();

  void define(const DisjointBoxLayout& a_gridsCoar,
	      const DisjointBoxLayout& a_gridsFine,
	      const EBLevelGrid&       a_eblgCoar,
	      const EBLevelGrid&       a_eblgFine,
	      const ProblemDomain&     a_domainCoar,
	      const ProblemDomain&     a_domainFine,
	      const int                a_ref_rat,
	      const int                a_nComp,
	      const int                a_ghost);

  void addGhostsToCoar(LevelData<EBCellFAB>& a_coarData, const LevelData<EBCellFAB>& a_fineData);

protected:

  bool                     m_isDefined;
  BoxLayout                m_gridsCoFi;
  BoxLayoutData<EBCellFAB> m_dataCoFi;
  DisjointBoxLayout        m_gridsFine;
  EBLevelGrid              m_eblgCoFi;
  LevelData<EBCellFAB>     m_scratchFine;
  ProblemDomain            m_domainCoar;
  int                      m_refRat;
  int                      m_nComp;
  int                      m_ghost;
};
#endif
