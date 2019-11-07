/*!
  @file   EBLevelConcentrationRedist.H
  @brief  Declaration file for EBLevelConcentrationRedist
  @author Robert Marskar
*/

#ifndef _EBLEVELCONCENTRATIONREDIST_
#define _EBLEVELCONCENTRATIONREDIST_

#include "DisjointBoxLayout.H"
#include "ProblemDomain.H"
#include "EBISLayout.H"
#include "ConcentrationStencil.H"


/*!
  @brief Strange class that does local redistribution from irregular cells into regular cells in order 
         to maintain non-negative densities in cut cells
*/
class EBLevelConcentrationRedist {
public:

  /*!
    @brief Default constructor. User must call define
  */
  EBLevelConcentrationRedist();

  /*!
    @brief Deconstructor
  */  
  ~EBLevelConcentrationRedist();

  /*!
    @brief Defining constructor
  */
  EBLevelConcentrationRedist(const DisjointBoxLayout& a_dbl,
			     const EBISLayout&        a_ebisl,
			     const ProblemDomain&     a_domain,
			     const int                a_ncomp,
			     const int                a_redistRad);

  /*!
    @brief Defines object. This is called by the defining constructor
  */
  void define(const DisjointBoxLayout& a_dbl,
	      const EBISLayout&        a_ebisl,
	      const ProblemDomain&     a_domain,
	      const int                a_ncomp,
	      const int                a_redistRad);

  /*!
    @brief Increment function
  */
  void increment(const BaseIVFAB<Real>& a_massDiff, const DataIndex& a_dit6, const Interval& a_variables);

  /*!
    @brief Redistribution
  */
  void redistribute(LevelData<EBCellFAB>& a_solution, const Interval& a_variables);
		 
  /*!
    @brief Redistribution
  */
  void redistribute(LevelData<EBCellFAB>& a_solution, const Interval& a_srcVariables, const Interval& a_dstVariables);

  /*!
    @brief Modify redistribution weights
  */
  void resetWeights(const LevelData<EBCellFAB>& a_modifier, const int a_ivar);

  /*!
    @brief Set to zero
  */
  void setToZero();

protected:

  bool m_isDefined;

  int m_redistRad;
  int m_ncomp;
  
  DisjointBoxLayout m_grids;
  ProblemDomain     m_domain;
  EBISLayout        m_ebisl;

  LevelData<BaseIVFAB<Real> > m_buffer;
  LayoutData<IntVectSet>      m_sets;

  ConcentrationStencil m_stencil;

private:

  // Disallowed for the usual reasons
  EBLevelConcentrationRedist(const EBLevelConcentrationRedist& a_input){
    MayDay::Abort("EBLevelConcentrationRedist(EBLevelRedist) - invalid operator");
  }

};

#endif
