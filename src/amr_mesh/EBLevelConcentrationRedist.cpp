/*!
  @file   EBLevelConcentrationRedist.cpp
  @brief  Declaration file for EBLevelConcentrationRedist
  @author Robert Marskar
*/


#include "EBLevelConcentrationRedist.H"

EBLevelConcentrationRedist::EBLevelConcentrationRedist(){
  m_isDefined = false;
}


EBLevelConcentrationRedist::~EBLevelConcentrationRedist(){

}

EBLevelConcentrationRedist::EBLevelConcentrationRedist(const DisjointBoxLayout& a_dbl,
						       const EBISLayout&        a_ebisl,
						       const ProblemDomain&     a_domain,
						       const int                a_ncomp,
						       const int                a_redistRad){
  CH_TIME("EBLevelConcentrationRedist::EBLevelConcentrationRedist");
  define(a_dbl, a_ebisl, a_domain, a_ncomp, a_redistRad);
}


void EBLevelConcentrationRedist::define(const DisjointBoxLayout& a_dbl,
					const EBISLayout&        a_ebisl,
					const ProblemDomain&     a_domain,
					const int                a_ncomp,
					const int                a_redistRad){
  CH_TIME("EBLevelConcentrationRedist::define");
  m_isDefined = true;
  m_grids     = a_dbl;
  m_ebisl     = a_ebisl;
  m_domain    = a_domain;

  
}

void EBLevelConcentrationRedist::resetWeights(const LevelData<EBCellFAB>& a_modifier, const int a_ivar){
  CH_TIME("EBLevelConcentrationRedist::resetWeights");
}

void EBLevelConcentrationRedist::increment(const BaseIVFAB<Real>& a_massDiff, const DataIndex& a_dit, const Interval& a_variables){
  CH_TIME("EBLevelConcentrationRedist::increment");
}

void EBLevelConcentrationRedist::redistribute(LevelData<EBCellFAB>& a_solution, const Interval& a_variables){
  CH_TIME("EBLevelConcentrationRedist::redistribute(LDEBCellFAB, Interval");

  this->redistribute(a_solution, a_variables, a_variables);
}
		 
void EBLevelConcentrationRedist::redistribute(LevelData<EBCellFAB>& a_solution,
					      const Interval& a_srcVariables,
					      const Interval& a_dstVariables){
    CH_TIME("EBLevelConcentrationRedist::redistribute(LDEBCellFAB, Interval, Interval");

}

void EBLevelConcentrationRedist::setToZero() {
  CH_TIME("EBLevelConcentrationRedist::setToZero");
  CH_assert(m_isDefined);
  
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    //    m_buffer[dit()].setVal(0.0);
  }
}
