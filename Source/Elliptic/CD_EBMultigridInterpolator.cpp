/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBMultigridInterpolator.cpp
  @brief  Implementation of CD_EBMultigridInterpolator.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBCellFactory.H>
#include <EBLevelDataOps.H>

// Our includes
#include <CD_EBMultigridInterpolator.H>
#include <CD_NamespaceHeader.H>

EBMultigridInterpolator::EBMultigridInterpolator(){
  m_isDefined = false;
}

EBMultigridInterpolator::EBMultigridInterpolator(const EBLevelGrid&            a_eblgFine,
						 const EBLevelGrid&            a_eblgCoar,
						 const int                     a_nRef,
						 const int                     a_nVar,
						 const LayoutData<IntVectSet>& a_ghostCells)
  : EBQuadCFInterp(a_eblgFine.getDBL(),
		   a_eblgCoar.getDBL(),
		   a_eblgFine.getEBISL(),
		   a_eblgCoar.getEBISL(),
		   a_eblgCoar.getDomain(),
		   a_nRef,
		   a_nVar,
		   a_ghostCells,
		   a_eblgFine.getEBIS(),
		   true){

  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;
  
  // Define a temp which is zero everywhere. Don't need ghost cells because the stencils should(!) only reach into valid cells AFAIK. 
  EBCellFactory cellFact(m_eblgCoar.getEBISL());
  m_zeroCoar.define(m_eblgCoar.getDBL(), a_nVar, IntVect::Zero, cellFact);
  EBLevelDataOps::setToZero(m_zeroCoar);
}

EBMultigridInterpolator::~EBMultigridInterpolator(){

}

void EBMultigridInterpolator::coarseFineInterp(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phiCoar, const Interval a_variables){
  EBQuadCFInterp::interpolate(a_phiFine, a_phiCoar, a_variables, false);
}


void EBMultigridInterpolator::coarseFineInterpH(LevelData<EBCellFAB>& a_phiFine, const Interval a_variables){
  EBQuadCFInterp::interpolate(a_phiFine, m_zeroCoar, a_variables, false);
}

#include <CD_NamespaceFooter.H>
