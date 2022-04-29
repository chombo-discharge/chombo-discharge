/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBFineInterp.cpp
  @brief  Implementation of EBFineInterp.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBFineInterp.H>
#include <CD_NamespaceHeader.H>

EBFineInterp::EBFineInterp() {
  CH_TIME("EBFineInterp::EBFineInterp(weak)");  
  m_isDefined = false;
}

EBFineInterp::EBFineInterp(const EBLevelGrid&        a_eblgFine,
			   const EBLevelGrid&        a_eblgCoar,
			   const int&                a_refRat,
			   const int&                a_nComp,
			   const EBIndexSpace* const a_ebisPtr) {
  CH_TIME("EBFineInterp::EBFineInterp(full)");

  this->define(a_eblgFine, a_eblgCoar, a_refRat, a_nComp, a_ebisPtr);
}

EBFineInterp::~EBFineInterp() {

}

void EBFineInterp::define(const EBLevelGrid&        a_eblgFine,
			  const EBLevelGrid&        a_eblgCoar,
			  const int&                a_refRat,
			  const int&                a_nComp,
			  const EBIndexSpace* const a_ebisPtr) {
  CH_TIME("EBFineInterp::define");

  CH_assert(a_refRat > 1);
  CH_assert(a_nComp  > 0);
  CH_assert(!(a_ebisPtr == nullptr));
  
  EBPWLFineInterp::define(a_eblgFine.getDBL(),
			  a_eblgCoar.getDBL(),
			  a_eblgFine.getEBISL(),
			  a_eblgCoar.getEBISL(),
			  a_eblgCoar.getDomain(),
			  a_refRat,
			  a_nComp,
			  a_ebisPtr);
}

void EBFineInterp::regridNoSlopes(LevelData<EBCellFAB>&       a_fineData,
				  const LevelData<EBCellFAB>& a_coarData,
				  const Interval&             a_variables) {
  CH_TIME("EBFineInterp::regridNoSlopes");

  MayDay::Abort("EBFineInterp::regridNoSlopes -- not implemented (yet)");
}

void EBFineInterp::regridMinMod(LevelData<EBCellFAB>&       a_fineData,
				const LevelData<EBCellFAB>& a_coarData,
				const Interval&             a_variables) {
  CH_TIME("EBFineInterp::regridMinMod");

  EBPWLFineInterp::interpolate(a_fineData, a_coarData, a_variables);
}

#include <CD_NamespaceFooter.H>
