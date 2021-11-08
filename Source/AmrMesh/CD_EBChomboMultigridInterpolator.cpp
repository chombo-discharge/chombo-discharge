/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBChomboMultigridInterpolator.cpp
  @brief  Implementation of CD_EBChomboMultigridInterpolator.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBChomboMultigridInterpolator.H>
#include <CD_NamespaceHeader.H>


EBChomboMultigridInterpolator::EBChomboMultigridInterpolator(const EBLevelGrid& a_eblgFine,
							     const EBLevelGrid& a_eblgCoar,
							     const int          a_refRat,
							     const int          a_nVar){
  CH_TIME("EBChomboMultigridInterpolator::EBChomboMultigridInterpolator");
}

EBChomboMultigridInterpolator::~EBChomboMultigridInterpolator(){
  CH_TIME("EBChomboMultigridInterpolator::~EBChomboMultigridInterpolator");
}


int EBChomboMultigridInterpolator::getGhostCF() const {
  return 1;
}


void EBChomboMultigridInterpolator::coarseFineInterp(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phiCoar, const Interval a_variables) {
  CH_TIME("EBChomboMultigridInterpolator::coarseFineInterp");

  MayDay::Error("EBChomboMultigridInterpolator::coarseFineInterp -- not implemented");
}

void EBChomboMultigridInterpolator::coarseFineInterpH(LevelData<EBCellFAB>& a_phiFine, const Interval a_variables) const {
  CH_TIME("EBChomboMultigridInterpolator::coarseFineInterpH(patch)");

  MayDay::Error("EBChomboMultigridInterpolator::coarseFineInterpH -- not implemented");

}

void EBChomboMultigridInterpolator::coarseFineInterpH(EBCellFAB& a_phiFine, const Interval a_variables, const DataIndex& a_dit) const {
  CH_TIME("EBChomboMultigridInterpolator::coarseFineInterpH(patch)");

  MayDay::Error("EBChomboMultigridInterpolator::coarseFineInterpH -- not implemented");    
}

#include <CD_NamespaceFooter.H>
