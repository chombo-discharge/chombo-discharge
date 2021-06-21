/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzOp.cpp
  @brief  Implementation of CD_MFHelmholtzOp.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzOp.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzOp::MFHelmholtzOp(){

}

MFHelmholtzOp::~MFHelmholtzOp(){

}

void MFHelmholtzOp::define(){

}

int MFHelmholtzOp::refToCoarser() {

}

unsigned int MFHelmholtzOp::orderOfAccuracy(void) const {
  return 99;
}

void MFHelmholtzOp::enforceCFConsistency(LevelData<MFCellFAB>& a_coarCorr, const LevelData<MFCellFAB>& a_fineCorr) {

}


#include <CD_NamespaceFooter.H>
