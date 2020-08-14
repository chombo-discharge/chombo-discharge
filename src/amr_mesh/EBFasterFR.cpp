/*!
  @brief   EBFasterFR.cpp
  @details Implementation of EBFasterFR.H
  @author  Robert Maskar
  @date    Aug. 2020
*/

#include "EBFasterFR.H"

EBFasterFR::EBFasterFR(){

}

EBFasterFR::EBFasterFR(const EBLevelGrid&       a_eblgFine,
		       const EBLevelGrid&       a_eblgCoar,
		       const int&               a_nref,
		       const int&               a_nvar,
		       const bool               a_forceNoEBCF){
  setDefaultValues();
  define(a_eblgFine, a_eblgCoar, a_nref, a_nvar, a_forceNoEBCF);
}

EBFasterFR::~EBFasterFR(){

}

void EBFasterFR::define(const EBLevelGrid&       a_eblgFine,
			const EBLevelGrid&       a_eblgCoar,
			const int&               a_nref,
			const int&               a_nvar,
			const bool               a_forceNoEBCF){
  EBFastFR::define(a_eblgFine, a_eblgCoar, a_nref, a_nvar, a_forceNoEBCF);
}

