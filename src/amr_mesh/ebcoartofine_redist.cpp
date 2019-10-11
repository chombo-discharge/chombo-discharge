/*!
  @file   ebcoartofine_redist.cpp
  @brief  Implementation of ebcoartofine_redist
  @author Robert Marskar
*/

#include "ebcoartofine_redist.H"

#include <EBArith.H>


ebcoartofine_redist::ebcoartofine_redist() : EBCoarToFineRedist(){

}

ebcoartofine_redist::~ebcoartofine_redist(){

}

void ebcoartofine_redist::new_define(const EBLevelGrid&  a_eblgFine,
				     const EBLevelGrid&  a_eblgCoar,
				     const int&          a_nref,
				     const int&          a_nvar,
				     const int&          a_redistRad){

  return EBCoarToFineRedist::define(a_eblgFine, a_eblgCoar, a_nref, a_nvar, a_redistRad);
}
