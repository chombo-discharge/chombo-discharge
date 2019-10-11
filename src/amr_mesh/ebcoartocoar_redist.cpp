/*!
  @file   ebcoartocoar_redist.cpp
  @brief  Implementation of ebcoartocoar_redist
  @author Robert Marskar
*/

#include "ebcoartocoar_redist.H"

#include <EBArith.H>


ebcoartocoar_redist::ebcoartocoar_redist() : EBCoarToCoarRedist(){

}

ebcoartocoar_redist::~ebcoartocoar_redist(){

}

void ebcoartocoar_redist::new_define(const EBLevelGrid&  a_eblgFine,
				     const EBLevelGrid&  a_eblgCoar,
				     const int&          a_nref,
				     const int&          a_nvar,
				     const int&          a_redistRad){

  return EBCoarToCoarRedist::define(a_eblgFine, a_eblgCoar, a_nref, a_nvar, a_redistRad);
}
