/*!
  @file   ebcoartocoar_redist.H
  @brief  New EBFineToCoarRedist class with faster defines
  @author Robert Marskar
*/

#ifndef _ebcoartocoar_redist_
#define _ebcoartocoar_redist_

#include "EBCoarToCoarRedist.H"

/*!
  @brief EBFineToCoarRedist with faster define function
*/
class ebcoartocoar_redist : public EBCoarToCoarRedist {
public:

  ebcoartocoar_redist();
  ~ebcoartocoar_redist();

  void new_define(const EBLevelGrid&  a_eblgFine,
		  const EBLevelGrid&  a_eblgCoar,
		  const int&          a_nref,
		  const int&          a_nvar,
		  const int&          a_redistRad);
};

#endif
