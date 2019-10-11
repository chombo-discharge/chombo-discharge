/*!
  @file   ebcoartofine_redist.H
  @brief  New EBFineToCoarRedist class with faster defines
  @author Robert Marskar
*/

#ifndef _ebcoartofine_redist_
#define _ebcoartofine_redist_

#include "EBCoarToFineRedist.H"

/*!
  @brief EBFineToCoarRedist with faster define function
*/
class ebcoartofine_redist : public EBCoarToFineRedist {
public:

  ebcoartofine_redist();
  ~ebcoartofine_redist();

  void new_define(const EBLevelGrid&  a_eblgFine,
		  const EBLevelGrid&  a_eblgCoar,
		  const int&          a_nref,
		  const int&          a_nvar,
		  const int&          a_redistRad);
};

#endif
