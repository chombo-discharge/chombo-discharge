/*!
  @file   ebfinetocoar_redist.H
  @brief  New EBFineToCoarRedist class with faster defines
  @author Robert Marskar
*/

#ifndef _ebfinetocoar_redist_
#define _ebfinetocoar_redist_

#include "EBFineToCoarRedist.H"

/*!
  @brief EBFineToCoarRedist with faster define function
*/
class ebfinetocoar_redist : public EBFineToCoarRedist {
public:

  ebfinetocoar_redist();
  ~ebfinetocoar_redist();

  void new_define(const EBLevelGrid&  a_eblgFine,
		  const EBLevelGrid&  a_eblgCoar,
		  const int&          a_nref,
		  const int&          a_nvar,
		  const int&          a_redistRad);
};

#endif
