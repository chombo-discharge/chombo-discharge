/*!
  @file   EBFastCoarToCoarRedist.H
  @brief  New EBFineToCoarRedist class with faster define function
  @author Robert Marskar
*/

#ifndef _EBFASTCOARTOCOARREDIST_H_
#define _EBFASTCOARTOCOARREDIST_H_

#include "EBCoarToCoarRedist.H"

#include "CD_NamespaceHeader.H"
/*!
  @brief EBFineToCoarRedist with faster define function
*/
class EBFastCoarToCoarRedist : public EBCoarToCoarRedist {
public:

  EBFastCoarToCoarRedist();
  ~EBFastCoarToCoarRedist();
  
  void define(const EBLevelGrid&                      a_eblgFine,
	      const EBLevelGrid&                      a_eblgCoar,
	      const LayoutData<Vector<LayoutIndex> >& a_neighborsFine,
	      const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar,
	      const int&                              a_nref,
	      const int&                              a_nvar,
	      const int&                              a_redistRad);

protected:

  // Debugging functions
  void gatherBroadcast(IntVectSet& a_set);
  void gatherCoarSet(IntVectSet& coarSet); // Debugging function
};
#include "CD_NamespaceFooter.H"

#endif
