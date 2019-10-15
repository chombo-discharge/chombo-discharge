/*!
  @file   EBFastCoarToFineRedist.H
  @brief  New EBFineToCoarRedist class with faster defines
  @author Robert Marskar
*/

#ifndef _EBFASTCOARTOFINEREDIST_H_
#define _EBFASTCOARTOFINEREDIST_H_

#include "EBCoarToFineRedist.H"

/*!
  @brief EBFineToCoarRedist with faster define function
*/
class EBFastCoarToFineRedist : public EBCoarToFineRedist {
public:

  EBFastCoarToFineRedist();
  ~EBFastCoarToFineRedist();

  void define(const EBLevelGrid&                      a_eblgFine,
	      const EBLevelGrid&                      a_eblgCoar,
	      const LayoutData<Vector<LayoutIndex> >& a_neighborsFine,
	      const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar,
	      const int&                              a_nref,
	      const int&                              a_nvar,
	      const int&                              a_redistRad);

protected:

  void gatherBroadcast(   IntVectSet& a_set);
  void makeCedFineSets(   IntVectSet& a_cfivs, const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar);
  void makeCoarSets(const IntVectSet& a_cfivs, const LayoutData<Vector<LayoutIndex> >& a_neighborsFine);

  // Debugging functions
  void printCoarSets();
  void printCedFineSets();


};

#endif
