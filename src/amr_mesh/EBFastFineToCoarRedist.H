/*!
  @file   EBFastFineToCoarRedist.H
  @brief  New EBFineToCoarRedist class with faster defines
  @author Robert Marskar
*/

#ifndef _EBFASTFINETOCOARREDIST_H_
#define _EBFASTFINETOCOARREDIST_H_

#include <EBFineToCoarRedist.H>

/*!
  @brief EBFineToCoarRedist with faster define function
*/
class EBFastFineToCoarRedist : public EBFineToCoarRedist {
public:

  EBFastFineToCoarRedist();
  ~EBFastFineToCoarRedist();

  void define(const EBLevelGrid&                      a_eblgFine,
	      const EBLevelGrid&                      a_eblgCoar,
	      const LayoutData<Vector<LayoutIndex> >& a_neighborsFine,
	      const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar,
	      const int&                              a_nref,
	      const int&                              a_nvar,
	      const int&                              a_redistRad);

protected:

  void makeFineMask(LevelData<BaseFab<bool> >& a_fineShellMask, const LayoutData<Vector<LayoutIndex> >& a_neighborsFine);
  void makeFineSets(const LevelData<BaseFab<bool> >& a_fineMask);
  void makeCoarSets(const LevelData<BaseFab<bool> >& a_refCoarMask);

  // These are debugging functions, these are never called unless you activate debugging in the .cpp file
  void gatherBroadcast(  IntVectSet& a_set);
  void gatherSetsFine(   IntVectSet& a_setsFine);   
  void gatherSetsRefCoar(IntVectSet& a_setsRefCoar);
  
  
};

#endif
