/*!
  @file   EBFastFineToCoarRedist.H
  @brief  New EBFineToCoarRedist class with faster defines
  @author Robert Marskar
*/

#ifndef _EBFastFineToCoarRedist_
#define _EBFastFineToCoarRedist_

#include <EBFineToCoarRedist.H>

/*!
  @brief EBFineToCoarRedist with faster define function
*/
class EBFastFineToCoarRedist : public EBFineToCoarRedist {
public:

  EBFastFineToCoarRedist();
  ~EBFastFineToCoarRedist();

  void new_define(const EBLevelGrid&                      a_eblgFine,
		  const EBLevelGrid&                      a_eblgCoar,
		  const LayoutData<Vector<LayoutIndex> >& a_neighborsFine,
		  const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar,
		  const int&                              a_nref,
		  const int&                              a_nvar,
		  const int&                              a_redistRad);

  void define(const EBLevelGrid&                      a_eblgFine,
	      const EBLevelGrid&                      a_eblgCoar,
	      const LayoutData<Vector<LayoutIndex> >& a_neighborsFine,
	      const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar,
	      const int&                              a_nref,
	      const int&                              a_nvar,
	      const int&                              a_redistRad);

protected:

  LevelData<BaseFab<bool> >* m_fineShell;    // Transitory memory
  LevelData<BaseFab<bool> >* m_refCoarShell; // Transitory memory

  void makeEmptyFineShell();
  void makeEmptyRefCoarShell();
  
  void makeFineSets(const LayoutData<Vector<LayoutIndex> >& a_neighborsFine);
  void makeCoarSets(const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar);
};

#endif
