/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbFastFineToCoarRedist.H
  @brief  New EBFineToCoarRedist class with faster define function.
  @author Robert Marskar
*/

#ifndef CD_EbFastFineToCoarRedist_H
#define CD_EbFastFineToCoarRedist_H

// Chombo includes
#include <EBFineToCoarRedist.H>

// Our includes
#include <CD_NamespaceHeader.H>
  
/*!
  @brief EBFineToCoarRedist with faster define function
*/
class EbFastFineToCoarRedist : public EBFineToCoarRedist {
public:

  /*!
    @brief Weak constructor.
  */
  EbFastFineToCoarRedist();

  /*!
    @brief Destructor.
  */
  ~EbFastFineToCoarRedist();

    /*!
    @brief New define function for EBFastCoarToFine
    @param[in] a_eblgFine      Fine level EBLevelGrid
    @param[in] a_eblg          Coarse level EBLevelGrid
    @param[in] a_neighborsFine Neigbor boxes on the fine grid
    @param[in] a_neighborsFine Neigbor boxes on the coarse grid
    @param[in] a_nRef          Refinement ratio between levels
    @param[in] a_nVar          Number of variables
    @param[in] a_redistRad     Redistribution radius
  */
  void define(const EBLevelGrid&                      a_eblgFine,
	      const EBLevelGrid&                      a_eblgCoar,
	      const LayoutData<Vector<LayoutIndex> >& a_neighborsFine,
	      const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar,
	      const int&                              a_nRef,
	      const int&                              a_nVar,
	      const int&                              a_redistRad);

protected:

  /*!
    @brief Make a mask which holds true for cells within the specified radius on the inside of the refinement boundary.
    @param[in] a_fineShellMask Fine grid shell mask, holding the cells. 
    @param[in] a_neighborsFine Neighbors on the fine grid level. 
  */
  void makeFineMask(LevelData<BaseFab<bool> >& a_fineShellMask, const LayoutData<Vector<LayoutIndex> >& a_neighborsFine);

  /*!
    @brief Fine sets, as per EBFineToCoarRedist.
    @param[in] a_fineMask Mask which is true on the inside of the CFIVS. 
  */
  void makeFineSets(const LevelData<BaseFab<bool> >& a_fineMask);

  /*!
    @brief Make coarse sets, as per EBFineToCoarRedist
    @param[in] a_refCoarMask Mask defined on the refined version of the coarse grid. 
    @details This 
  */
  void makeCoarSets(const LevelData<BaseFab<bool> >& a_refCoarMask);

  // These are debugging functions, these are never called unless you activate debugging in the .cpp file
  void gatherBroadcast(  IntVectSet& a_set);
  void gatherSetsFine(   IntVectSet& a_setsFine);   
  void gatherSetsRefCoar(IntVectSet& a_setsRefCoar);
};

#include <CD_NamespaceFooter.H>

#endif
