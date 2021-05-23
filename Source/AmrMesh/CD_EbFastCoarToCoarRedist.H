/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbFastCoarToCoarRedist.H
  @brief  New EBFineToCoarRedist class with faster define function
  @author Robert Marskar
*/

#ifndef CD_EbFastCoarToCoarRedist_H
#define CD_EbFastCoarToCoarRedist_H

// Chombo includes
#include "EBCoarToCoarRedist.H"

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief EBFineToCoarRedist with faster define function
*/
class EbFastCoarToCoarRedist : public EBCoarToCoarRedist {
public:

  /*!
    @brief Weak constructor
  */
  EbFastCoarToCoarRedist();

  /*!
    @brief Destructor
  */
  ~EbFastCoarToCoarRedist();

  /*!
    @brief New define function for EBFastCoarToCoar
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

  // Debugging functions
  void gatherBroadcast(IntVectSet& a_set);
  void gatherCoarSet(IntVectSet& coarSet); // Debugging function
};

#include <CD_NamespaceFooter.H>

#endif
