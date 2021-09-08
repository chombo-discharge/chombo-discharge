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
  @brief EBFineToCoarRedist with faster define function.
  @details This is just like Chombo's EBCoarToCoarRedist except that it is much faster to define. 
  @note To use run-time debugging and profiling, use EbFastCoarToCoarRedist.debug=true and
  EbFastCoarToCoarRedist.profile=true
*/
class EbFastCoarToCoarRedist : public EBCoarToCoarRedist {
public:

  /*!
    @brief Weak constructor. Must subsequently call the define function. 
  */
  EbFastCoarToCoarRedist();

  /*!
    @brief Destructor (does nothing).
  */
  ~EbFastCoarToCoarRedist();

  /*!
    @brief New define function for EBFastCoarToCoar
    @param[in] a_eblgFine      Fine level EBLevelGrid
    @param[in] a_eblg          Coarse level EBLevelGrid
    @param[in] a_neighborsFine Neigbor boxes on the fine grid
    @param[in] a_neighborsFine Neigbor boxes on the coarse grid
    @param[in] a_refRat        Refinement ratio between levels
    @param[in] a_nComp         Number of variables
    @param[in] a_redistRad     Redistribution radius
  */
  void fastDefine(const EBLevelGrid& a_eblgFine,
		  const EBLevelGrid& a_eblgCoar,
		  const int&         a_refRat,
		  const int&         a_nComp,
		  const int&         a_redistRad);

protected:

  /*!
    @brief Debugging function. Used when we want to ensure that this class and EbCoarToCoarRedist produce the same sets
    @param[inout] a_set Set to gather and broadcast. 
  */
  void gatherBroadcast(IntVectSet& a_set);

  /*!
    @brief Debugging function. Used when we want to ensure that this class and EbCoarToCoarRedist produce the same sets.
    @param[inout] a_coarSet Coarse set.
    @details This gathers m_setsCoar over all ranks. 
  */  
  void gatherCoarSet(IntVectSet& coarSet); 
};

#include <CD_NamespaceFooter.H>

#endif
