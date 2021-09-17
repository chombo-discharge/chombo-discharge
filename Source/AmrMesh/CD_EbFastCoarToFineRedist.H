/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbFastCoarToFineRedist.H
  @brief  New EBCoarToFineRedist with much faster define functions
  @author Robert Marskar
*/

#ifndef CD_EbFastCoarToFineRedist_H
#define CD_EbFastCoarToFineRedist_H

// Chombo includes
#include <EBCoarToFineRedist.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief EBFineToCoarRedist with faster define function.
  @details This is just like Chombo's EBCoarToFineRedist except that it is much faster to define. 
  @note To use run-time debugging and profiling, use EbFastCoarToFineRedist.debug=true and
  EbFastCoarToFineRedist.profile=true
*/
class EbFastCoarToFineRedist : public EBCoarToFineRedist {
public:

  /*!
    @brief Weak constructor. Must subsequently call the define function. 
  */
  EbFastCoarToFineRedist();

  /*!
    @brief Destructor.
  */
  ~EbFastCoarToFineRedist();

  /*!
    @brief New define function for EBFastCoarToFine
    @param[in] a_eblgFine      Fine level EBLevelGrid
    @param[in] a_eblgCoar      Coarse level EBLevelGrid
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
    @brief Make the coarsened fine sets
  */
  void makeCedFineSets();

  /*!
    @brief Make the coarse set
  */  
  void makeCoarSets();

  /*!
    @brief Debugging function. Gathers and broadcasts a set
    @param[inout] a_set Set to gather and broadcast
  */
  void gatherBroadcast(IntVectSet& a_set);

  /*!
    @brief Debugging function. Gathers and broadcasts m_setsCedFine onto a_setsCedFine
    @param[inout] a_setsCedFine Fine set
  */  
  void gatherSetsCedFine(IntVectSet& a_setsCedFine);

  /*!
    @brief Debugging function. Gathers and broadcasts m_setsCoar onto a_setsCoar. 
    @param[inout] a_setsCoar Fine set
  */  
  void gatherSetsCoar(IntVectSet& a_setsCoar);
};

#include <CD_NamespaceFooter.H>

#endif
