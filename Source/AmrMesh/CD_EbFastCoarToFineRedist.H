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
  @brief EBFineToCoarRedist with faster define function
*/
class EbFastCoarToFineRedist : public EBCoarToFineRedist {
public:

  /*!
    @brief Weak constructor
  */
  EbFastCoarToFineRedist();

  /*!
    @brief Destructor.
  */
  ~EbFastCoarToFineRedist();

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
    @brief Create a global view of the intersection of the irregular cells on the coarse side of the refinement boundary. 
    @param[out] a_cfivs Set
    @param[in]  a_neighborsCoar Coarse grid neighboring boxes. 
  */
  void makeCedFineSets(   IntVectSet& a_cfivs, const LayoutData<Vector<LayoutIndex> >& a_neighborsCoar);

  /*!
    @brief Make coarse sets, as per EBCoarToFineRedist
    @param[in] a_cfivs Set
    @param[in] a_neighborsFine Neighbors on fine level.
  */
  void makeCoarSets(const IntVectSet& a_cfivs, const LayoutData<Vector<LayoutIndex> >& a_neighborsFine);

  // Debugging functions
  void gatherBroadcast(  IntVectSet& a_set);
  void gatherSetsCedFine(IntVectSet& a_setsCedFine);
  void gatherSetsCoar(   IntVectSet& a_setsCoar);
};

#include <CD_NamespaceFooter.H>

#endif
