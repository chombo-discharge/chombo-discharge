/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbCoarAve.H
  @brief  Declaration of conservative coarsening utility 
  @author Robert Marskar
*/

#ifndef CD_EbCoarAve_H
#define CD_EbCoarAve_H

// Chombo includes
#include <EBCoarseAverage.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
   @brief Class which replaces data at coarse level of refinement with average
   at fine level of refinement. Just like Chombo's EBCoarseAverage, but with an extra
   function. 
*/
class EbCoarAve : public EBCoarseAverage {
public:

  /*!
     @brief Default constructor. Must call define afterwards. 
  */
  EbCoarAve();

  /*!
    @brief Copy constructor not allowed
  */
  EbCoarAve(const EbCoarAve& a_other) = delete;

  /*!
     @brief Defining constructor. Calls the define function. 
     @param[in] a_dblFine    Fine grids
     @param[in] a_dblCoar    Coarse grids
     @param[in] a_ebislFine  Fine EBISLayout
     @param[in] a_ebislCoar  Coarse EBISLayout
     @param[in] a_domainCoar Coarse-level domain
     @param[in] a_refRat     Refinement ratio between fine and coarse
     @param[in] a_nComp      Maximum number of variables to coarsen
     @param[in] a_ebisPtr    EBIndexSpace
  */
  EbCoarAve(const DisjointBoxLayout&       a_dblFine,
		  const DisjointBoxLayout& a_dblCoar,
		  const EBISLayout&        a_ebislFine,
		  const EBISLayout&        a_ebislCoar,
		  const ProblemDomain&     a_domainCoar,
		  const int&               a_refRat,
		  const int&               a_nComp,
		  const EBIndexSpace*      a_ebisPtr);

  /*!
    @brief Defining constructor. Calls the define function. 
    @param[in] a_eblgFine Fine grids
    @param[in] a_eblgCoar Coarse grids
    @param[in] a_eblgCoFi Coarsened fine grids
    @param[in] a_refRat   Refinement ratio between coarse and level
    @param[in] a_nComp    Maximum number of variables to coarsen
  */
  EbCoarAve(const EBLevelGrid& a_eblgFine,
	    const EBLevelGrid& a_eblgCoar,
	    const EBLevelGrid& a_eblgCoFi,
	    const int&         a_refRat,
	    const int&         a_nComp);

  /*!
    @brief Destructor (does nothing)
  */
  ~EbCoarAve();

  /*!
    @brief Assignement not allowed.
  */
  void operator=(const EbCoarAve& fabin) = delete;

  /*!
    @brief Conservatively average data defined on irregular cells
    @details a_coarData will be replaced by a conservative coarsening of a_fineData in regions where the 
    fine grid covers the coarse grid. 
    @param[inout] a_coarData  Coarse data
    @param[in]    a_fineData  Fine data
    @param[in]    a_variables Variables to coarsen
  */
  void conservativeAverage(LevelData<BaseIVFAB<Real> >&       a_coarData,
			   const LevelData<BaseIVFAB<Real> >& a_fineData,
			   const Interval&                    a_variables);

protected:

  /*!
    @brief Conservatively average data defined on irregular cells
    @details a_coarData will be replaced by a conservative coarsening of a_fineData in regions where the 
    fine grid covers the coarse grid. 
    @param[inout] a_coarData  Coarse data
    @param[in]    a_fineData  Fine data
    @param[in]    a_datInd    Grid index
    @param[in]    a_variables Variables to coarsen
  */
  void conservativeAverageFAB(BaseIVFAB<Real>&       a_coarData,
			      const BaseIVFAB<Real>& a_fineData,
			      const DataIndex&       a_datInd,
			      const Interval&        a_variables) const;



};

#include <CD_NamespaceFooter.H>

#endif
