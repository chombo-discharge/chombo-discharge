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
   at fine level of refinement. Just like Chombo's EBCoarseAverage, with an extra
   function. 
*/
class EbCoarAve : public EBCoarseAverage {
public:

  /*!
     @brief Default constructor but must call defined afterwards. 
  */
  EbCoarAve();

  /*!
    @brief Copy constructor not allowed
  */
  EbCoarAve(const EbCoarAve& ebcin) = delete;

  /*!
     @brief Defining constructor. 
  */
  EbCoarAve(const DisjointBoxLayout& dblFine,
		  const DisjointBoxLayout& dblCoar,
		  const EBISLayout& ebislFine,
		  const EBISLayout& ebislCoar,
		  const ProblemDomain& domainCoar,
		  const int& nref,
		  const int& nvar,
		  const EBIndexSpace* ebisPtr);

  /*!
    @brief Defining constructor
  */
  EbCoarAve(const EBLevelGrid& a_eblgFine,
		  const EBLevelGrid& a_eblgCoar,
		  const EBLevelGrid& a_eblgCoFi,
		  const int& nref,
		  const int& nvar);

  /*!
    @brief Destructor
  */
  ~EbCoarAve();



  /*!
    @brief Assignement not allowed.
  */
  void operator=(const EbCoarAve& fabin) = delete;


  /*!
    @brief Conservatively average data defined on irregular cells
  */
  void conservativeAverage(LevelData<BaseIVFAB<Real> >&       a_coarData,
			   const LevelData<BaseIVFAB<Real> >& a_fineData,
			   const Interval&                    a_variables);

protected:

  /*!
    @brief Conservatively average data defined on irregular cells
  */
  void conservativeAverageFAB(BaseIVFAB<Real>&       a_coar,
			      const BaseIVFAB<Real>& a_fine,
			      const DataIndex&       a_datInd,
			      const Interval&        a_variables) const;



};

#include <CD_NamespaceFooter.H>

#endif
