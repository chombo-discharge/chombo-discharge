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
class ebcoarseaverage : public EBCoarseAverage {
public:

  /*!
     @brief Default constructor but must call defined afterwards. 
  */
  ebcoarseaverage();

  /*!
    @brief Copy constructor not allowed
  */
  ebcoarseaverage(const ebcoarseaverage& ebcin) = delete;

  /*!
     @brief Defining constructor. 
  */
  ebcoarseaverage(const DisjointBoxLayout& dblFine,
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
  ebcoarseaverage(const EBLevelGrid& a_eblgFine,
		  const EBLevelGrid& a_eblgCoar,
		  const EBLevelGrid& a_eblgCoFi,
		  const int& nref,
		  const int& nvar);

  /*!
    @brief Destructor
  */
  ~ebcoarseaverage();



  /*!
    @brief Assignement not allowed.
  */
  void operator=(const ebcoarseaverage& fabin) = delete;


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
